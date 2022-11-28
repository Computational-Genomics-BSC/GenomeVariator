# Copyright 2022 - Barcelona Supercomputing Center
# Author: Rodrigo Martin
# BSC Dual License

from collections import defaultdict
import logging
import subprocess
import pysam

from variant_extractor import VariantExtractor


def _extract_header(vcf_files):
    vcf_file = vcf_files[0]
    with open(vcf_file) as f:
        with pysam.VariantFile(f) as vcf:
            main_header = vcf.header
    for vcf_file in vcf_files[1:]:
        with open(vcf_file) as f:
            with pysam.VariantFile(f) as vcf:
                main_header.merge(vcf.header)
    return main_header


def _chr_to_num(chr_name):
    return int(chr_name) if chr_name.isdigit() else 22 + ord(chr_name[0])


def concat_vcf_files(vcf_list, output_vcf, fasta_ref):
    header = _extract_header(vcf_list)
    all_variants = []
    for vcf_file in vcf_list:
        extractor = VariantExtractor(vcf_file, fasta_ref=fasta_ref)
        for variant in extractor:
            all_variants.append(variant)

    all_variants.sort(key=lambda x: (_chr_to_num(x.contig), x.pos))
    with open(output_vcf, 'w') as f:
        f.write(str(header))
        for variant in all_variants:
            f.write(str(variant))
            f.write('\n')


def index_alignment(alignment_path, processes):
    """
    Index an alignment file.
    """
    index_args = ['samtools', 'index', '-@', str(processes), alignment_path]
    logging.info(f'Indexing with: {" ".join(index_args)}')
    return subprocess.run(index_args, check=True)


def _get_read_mode(filename):
    return 'rb' if filename.endswith('.bam') else 'rc' if filename.endswith('.cram') else 'r'


def _get_write_mode(filename):
    return 'wb' if filename.endswith('.bam') else 'wc' if filename.endswith('.cram') else 'w'


def replace_reads(original_file, replace_file, output_file, exclude_file=None, threads=1):
    # Open the original file
    original_alignment = pysam.AlignmentFile(original_file, _get_read_mode(original_file), threads=threads)

    # Get excluded reads
    exclude_file_reads = set()
    if exclude_file:
        with open(exclude_file, 'r') as exclude_file_handle:
            exclude_file_reads = set([line.strip() for line in exclude_file_handle])

    exclude_reads_count = len(exclude_file_reads)
    added_reads_count = 0
    non_modified_reads_count = 0

    # Create replace reads dictionary with contig as keys and list of reads as values
    replace_reads = dict()
    for contig in original_alignment.references:
        replace_reads[contig] = []
    replace_reads['unmapped'] = []

    # Setup header and SM tag
    sm_name = 'TUMOR'
    header = original_alignment.header
    new_header = header.to_dict()
    rg_name = new_header['RG'][0]['ID'] if 'RG' in new_header else 'VARIATED_GENOME'
    new_header['RG'] = [{'ID': rg_name, 'SM': sm_name}]

    exclude_reads = exclude_file_reads.copy()
    already_used_read_flags = defaultdict(set)
    # Get reads to replace
    replace_alignment = pysam.AlignmentFile(replace_file, _get_read_mode(replace_file), threads=threads)
    for read in replace_alignment.fetch(until_eof=True):
        if read.query_name in exclude_file_reads:
            continue
        # Avoid duplicate reads
        read_flag = read.flag & 0xfc0
        read_flag_set = already_used_read_flags[read.query_name]
        if read_flag in read_flag_set:
            continue
        read_flag_set.add(read_flag)
        read.set_tag('RG', rg_name)
        # Add read to replace list in its contig
        if read.reference_name is not None:
            replace_reads[read.reference_name].append(read)
        else:
            replace_reads['unmapped'].append(read)
        # Add also read to exclude list
        added_reads_count += 1
        exclude_reads.add(read.query_name)
    replace_alignment.close()
    del exclude_file_reads
    del already_used_read_flags
    # Sort reads to replace for better performance
    for contig in replace_reads:
        replace_reads[contig] = sorted(replace_reads[contig], key=lambda x: x.reference_start)

    # Open output file
    output_alignment = pysam.AlignmentFile(output_file, _get_write_mode(output_file), header=new_header, threads=32)

    contig = None
    replace_index = 0
    # Iterate over original file
    for read in original_alignment.fetch(until_eof=True):
        if read.reference_name != contig:
            # Write remaining reads
            while contig in replace_reads and replace_index < len(replace_reads[contig]):
                output_alignment.write(replace_reads[contig][replace_index])
                replace_index += 1
            contig = read.reference_name
            if contig is None:  # Unmapped reads
                non_modified_reads_count += 1
                output_alignment.write(read)
                continue
            # Index and position of the next read to replace
            replace_index = 0
            replace_pos = replace_reads[contig][replace_index].reference_start \
                if len(replace_reads[contig]) > 0 \
                else original_alignment.get_reference_length(contig)
        # Ignore secondary and supplementary reads
        if read.query_name not in exclude_reads:
            while read.reference_start > replace_pos:
                output_alignment.write(replace_reads[contig][replace_index])
                replace_index += 1
                if replace_index < len(replace_reads[contig]):
                    replace_pos = replace_reads[contig][replace_index].reference_start
                else:
                    replace_pos = original_alignment.get_reference_length(contig)
            # Write original read
            non_modified_reads_count += 1
            output_alignment.write(read)

    # Write unmapped reads
    for read in replace_reads['unmapped']:
        output_alignment.write(read)

    original_alignment.close()
    output_alignment.close()

    logging.info(f'Replaced reads in {original_file} with {replace_file} and saved in {output_file}:\n'
                 f'\tRemoved reads: {exclude_reads_count}\n'
                 f'\tAdded reads: {added_reads_count}\n'
                 f'\tNon modified reads: {non_modified_reads_count}\n')
