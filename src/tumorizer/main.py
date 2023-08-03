# Copyright 2022 - Barcelona Supercomputing Center
# Author: Rodrigo Martin
# BSC Dual License

import os
import sys
import argparse
import logging
import random
import numpy as np
import pysam


from common import index_alignment, replace_reads, concat_vcf_files  # noqa
from bamsurgeon_tools import run_bamsurgeon, vct_to_bamsurgeon  # noqa


def _stats_from_alignment(alignment_file, fasta_ref):
    """
    Calculate average read size from an alignment file.
    """
    f = pysam.AlignmentFile(alignment_file, reference_filename=fasta_ref)
    # Aproximate by the middle 1000 reads for each contig
    num_reads = 1000
    contigs = f.references
    read_sizes = []
    insert_sizes = []
    for contig in contigs:
        contig_len = f.get_reference_length(contig)
        for i, read in enumerate(f.fetch(contig=contig, start=contig_len // 2)):
            if i >= num_reads:
                break
            read_sizes.append(read.query_length)
            insert_sizes.append(abs(read.template_length))
    f.close()
    data_stats = dict()
    data_stats['read_size_median'] = np.median(read_sizes)
    data_stats['insert_size_median'] = np.median(insert_sizes)
    data_stats['insert_size_iqr'] = np.percentile(insert_sizes, 75) - np.percentile(insert_sizes, 25)
    return data_stats


def variate_alignment(args):
    """
    Add variants from a VCF file to an alignment file.
    """

    # Get original alignment stats
    original_stats = _stats_from_alignment(args.input_alignment, args.fasta_ref)
    # Get indel threshold based on median read size (as recommended by bamsurgeon)
    indel_threshold = int(original_stats['read_size_median'] * 0.9)
    # Convert VCF to BAMSurgeon input
    bamsurgeon_input_snv, bamsurgeon_input_indel, bamsurgeon_input_sv = vct_to_bamsurgeon(
        args.vcf_file, args.tmp_dir, args.vaf, indel_threshold)

    # Check commands to run
    run_snv = os.stat(bamsurgeon_input_snv).st_size > 0
    run_indel = os.stat(bamsurgeon_input_indel).st_size > 0
    run_sv = os.stat(bamsurgeon_input_sv).st_size > 0
    command_list = []
    if run_snv:
        command_list.append(('addsnv.py', bamsurgeon_input_snv))
    if run_indel:
        command_list.append(('addindel.py', bamsurgeon_input_indel))
    if run_sv:
        command_list.append(('addsv.py', bamsurgeon_input_sv))

    # BAMSurgeon parameters
    bamsurgeon_extra_params = {
        'vaf': str(args.vaf),
        'maxlibsize': str(int(original_stats['insert_size_median'] * 0.9)),
        'ismean': str(original_stats['insert_size_median']),
        'issd': str(original_stats['insert_size_iqr']),
        'seed': args.seed,
        'processes': args.processes,
        'max_memory': args.max_memory,
        'donor': args.donor
    }

    # Get original alignment file size in GB
    original_alignment_size = os.stat(args.input_alignment).st_size / (1024**3)

    # Run BAMSurgeon
    file_format = 'bam' if args.input_alignment.endswith('.bam') else 'cram'
    output_vcfs = []
    mutated_alignment = args.input_alignment
    for i, command_tuple in enumerate(command_list):
        command, bamsurgeon_input = command_tuple
        modified_reads_bam_file, output_vcf_file, exclude_reads_file = run_bamsurgeon(
            command, mutated_alignment, args.fasta_ref, bamsurgeon_input, args.tmp_dir, bamsurgeon_extra_params)
        output_vcfs.append(output_vcf_file)
        # Get output filename
        if i + 1 < len(command_list):
            if original_alignment_size * 1.05 > args.tmp_dir_size:
                new_mutated_alignment = args.output_alignment + command.replace('.py', '') + '_temp.' + file_format
            else:
                new_mutated_alignment = os.path.join(args.tmp_dir, command.replace('.py', '') + '_temp.' + file_format)
        else:
            new_mutated_alignment = args.output_alignment

        # Replace reads
        replace_reads(mutated_alignment, modified_reads_bam_file,
                      new_mutated_alignment, args.fasta_ref, exclude_reads_file, args.processes)
        # Index new alignment file
        index_alignment(new_mutated_alignment, args.processes)

        # Remove temporal files
        os.remove(modified_reads_bam_file)
        if exclude_reads_file:
            os.remove(exclude_reads_file)
        if mutated_alignment != args.input_alignment:
            os.remove(mutated_alignment)
            if file_format == 'cram':
                os.remove(mutated_alignment + '.crai')
            else:
                os.remove(mutated_alignment + '.bai')
        mutated_alignment = new_mutated_alignment

    # Concat VCF files
    out_vcf_path = args.output_alignment + '.vcf'
    concat_vcf_files(output_vcfs, out_vcf_path, args.fasta_ref)
    logging.info(f'Concatenated BAMSurgeon\'s VCF files to {out_vcf_path}')

    # Remove BAMSurgeon's VCF files
    [os.remove(vcf) for vcf in output_vcfs]

    # Remove bamsurgeon input files
    os.remove(bamsurgeon_input_snv)
    os.remove(bamsurgeon_input_indel)
    os.remove(bamsurgeon_input_sv)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-f', '--fasta-ref', type=str, required=True, help='Reference FASTA file.')
    parser.add_argument('-v', '--vcf-file', type=str, required=True, help='VCF file with variants.')
    parser.add_argument('-i', '--input-alignment', type=str, required=True, help='Input alignment file (SAM/BAM/CRAM).')
    parser.add_argument('-o', '--output-alignment', type=str, required=True,
                        help='Output alignment file (SAM/BAM/CRAM).')
    parser.add_argument('--vaf', type=float, help='Variant allele frequency (0.0-1.0).')
    parser.add_argument('--donor', type=str, help='Donor alignment file (for BAMSurgeon).')
    parser.add_argument('-td', '--tmp-dir', type=str, default='.', help='Directory where temporal files are stored.')
    parser.add_argument('-tds', '--tmp-dir-size', type=float, default=200,
                        help='Maximum size of temporal directory (in GB).')
    parser.add_argument('-p', '--processes', type=int, default=1, help='Maximum number of processes.')
    parser.add_argument('-mm', '--max-memory', type=int, default=96, help='Maximum memory usage (in GB).')
    parser.add_argument('-s', '--seed', type=int, default=random.randint(0, 2**32), help='Random seed.')

    args = parser.parse_args()
    logging.basicConfig(format='%(levelname)s %(asctime)s %(message)s', level='INFO')

    if not args.donor:
        args.donor = args.input_alignment

    # Set random seed
    random.seed(args.seed)

    # Mutate alignment file
    variate_alignment(args)
