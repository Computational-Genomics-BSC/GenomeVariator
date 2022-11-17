# Copyright 2022 - Barcelona Supercomputing Center
# Author: Rodrigo Martin
# BSC AS IS License

import os
import glob
import random
import subprocess
import logging
import time
import shutil

from variant_extractor.variants import VariantType
from variant_extractor import VariantExtractor

BAMSURGEON_DIR = os.path.abspath(os.path.join(os.path.dirname(
    os.path.realpath(__file__)), '..', '..', 'dependencies', 'bamsurgeon'))
PICARD_JAR = os.path.abspath(os.path.join(os.path.dirname(
    os.path.realpath(__file__)), '..', '..', 'dependencies', 'picard.jar'))


def _generate_random_dna(length):
    '''
    Generates a random DNA sequence of the given length
    '''
    return ''.join(random.choices(['A', 'C', 'G', 'T'], weights=[0.3, 0.2, 0.2, 0.3], k=length))


def vct_to_bamsurgeon(vcf_file, output_dir, vaf, indel_threshold):
    """
    Convert a VCF file to BAMSurgeon input files.
    """

    # Open output files
    output_snv_path = os.path.join(output_dir, 'bamsurgeon_input_snv.txt')
    output_indel_path = os.path.join(output_dir, 'bamsurgeon_input_indel.txt')
    output_sv_path = os.path.join(output_dir, 'bamsurgeon_input_sv.txt')
    output_file_snv = open(output_snv_path, 'w')
    output_file_indel = open(output_indel_path, 'w')
    output_file_sv = open(output_sv_path, 'w')

    extractor = VariantExtractor(vcf_file, pass_only=True)
    # Read VCF file
    for variant_record in extractor:
        var_vaf = vaf
        if 'VAF' in variant_record.info:
            var_vaf = variant_record.info['VAF']

        if variant_record.variant_type == VariantType.SNV:
            output_file_snv.write(
                f'{variant_record.contig} {variant_record.pos} {variant_record.pos} {var_vaf} {variant_record.alt}\n')
        else:
            if variant_record.variant_type == VariantType.TRN or variant_record.variant_type == VariantType.INV:
                # Convert INV to TRN since most of them are not complete
                alt_contig = variant_record.alt_sv_bracket.contig
                alt_pos = variant_record.alt_sv_bracket.pos
                # Calculate strand notation
                if variant_record.alt_sv_bracket.bracket == '[' and variant_record.alt_sv_bracket.prefix:
                    strand_notation = '++'
                elif variant_record.alt_sv_bracket.bracket == ']' and variant_record.alt_sv_bracket.prefix:
                    strand_notation = '+-'
                elif variant_record.alt_sv_bracket.bracket == '[' and not variant_record.alt_sv_bracket.prefix:
                    strand_notation = '-+'
                else:
                    strand_notation = '--'

                op = f'TRN {alt_contig} {alt_pos} {alt_pos} {strand_notation} {var_vaf}'
                output_file_sv.write(
                    f'{variant_record.contig} {variant_record.pos} {variant_record.pos} {op}\n')
            elif variant_record.variant_type == VariantType.DUP:
                op = f'DUP 1 {var_vaf}'
                output_file_sv.write(
                    f'{variant_record.contig} {variant_record.pos} {variant_record.end} {op}\n')
            elif variant_record.variant_type == VariantType.DEL:
                # Check if INDEL
                if variant_record.end - variant_record.pos < indel_threshold:
                    output_file_indel.write(
                        f'{variant_record.contig} {variant_record.pos} {variant_record.end} {var_vaf} DEL\n')
                else:
                    op = f'DEL {var_vaf}'
                    output_file_sv.write(
                        f'{variant_record.contig} {variant_record.pos} {variant_record.end} {op}\n')
            elif variant_record.variant_type == VariantType.INS:
                if variant_record.alt_sv_shorthand:
                    insert_length = variant_record.length
                    dna_sequence = _generate_random_dna(insert_length)
                else:
                    dna_sequence = variant_record.alt[1:]
                    insert_length = len(dna_sequence)
                # Check if INDEL
                if insert_length < indel_threshold:
                    output_file_indel.write(
                        f'{variant_record.contig} {variant_record.pos} {variant_record.pos+1} {var_vaf} INS {dna_sequence}\n')
                else:
                    # Cannot set vaf for insertions
                    op = f'INS {dna_sequence}'
                    output_file_sv.write(
                        f'{variant_record.contig} {variant_record.pos} {variant_record.pos} {op}\n')

    logging.info(
        f'Converted {vcf_file} to BAMSurgeon input files {output_snv_path} {output_indel_path} {output_sv_path}.')
    output_file_sv.close()
    output_file_snv.close()
    output_file_indel.close()

    return output_snv_path, output_indel_path, output_sv_path


def run_bamsurgeon(bamsurgeon_command, original_alignment, fasta_ref, bamsurgeon_input, output_dir, extra_params):
    """
    Run BAMSurgeon.
    """
    # Add BAMSurgeon to PYTHONPATH
    env_dict = os.environ.copy()
    if 'PYTHONPATH' in env_dict:
        env_dict['PYTHONPATH'] += ':' + BAMSURGEON_DIR
    else:
        env_dict['PYTHONPATH'] = BAMSURGEON_DIR

    # Change all paths to absolute
    bamsurgeon_input = os.path.abspath(bamsurgeon_input)
    original_alignment = os.path.abspath(original_alignment)
    fasta_ref = os.path.abspath(fasta_ref)
    bamsurgeon_command_path = os.path.abspath(os.path.join(BAMSURGEON_DIR, 'bin', bamsurgeon_command))
    output_dir = os.path.abspath(output_dir)
    output_alignment_name = f'{bamsurgeon_command}_{int(time.time())}_dummy_.bam'
    output_alignment_path = os.path.abspath(os.path.join(output_dir, output_alignment_name))
    donor = os.path.abspath(extra_params['donor'])

    seed = extra_params['seed']

    # Calculate number of processes
    max_bamsurgeon_processes = extra_params['max_memory'] // 5  # 4GB per process, more just in case
    bamsurgeon_processes = min(max_bamsurgeon_processes, extra_params['processes'])
    bamsurgeon_processes = max(bamsurgeon_processes, 1)
    aligner_threads = extra_params['processes'] // bamsurgeon_processes
    aligner_threads = max(aligner_threads, 1)

    vaf = extra_params['vaf']
    extra_args = []
    if bamsurgeon_command == 'addsv.py':
        extra_args = ['--allowN', '--maxdfrac', '0.2', '-l', extra_params['maxlibsize'], '--ismean', extra_params['ismean'],
                      '--issd', extra_params['issd'], '--donorbam', donor, '--require_exact']
    else:
        extra_args = ['--force', '--insane', '--picardjar', PICARD_JAR, '--minmutreads', '0']

    # Run BAMSurgeon
    bamsurgeon_args = ['python3', '-O', bamsurgeon_command_path, '-v', bamsurgeon_input, '--mindepth', '8', '--maxdepth', '600',
                       '-s', vaf, '-r', fasta_ref, '-f', original_alignment, '-o', output_alignment_path, '-p',
                       str(bamsurgeon_processes), '--skipmerge', '--aligner', 'mem', '--alignerthreads', str(aligner_threads),
                       '--tmpdir', output_dir, '--seed', str(seed)] + extra_args

    logging.info(f'Running BAMSurgeon: { " ".join(bamsurgeon_args)}')
    subprocess.run(bamsurgeon_args, check=True, env=env_dict, cwd=output_dir)
    command_raw = bamsurgeon_command.replace('.py', '')
    # Remove logs
    shutil.rmtree(os.path.join(output_dir, f'{command_raw}_logs_{output_alignment_name}'))

    # Get output files
    # Get alignment file
    list_of_files = glob.glob(os.path.join(output_dir, command_raw+'.*.bam'))
    modified_reads_bam_file = max(list_of_files, key=os.path.getctime)
    # Get VCF file
    list_of_files = glob.glob(os.path.join(output_dir, f'*.{command_raw}.*.vcf'))
    output_vcf_file = max(list_of_files, key=os.path.getctime)
    # Get excluded reads file
    list_of_files = glob.glob(os.path.join(output_dir, command_raw+'.exclude.final.*.txt'))
    if list_of_files:
        excluded_reads = max(list_of_files, key=os.path.getctime)
    else:
        excluded_reads = None
    return modified_reads_bam_file, output_vcf_file, excluded_reads
