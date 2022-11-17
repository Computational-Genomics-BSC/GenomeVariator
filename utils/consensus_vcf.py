# Copyright 2022 - Barcelona Supercomputing Center
# Author: Rodrigo Martin
# BSC AS IS License

import os
import sys
import argparse
import heapq
import bisect
import re
from concurrent.futures import ThreadPoolExecutor
import pandas as pd
import numpy as np
import pysam

# Add variant_extractor to PYTHONPATH
VARIANT_EXTRACTOR_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__),
                                                     '..', 'dependencies', 'variant-extractor', 'src'))
sys.path.insert(0, VARIANT_EXTRACTOR_DIR)

from variant_extractor import VariantExtractor  # noqa

SEED = 1
COMPLETE_INVERSION_THRESHOLD = 50
RANDOM_VAFS = [0.05, 0.1, 0.2, 0.3, 0.4, 0.5]
RANDOM_VAFS_WEIGHTS = [0.05, 0.2, 0.3, 0.25, 0.1, 0.1]

np.random.seed(SEED)


def _contig_to_int(contig):
    contig = contig.lower().replace('chr', '')
    if contig.isdigit():
        return int(contig)
    else:
        return 22 + ord(contig[0])

def _is_sorted(l):
    return all(l[i] <= l[i+1] for i in range(len(l) - 1))


def _subsample_df(df, n):
    if n > len(df):
        return df
    return df.sample(n, random_state=SEED)


def _read_vcf(vcf_file):
    extractor = VariantExtractor(vcf_file, pass_only=True)
    variant_data = []
    for variant_record in extractor:
        start_chrom = variant_record.contig.replace('chr', '')
        start = variant_record.pos
        ref = variant_record.ref
        alt = variant_record.alt
        length = variant_record.length
        end = variant_record.end
        if variant_record.alt_sv_bracket:
            end_chrom = variant_record.alt_sv_bracket.contig.replace('chr', '')
            if start_chrom != end_chrom:
                end = variant_record.alt_sv_bracket.pos
        else:
            end_chrom = start_chrom
        type_inferred = variant_record.variant_type.name

        # Brackets
        brackets = ''
        if type_inferred == 'DEL':
            brackets = 'N['
        elif type_inferred == 'DUP':
            brackets = ']N'
        elif type_inferred == 'INV':
            prefix = 'N' if variant_record.alt_sv_bracket.prefix else ''
            suffix = 'N' if variant_record.alt_sv_bracket.suffix else ''
            brackets = prefix + variant_record.alt_sv_bracket.bracket + suffix
        elif type_inferred == 'TRN':
            prefix = 'N' if variant_record.alt_sv_bracket.prefix else ''
            suffix = 'N' if variant_record.alt_sv_bracket.suffix else ''
            brackets = prefix + variant_record.alt_sv_bracket.bracket + variant_record.alt_sv_bracket.bracket + suffix

        if length == 0 and not type_inferred == 'TRN' and not type_inferred == 'SNV' and not type_inferred == 'SGL':
            print('Warning: Skipped variant with 0 length.')
            continue

        # Remove FORMAT and SAMPLE fields
        variant_record = variant_record._replace(format=[], samples=dict())

        if 'VAF' not in variant_record.info:
            # if len(RANDOM_VAFS) > 0:
            #     vaf = np.random.choice(RANDOM_VAFS, p=RANDOM_VAFS_WEIGHTS)
            #     variant_record.info['VAF'] = vaf
            ev_reads = variant_record.info['TumourEvidenceReads']
            ev_reads = int(ev_reads) if isinstance(ev_reads, str) else sum([int(x) for x in ev_reads])
            # if 'TumourBp2ClipEvidence' in variant_record.info:
            #     ev_reads += (int(variant_record.info['TumourBp2ClipEvidence']) + int(variant_record.info['TumourBp1ClipEvidence'])) / 2
            total_reads = int(variant_record.info['TumourReads'])
            vaf = ev_reads / total_reads
            variant_record.info['VAF'] = vaf
            variant_record.info['file'] = vcf_file.split('/')[-1]
            if length > 100 and length < 550: # Remove small indels
                continue

        variant_data.append([start_chrom, start, ref, alt, length, end_chrom,
                            end, brackets, type_inferred, variant_record])

    print('Read {} variants from {}'.format(len(variant_data), vcf_file))
    return pd.DataFrame(variant_data, columns=['start_chrom', 'start', 'ref', 'alt', 'length', 'end_chrom', 'end', 'brackets', 'type_inferred', 'variant_obj'])


def extract_variants(vcf_files):
    with ThreadPoolExecutor() as executor:
        variants_dfs = executor.map(_read_vcf, vcf_files)
    return pd.concat(variants_dfs, ignore_index=True)


def _maximum_non_overlapping_variants(variants, compare_values, intervals_used, padding):
    def check_in_used_intervals(val):
        index_start = bisect.bisect_right(intervals_used['start'], val[start]-padding)
        index_end = bisect.bisect_right(intervals_used['start'], val[end]+padding)
        index_end_end = bisect.bisect_right(intervals_used['end'], val[end]+padding)
        return index_start != len(intervals_used['start']) and (index_start != index_end or index_start != index_end_end)

    selected_variants = [False] * len(variants)
    start, end = compare_values
    i = 0
    last_end = 0
    while i < len(variants):
        curr_var = variants.iloc[i]
        if curr_var[start] <= last_end + padding or check_in_used_intervals(curr_var):
            i += 1
        elif i == len(variants) - 1:
            selected_variants[i] = True
            i += 1
        else:
            next_var = variants.iloc[i+1]
            if curr_var[end] + padding < next_var[start] or check_in_used_intervals(next_var):
                selected_variants[i] = True
                last_end = curr_var[end]
                i += 1
            else:
                # Overlapping variants
                if curr_var[end] <= next_var[end]:
                    selected_variants[i] = True
                    last_end = curr_var[end]
                else:
                    selected_variants[i+1] = True
                    last_end = next_var[end]
                i += 2
    return variants.loc[selected_variants]


def _retrieve_complete_inversions(selected_inversions, all_inversions):
    inversions = []
    for idx, inversion in selected_inversions.iterrows():
        inversions.append(inversion)
        posible_inversions = []
        if idx - 1 in all_inversions.index:
            posible_inversions.append(all_inversions.loc[idx - 1])
        if idx + 1 in all_inversions.index:
            posible_inversions.append(all_inversions.loc[idx + 1])
        for close_inversion in posible_inversions:
            if np.abs(inversion['start'] - close_inversion['start']) < COMPLETE_INVERSION_THRESHOLD and \
                    np.abs(inversion['end'] - close_inversion['end']) < COMPLETE_INVERSION_THRESHOLD and \
                    inversion['brackets'] != close_inversion['brackets']:
                inversions.append(close_inversion)
                # Must have the same VAF
                # close_inversion['variant_obj'].info['VAF'] = inversion['variant_obj'].info['VAF']
                break

    return pd.DataFrame.from_records(inversions, columns=selected_inversions.columns)


def select_variants(variants_df, indel_threshold, num_variants, padding):
    contigs = variants_df['start_chrom'].unique()
    variants_df['start_length'] = variants_df['start'] + variants_df['length']
    variants_df.sort_values(by=['start'], inplace=True)
    variants_df.reset_index(drop=True, inplace=True)
    used_intervals_by_contig = dict()
    for chrom in contigs:
        used_intervals_by_contig[chrom] = {'start': [], 'end': []}

    selected_variants = []

    # TRN and INV
    for chrom in contigs:
        # Start chrom TRN
        curr_variants = variants_df[variants_df['start_chrom'] == chrom]
        trn_base_variants = curr_variants[(curr_variants['type_inferred'] == 'TRN') |
                                          (curr_variants['type_inferred'] == 'INV')]
        trn_start_variants = _maximum_non_overlapping_variants(trn_base_variants, ['start', 'start'],
                                                               used_intervals_by_contig[chrom], padding)
        trn_selected_variants_list = []
        # End chrom TRN
        for end_chrom in trn_start_variants['end_chrom'].unique():
            trn_end_base_variants = trn_start_variants[trn_start_variants['end_chrom'] == end_chrom] \
                .sort_values(by=['end'])
            trn_end_variants = _maximum_non_overlapping_variants(trn_end_base_variants, ['end', 'end'],
                                                                 used_intervals_by_contig[end_chrom], padding)
            trn_selected_variants_list.append(trn_end_variants)
            used_intervals_by_contig[end_chrom]['start'] = \
                list(heapq.merge(used_intervals_by_contig[end_chrom]['start'], trn_end_variants['end']))
            used_intervals_by_contig[end_chrom]['end'] = \
                list(heapq.merge(used_intervals_by_contig[end_chrom]['end'], trn_end_variants['end']))
        if len(trn_selected_variants_list) == 0:
            continue
        trn_selected_variants = pd.concat(trn_selected_variants_list)
        trn_selected_subsampled = _subsample_df(
            trn_selected_variants[trn_selected_variants['type_inferred'] == 'TRN'], num_variants)
        inv_selected_subsampled = _subsample_df(
            trn_selected_variants[trn_selected_variants['type_inferred'] == 'INV'], num_variants)
        trn_selected_variants = pd.concat([trn_selected_subsampled, inv_selected_subsampled])
        trn_selected_variants.sort_values(by=['start'], inplace=True)
        used_intervals_by_contig[chrom]['start'] = \
            list(heapq.merge(used_intervals_by_contig[chrom]['start'], trn_selected_variants['start']))
        used_intervals_by_contig[chrom]['end'] = \
            list(heapq.merge(used_intervals_by_contig[chrom]['end'], trn_selected_variants['start']))

        complete_inversions = _retrieve_complete_inversions(
            trn_selected_variants[trn_selected_variants['type_inferred'] == 'INV'], curr_variants[curr_variants['type_inferred'] == 'INV'])
        trn_selected_variants = pd.concat(
            [trn_selected_variants[trn_selected_variants['type_inferred'] == 'TRN'], complete_inversions])

        selected_variants.append(trn_selected_variants)

    def non_overlapping(variants, padding, compare_values=['start', 'end'], max_variants=num_variants):
        chrom_selected_variants_list = []
        for chrom in contigs:
            chrom_variants = variants[variants['start_chrom'] == chrom]
            chrom_selected_variants = _maximum_non_overlapping_variants(
                chrom_variants, compare_values, used_intervals_by_contig[chrom], padding)
            chrom_selected_variants_list.append(chrom_selected_variants)
        # Concatenate all variants and select the maximum number of variants
        curr_selected_variants = pd.concat(chrom_selected_variants_list)
        curr_selected_variants = _subsample_df(curr_selected_variants, max_variants)
        # Update used intervals
        for chrom in contigs:
            curr_selected_variants_chrom = curr_selected_variants[curr_selected_variants['start_chrom'] == chrom]
            used_intervals_by_contig[chrom]['start'] = \
                list(heapq.merge(used_intervals_by_contig[chrom]['start'], curr_selected_variants_chrom[compare_values[0]]))
            used_intervals_by_contig[chrom]['end'] = \
                list(heapq.merge(used_intervals_by_contig[chrom]['end'], curr_selected_variants_chrom[compare_values[1]]))
        selected_variants.append(curr_selected_variants)

    # INS
    ins_base_variants = variants_df[(variants_df['type_inferred'] == 'INS') &
                                    (variants_df['length'] > indel_threshold)]
    non_overlapping(ins_base_variants, padding, ['start', 'start_length'])
    # DUP
    dup_base_variants = variants_df[(variants_df['type_inferred'] == 'DUP')]
    non_overlapping(dup_base_variants, padding)
    # # DEL (> 500)
    # del_base_variants = variants_df[(variants_df['type_inferred'] == 'DEL') &
    #                                 (variants_df['length'] >= 500)]
    # non_overlapping(del_base_variants, padding)
    # # DEL (< 500)
    # del_base_variants = variants_df[(variants_df['type_inferred'] == 'DEL') &
    #                                 (variants_df['length'] > indel_threshold) &
    #                                 (variants_df['length'] < 500)]
    # non_overlapping(del_base_variants, padding, max_variants=40)
    # DEL
    del_base_variants = variants_df[(variants_df['type_inferred'] == 'DEL') &
                                    (variants_df['length'] > indel_threshold)]
    non_overlapping(del_base_variants, padding)
    # Indel INS
    indel_ins_base_variants = variants_df[(variants_df['type_inferred'] == 'INS') &
                                          (variants_df['length'] <= indel_threshold)]
    non_overlapping(indel_ins_base_variants, padding, ['start', 'start_length'], max_variants=num_variants // 2)
    # Indel DEL
    indel_del_base_variants = variants_df[(variants_df['type_inferred'] == 'DEL') &
                                          (variants_df['length'] <= indel_threshold)]
    non_overlapping(indel_del_base_variants, padding, max_variants=num_variants // 2)
    # SNV
    snv_base_variants = variants_df[variants_df['type_inferred'] == 'SNV']
    non_overlapping(snv_base_variants, padding)

    all_variants_df = pd.concat(selected_variants)
    return all_variants_df


def _extract_header(vcf_files):
    with open(vcf_files[0], 'r') as f:
        with pysam.VariantFile(f) as input_file:
            main_header = input_file.header
    for vcf_file in vcf_files[1:]:
        with open(vcf_file, 'r') as f:
            with pysam.VariantFile(f) as input_file:
                main_header.merge(input_file.header)
    return main_header


def write_vcf(variants_df, output_file, template_vcfs):
    # Extract header
    header = _extract_header(template_vcfs)
    header.add_line('##INFO=<ID=VAF,Number=1,Type=Float,Description="VAF">')
    header_str = re.sub(r'\tFORMAT.*', '', str(header))  # Remove FORMAT and SAMPLE fields
    # Sort variants
    variants_df['start_chrom'] = variants_df['start_chrom'].map(_contig_to_int)
    variants_df.sort_values(by=['start_chrom', 'start'], inplace=True)
    with open(output_file, 'w') as output_vcf:
        output_vcf.write(header_str)
        for _, variant_record_obj in variants_df['variant_obj'].iteritems():
            output_vcf.write(str(variant_record_obj)+'\n')


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--inputs', required=True, nargs='+', type=str, help='Input VCF files')
    parser.add_argument('-o', '--output', required=True, type=str, help='Output VCF file')
    parser.add_argument('-n', '--num-variants', type=int, default=200,
                        help='Maximum number of variants to extract of each type')
    parser.add_argument('--padding', type=int, default=550, help='Minumum padding (bp) between variants')
    parser.add_argument('--indel-threshold', type=int, default=100, help='Maximum length of an indel')

    args = parser.parse_args()

    variants_df = extract_variants(args.inputs)
    selected_variants_df = select_variants(variants_df,  args.indel_threshold, args.num_variants, args.padding)
    write_vcf(selected_variants_df, args.output, args.inputs)
