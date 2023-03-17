# Copyright 2022 - Barcelona Supercomputing Center
# Author: Rodrigo Martin
# BSC Dual License

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

from variant_extractor import VariantExtractor

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
    filter_array = np.full(len(df), False)
    # Set n random positions to True
    filter_array[np.random.choice(len(df), n, replace=False)] = True
    return df[filter_array]


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
        if variant_record.alt_sv_breakend:
            end_chrom = variant_record.alt_sv_breakend.contig.replace('chr', '')
            if start_chrom != end_chrom:
                end = variant_record.alt_sv_breakend.pos
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
            prefix = 'N' if variant_record.alt_sv_breakend.prefix else '' # type: ignore
            suffix = 'N' if variant_record.alt_sv_breakend.suffix else '' # type: ignore
            brackets = prefix + variant_record.alt_sv_breakend.bracket + suffix # type: ignore
        elif type_inferred == 'TRA':
            prefix = 'N' if variant_record.alt_sv_breakend.prefix else '' # type: ignore
            suffix = 'N' if variant_record.alt_sv_breakend.suffix else '' # type: ignore
            brackets = prefix + variant_record.alt_sv_breakend.bracket + variant_record.alt_sv_breakend.bracket + suffix # type: ignore

        if length == 0 and not type_inferred == 'TRA' and not type_inferred == 'SNV' and not type_inferred == 'SGL':
            print('Warning: Skipped variant with 0 length.')
            continue

        # Remove FORMAT and SAMPLE fields
        variant_record = variant_record._replace(format=[], samples=dict())

        # if 'VAF' not in variant_record.info:
        #     # if len(RANDOM_VAFS) > 0:
        #     #     vaf = np.random.choice(RANDOM_VAFS, p=RANDOM_VAFS_WEIGHTS)
        #     #     variant_record.info['VAF'] = vaf
        #     ev_reads = variant_record.info['TumourEvidenceReads']
        #     ev_reads = int(ev_reads) if isinstance(ev_reads, str) else sum([int(x) for x in ev_reads])
        #     # if 'TumourBp2ClipEvidence' in variant_record.info:
        #     #     ev_reads += (int(variant_record.info['TumourBp2ClipEvidence']) + int(variant_record.info['TumourBp1ClipEvidence'])) / 2
        #     total_reads = int(variant_record.info['TumourReads'])
        #     vaf = ev_reads / total_reads
        #     variant_record.info['VAF'] = vaf
        #     # variant_record.info['file'] = vcf_file.split('/')[-1]
        #     # if length > 100 and length < 550:  # Remove small indels
        #     #     continue

        variant_data.append([start_chrom, start, ref, alt, length, end_chrom,
                            end, brackets, type_inferred, variant_record])

    print('Read {} variants from {}'.format(len(variant_data), vcf_file))
    return pd.DataFrame(variant_data, columns=['start_chrom', 'start', 'ref', 'alt', 'length', 'end_chrom', 'end', 'brackets', 'type_inferred', 'variant_obj'])


def extract_variants(vcf_files):
    with ThreadPoolExecutor() as executor:
        variants_dfs = executor.map(_read_vcf, vcf_files)
    return pd.concat(variants_dfs, ignore_index=True)


def _maximum_non_overlapping_variants(variants, compare_values, intervals_used, padding):
    def in_used_intervals(val):
        index_start = bisect.bisect_right(intervals_used, (val[start], -1))
        in_interval = (index_start > 0 and intervals_used[index_start-1][1]+padding >= val[start]) \
            or (index_start != len(intervals_used) and intervals_used[index_start][0]-padding <= val[end])
        return in_interval

    selected_variants = np.full(len(variants), False)
    start, end = compare_values
    i = 0
    last_end = 0
    while i < len(variants):
        curr_var = variants.iloc[i]
        if curr_var[start] <= last_end + padding or in_used_intervals(curr_var):
            i += 1
        elif i == len(variants) - 1:
            selected_variants[i] = True
            i += 1
        else:
            next_var = variants.iloc[i+1]
            if curr_var[end] + padding < next_var[start] or in_used_intervals(next_var):
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


def select_variants(variants_df, indel_threshold, num_variants_list, padding):
    snv_max, indel_max, trn_max, inv_max, dup_max, del_max, ins_max = num_variants_list
    contigs = variants_df['start_chrom'].unique()
    variants_df['start_length'] = variants_df['start'] + variants_df['length']
    variants_df.sort_values(by=['start'], inplace=True)
    variants_df.reset_index(drop=True, inplace=True)
    used_intervals_by_contig = dict()
    for chrom in contigs:
        used_intervals_by_contig[chrom] = []

    selected_variants = []

    def non_overlapping_trn(variants, padding, max_variants):
        chrom_selected_variants_list = []
        # Copy used intervals for temporal use
        temp_used_intervals_by_contig = dict()
        for chrom in contigs:
            temp_used_intervals_by_contig[chrom] = used_intervals_by_contig[chrom]

        for chrom in contigs:
            chrom_variants = variants[variants['start_chrom'] == chrom]
            chrom_selected_variants_start = _maximum_non_overlapping_variants(
                chrom_variants, ['start', 'start'], temp_used_intervals_by_contig[chrom], padding)
            for end_chrom in chrom_selected_variants_start['end_chrom'].unique():
                if end_chrom == chrom:
                    temp_used_intervals = list(heapq.merge(temp_used_intervals_by_contig[chrom],
                                                           chrom_selected_variants_start[['start', 'start']].itertuples(index=False, name=None)))
                else:
                    temp_used_intervals = temp_used_intervals_by_contig[end_chrom]
                end_base_variants = chrom_selected_variants_start[chrom_selected_variants_start['end_chrom'] == end_chrom].sort_values(by=[
                                                                                                                                       'end'])
                chrom_selected_variants_end = _maximum_non_overlapping_variants(
                    end_base_variants, ['end', 'end'], temp_used_intervals, padding)
                # Update used intervals
                temp_used_intervals_by_contig[end_chrom] = list(heapq.merge(temp_used_intervals_by_contig[end_chrom],
                                                                            chrom_selected_variants_end[['end', 'end']].itertuples(index=False, name=None)))
                # Get the variants keeping order
                chrom_selected_variants = chrom_selected_variants_start[chrom_selected_variants_start.index.isin(
                    chrom_selected_variants_end.index)]

                # Update used intervals
                temp_used_intervals_by_contig[chrom] = list(heapq.merge(temp_used_intervals_by_contig[chrom],
                                                                        chrom_selected_variants[['start', 'start']].itertuples(index=False, name=None)))
                chrom_selected_variants_list.append(chrom_selected_variants)

        # Concatenate all variants and select the maximum number of variants
        curr_selected_variants = pd.concat(chrom_selected_variants_list)
        curr_selected_variants = _subsample_df(curr_selected_variants, max_variants)
        # Update used intervals
        end_sorted_variants = curr_selected_variants.sort_values(by=['end'])
        curr_selected_variants = curr_selected_variants.sort_values(by=['start'])
        for chrom in contigs:
            curr_selected_variants_chrom = curr_selected_variants[curr_selected_variants['start_chrom'] == chrom]
            end_sorted_variants_chrom = end_sorted_variants[end_sorted_variants['end_chrom'] == chrom]
            used_intervals_by_contig[chrom] = list(heapq.merge(used_intervals_by_contig[chrom],
                                                               curr_selected_variants_chrom[['start', 'start']].itertuples(index=False, name=None)))
            used_intervals_by_contig[chrom] = list(heapq.merge(used_intervals_by_contig[chrom],
                                                               end_sorted_variants_chrom[['end', 'end']].itertuples(index=False, name=None)))
        selected_variants.append(curr_selected_variants)

    def non_overlapping(variants, padding, max_variants, compare_values=['start', 'end']):
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
            used_intervals_by_contig[chrom] = list(heapq.merge(used_intervals_by_contig[chrom],
                                                               curr_selected_variants_chrom[compare_values].itertuples(index=False, name=None)))
        selected_variants.append(curr_selected_variants)

    # TRA
    trn_base_variants = variants_df[variants_df['type_inferred'] == 'TRA']
    non_overlapping_trn(trn_base_variants, padding, trn_max)

    # INV
    inv_base_variants = variants_df[variants_df['type_inferred'] == 'INV']
    non_overlapping_trn(inv_base_variants, padding, inv_max)

    # INS
    ins_base_variants = variants_df[(variants_df['type_inferred'] == 'INS') &
                                    (variants_df['length'] > indel_threshold)]
    non_overlapping(ins_base_variants, padding, ins_max, ['start', 'start_length'])
    # DUP
    dup_base_variants = variants_df[(variants_df['type_inferred'] == 'DUP')]
    non_overlapping(dup_base_variants, padding, dup_max)

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
    non_overlapping(del_base_variants, padding, del_max)
    # Indel INS
    indel_ins_base_variants = variants_df[(variants_df['type_inferred'] == 'INS') &
                                          (variants_df['length'] <= indel_threshold)]
    non_overlapping(indel_ins_base_variants, padding, indel_max // 2, ['start', 'start_length'])
    # Indel DEL
    indel_del_base_variants = variants_df[(variants_df['type_inferred'] == 'DEL') &
                                          (variants_df['length'] <= indel_threshold)]
    non_overlapping(indel_del_base_variants, padding, indel_max // 2)
    # SNV
    snv_base_variants = variants_df[variants_df['type_inferred'] == 'SNV']
    non_overlapping(snv_base_variants, padding, snv_max)

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
        for _, variant_record_obj in variants_df['variant_obj'].items():
            output_vcf.write(str(variant_record_obj)+'\n')


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--inputs', required=True, nargs='+', type=str, help='Input VCF files')
    parser.add_argument('-o', '--output', required=True, type=str, help='Output VCF file')
    parser.add_argument('-n', '--num-variants', type=int, nargs=7, default=[200, 200, 200, 200, 200, 200, 200],
                        help='Maximum number of variants to extract of each type [SNV, INDEL, TRA, INV, DUP, DEL, INS]. Default: [200, 200, 200, 200, 200, 200, 200]')
    parser.add_argument('--padding', type=int, default=550, help='Minumum padding (bp) between variants')
    parser.add_argument('--indel-threshold', type=int, default=100, help='Maximum length of an indel')

    args = parser.parse_args()

    variants_df = extract_variants(args.inputs)
    selected_variants_df = select_variants(variants_df,  args.indel_threshold, args.num_variants, args.padding)
    write_vcf(selected_variants_df, args.output, args.inputs)
