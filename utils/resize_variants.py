# Copyright 2022 - Barcelona Supercomputing Center
# Author: Rodrigo Martin
# BSC Dual License

import argparse
import pandas as pd
import random
import pysam


from variant_extractor import VariantExtractor  # noqa


def _read_vcf(vcf_file):
    extractor = VariantExtractor(vcf_file, pass_only=True)
    variant_data = []
    for variant_record in extractor:
        length = variant_record.length
        type_inferred = variant_record.variant_type.name

        if length == 0 and not type_inferred == 'TRA' and not type_inferred == 'SNV' and not type_inferred == 'SGL':
            print('Warning: Skipped variant with 0 length.')
            continue

        variant_data.append([length, type_inferred, variant_record])

    print('Read {} variants from {}'.format(len(variant_data), vcf_file))
    return pd.DataFrame(variant_data, columns=['length', 'type_inferred', 'variant_obj'])


def _shorten_variants(variant_row, min_length, max_length):
    new_length = random.randint(min_length, max_length)
    variant_row['length'] = new_length
    old_variant_record = variant_row['variant_obj']
    if old_variant_record.variant_type.name == 'INS':
        if old_variant_record.alt_sv_shorthand is not None:
            new_variant_record = old_variant_record._replace(length=new_length)
        else:
            new_alt = old_variant_record.alt[:new_length]
            new_variant_record = old_variant_record._replace(length=new_length, alt=new_alt)
    else:
        new_end = old_variant_record.pos + new_length
        if old_variant_record.alt_sv_shorthand is not None:
            new_variant_record = old_variant_record._replace(length=new_length, end=new_end)
        elif old_variant_record.alt_sv_breakend is not None:
            new_alt = f'{old_variant_record.alt_sv_breakend.prefix}{old_variant_record.alt_sv_breakend.bracket}{old_variant_record.alt_sv_breakend.contig}:{old_variant_record.pos + new_length}{old_variant_record.alt_sv_breakend.bracket}{old_variant_record.alt_sv_breakend.suffix}'
            new_bracket_sv_record = old_variant_record.alt_sv_breakend._replace(pos=new_end)
            new_variant_record = old_variant_record._replace(
                length=new_length, alt=new_alt, end=new_end, alt_sv_breakend=new_bracket_sv_record)
        else:
            if old_variant_record.variant_type.name == 'DUP':
                new_alt = old_variant_record.alt[:new_length]
                new_variant_record = old_variant_record._replace(length=new_length, alt=new_alt, end=new_end)
            else:  # DEL
                new_ref = old_variant_record.ref[:new_length]
                new_variant_record = old_variant_record._replace(length=new_length, ref=new_ref, end=new_end)

    if 'SVLEN' in old_variant_record.info:
        old_variant_record.info['SVLEN'] = new_length

    new_variant_record.info['SHORTENED'] = True
    variant_row['variant_obj'] = new_variant_record
    return variant_row


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument(action='store', dest='input_vcf', help='Input VCF file')
    parser.add_argument(action='store', dest='output_vcf', help='Output VCF file')
    parser.add_argument('--target-size', type=int, nargs=2, default=[100, 600], help='Target size of the variants')
    parser.add_argument('--affected-size', type=int, nargs=2,
                        default=[2000, 999999999], help='Size of the affected variants')
    parser.add_argument('--num-modifications', type=int, default=4000, help='Number of modifications to perform')
    parser.add_argument('--seed', type=int, default=None, help='Random seed')

    args = parser.parse_args()

    if args.seed is not None:
        random.seed(args.seed)

    vcf_data = _read_vcf(args.input_vcf)

    # Filter variants by affected size
    affected_size_mask = (vcf_data['length'] >= args.affected_size[0]) & (vcf_data['length'] <= args.affected_size[1])
    affected_variants = vcf_data[affected_size_mask]

    # Get a subset of the affected variants
    affected_variants = affected_variants.sample(n=args.num_modifications, random_state=42)
    print(affected_variants)
    # Shorten the affected variants
    modified_variants = affected_variants.apply(_shorten_variants, axis=1, args=(args.target_size[0], args.target_size[1]))
    # Merge the affected variants with the original data
    vcf_data.loc[modified_variants.index] = modified_variants
    print(modified_variants)

    # Write the new VCF file
    # Extract header from the original VCF file
    with pysam.VariantFile(args.input_vcf) as vcf_file:
        header = vcf_file.header
    # Write the new VCF file
    with open(args.output_vcf, 'w') as vcf_file:
        vcf_file.write(str(header))
        for variant_row in vcf_data.itertuples():
            vcf_file.write(str(variant_row.variant_obj))
            vcf_file.write('\n')
