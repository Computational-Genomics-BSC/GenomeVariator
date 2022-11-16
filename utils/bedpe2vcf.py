# Copyright 2022 - Barcelona Supercomputing Center
# Author: Rodrigo Martin Posada
# BSC AS IS License

import argparse
import gzip

VCF_HEADER = ('##fileformat=VCFv4.2\n'
              '#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO\n')


def parse_bedpe_line(line):
    fields = line.split('\t')
    chrom_1 = fields[0]
    start_1 = int(fields[1]) + 1
    chrom_2 = fields[3]
    start_2 = int(fields[4]) + 1
    id_ = fields[6]
    ref = 'N'
    strand_1 = fields[8]
    strand_2 = fields[9]
    bracket = '[' if strand_2 == '-' else ']'
    prefix = ref if strand_1 == '+' else ''
    suffix = ref if strand_1 == '-' else ''
    alt = f'{prefix}{bracket}{chrom_2}:{start_2}{bracket}{suffix}'
    vcf_str = f'{chrom_1}\t{start_1}\t{id_}\t{ref}\t{alt}\t.\tPASS\t.'
    return vcf_str


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument(action='store', dest='input', help='Input BEDPE file')
    parser.add_argument(action='store', dest='output', help='Output VCF file')

    args = parser.parse_args()

    f_in = open(args.input, 'r') if not args.input.endswith('.gz') else gzip.open(args.input, 'rt')
    # Skip header
    try:
        _ = f_in.readline()
    except:
        raise Exception(f'Error reading input file {args.input}')
    
    f_out = open(args.output, 'w')
    f_out.write(VCF_HEADER)
    for line in f_in.readlines():
        vcf_line = parse_bedpe_line(line)
        f_out.write(vcf_line)
        f_out.write('\n')
