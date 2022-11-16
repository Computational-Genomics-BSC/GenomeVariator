# Copyright 2022 - Barcelona Supercomputing Center
# Author: Rodrigo Martin Posada
# BSC AS IS License

import os
import sys
import argparse
import bisect
import hashlib
import random
import pysam

# Add variant_extractor to PYTHONPATH
VARIANT_EXTRACTOR_DIR = os.path.abspath(os.path.join(
    os.path.dirname(__file__), '..', 'dependencies', 'variant-extractor', 'src'))
sys.path.insert(0, VARIANT_EXTRACTOR_DIR)
from variant_extractor import VariantExtractor  # noqa
from variant_extractor.variants import VariantType  # noqa


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', required=True, type=str, help='Input alignment file')
    parser.add_argument('-ic', '--input-coverage', required=True, type=int, help='Input file coverage')
    parser.add_argument('-o', '--output', required=True, type=str, help='Output prefix')
    parser.add_argument('-oc', '--output-coverages', required=True, nargs='+', type=int, help='Output coverages')
    parser.add_argument('-sc', '--sample-count', type=int, default=2,
                        help='Number of samples to generate for each output coverage')
    parser.add_argument('-t', '--threads', type=int, default=1, help='Number of compression threads')
    parser.add_argument('-s', '--seed', type=int, default=0, help='Random seed')
    args = parser.parse_args()

    random.seed(args.seed)

    # Sort output coverages
    args.output_coverages.sort()

    # Calculate read fraction for each output coverage
    read_fractions = [oc * args.sample_count / args.input_coverage for oc in args.output_coverages]
    if max(read_fractions) > 1:
        raise ValueError('Output coverages exceed input coverage')

    # Open input file
    input_file = pysam.AlignmentFile(args.input, threads=args.threads)
    # Get output file extension
    write_mode = 'wc' if input_file.is_cram else 'wb' if input_file.is_bam else 'w'
    extension_type = 'cram' if input_file.is_cram else 'bam' if input_file.is_bam else 'sam'
    # Open output files
    output_file_list = []
    for j in range(args.sample_count):
        out_files = []
        for i in range(len(read_fractions)):
            filename = f'{args.output}{i}_{args.output_coverages[i]}X_{j}.{extension_type}'
            alignment_file = pysam.AlignmentFile(filename, write_mode, template=input_file, threads=args.threads)
            out_files.append(alignment_file)
        output_file_list.append(out_files)

    random_salt = str(random.randint(0, 2**32))
    # Read input file
    for read in input_file.fetch(until_eof=True):
        # Select a random output file
        read_hash = int(hashlib.md5((random_salt + read.query_name).encode()).hexdigest(), 16)
        read_fraction = read_hash % 1000 / 1000
        frac_i = bisect.bisect_right(read_fractions, read_fraction)
        read_choice = read_hash % len(output_file_list)
        output_files = output_file_list[read_choice][frac_i:]
        for output_file in output_files:
            # Write read to output file
            output_file.write(read)
