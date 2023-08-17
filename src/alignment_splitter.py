# Copyright 2022 - Barcelona Supercomputing Center
# Author: Rodrigo Martin
# BSC Dual License

import argparse
import hashlib
import random
from concurrent.futures import ProcessPoolExecutor
import pysam


def write_output(input_path: str, output_file: str, cov_fraction: float, sample_num: int, total_samples: int, random_salt: str, fasta_ref: str = None, threads: int = 1):
    input_file = pysam.AlignmentFile(input_path, threads=threads, reference_filename=fasta_ref)
    write_mode = 'wc' if input_file.is_cram else 'wb' if input_file.is_bam else 'w'
    output_file = pysam.AlignmentFile(output_file, write_mode, template=input_file, threads=threads)
    # Read input file
    for read in input_file.fetch(until_eof=True):
        # Select a random output file
        read_hash = int(hashlib.md5((random_salt + read.query_name).encode()).hexdigest(), 16)
        # Check if this read is in its coverage fraction
        read_fraction = read_hash % 1000 / 1000
        if read_fraction > cov_fraction:
            continue
        # Check if the read is in this sample
        sample_fraction = read_hash % total_samples
        if sample_fraction != sample_num:
            continue
        output_file.write(read)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', required=True, type=str, help='Input alignment file')
    parser.add_argument('-f', '--fasta', type=str, help='Reference fasta file (required for CRAM files)')
    parser.add_argument('-ic', '--input-coverage', required=True, type=int, help='Input file coverage')
    parser.add_argument('-o', '--output', required=True, type=str, help='Output prefix')
    parser.add_argument('-oc', '--output-coverages', required=True, nargs='+', type=int, help='Output coverages')
    parser.add_argument('-sc', '--sample-count', type=int, default=2,
                        help='Number of samples to generate for each output coverage')
    parser.add_argument('-p', '--max-processes', type=int, default=1, help='Number of max processes')
    parser.add_argument('-s', '--seed', type=int, default=0, help='Random seed')
    args = parser.parse_args()

    random.seed(args.seed)

    # Sort output coverages
    args.output_coverages.sort()

    # Calculate read fraction for each output coverage
    coverage_fractions = [oc * args.sample_count / args.input_coverage for oc in args.output_coverages]
    if max(coverage_fractions) > 1:
        raise ValueError('Output coverages exceed input coverage')

    # Get extension type
    input_file = pysam.AlignmentFile(args.input, reference_filename=args.fasta)
    extension_type = 'cram' if input_file.is_cram else 'bam' if input_file.is_bam else 'sam'
    input_file.close()

    # Calculate random salt
    random_salt = str(random.randint(0, 2**32))

    available_threads = max(args.max_processes - len(coverage_fractions) * args.sample_count, 1)
    pool = ProcessPoolExecutor(max_workers=args.max_processes)
    futures = []
    for i in range(len(coverage_fractions)):
        for j in range(args.sample_count):
            futures.append(pool.submit(write_output, args.input,
                                       f'{args.output}{i}_{args.output_coverages[i]}X_{j}.{extension_type}',
                                       coverage_fractions[i], j, args.sample_count,
                                       random_salt, args.fasta, available_threads))

    for future in futures:
        future.result()
    pool.shutdown()
