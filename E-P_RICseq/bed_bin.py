#!/usr/bin/env python3

# MIT License
# Copyright (c) 2025 Eric Nels Pederson, University of Massachusetts Amherst
# Licensed under the MIT License. See LICENSE file in the project root for full license information.

import argparse
import math

def round_to_bin(value, bin_size):
    return round(value / bin_size) * bin_size

def round_start_to_bin(value, bin_size):
    return (value // bin_size) * bin_size

def round_end_to_bin(value, bin_size):
    return math.ceil(value / bin_size) * bin_size


def process_bed_file(input_bed, output_bed, bin_size):
    with open(input_bed, 'r') as infile, open(output_bed, 'w') as outfile:
        for line in infile:
            if line.startswith("#") or line.strip() == "":
                outfile.write(line)
                continue

            parts = line.strip().split('\t')
            if len(parts) < 6:
                continue  # Skip malformed lines

            chrom = parts[0]
            start = int(parts[1])
            end = int(parts[2])
            name = parts[3]
            score = parts[4]
            strand = parts[5]

            # Create bins from start to end at bin_size intervals
            for bin_start in range(start, end, bin_size):
                bin_end = min(bin_start + bin_size, end)

                bin_parts = parts[:]
                bin_parts[1] = str(bin_start)
                bin_parts[2] = str(bin_end)

                outfile.write('\t'.join(bin_parts) + '\n')

def process_bed_file_v1(input_bed, output_bed, bin_size):
    with open(input_bed, 'r') as infile, open(output_bed, 'w') as outfile:
        for line in infile:
            if line.startswith("#") or line.strip() == "":
                outfile.write(line)
                continue

            parts = line.strip().split('\t')
            #if len(parts) < 3:
            #    continue  # Skip malformed lines

            chrom = parts[0]
            start = int(parts[1])
            end = int(parts[2])
            name = str(parts[3])
            score = str(parts[4])
            strand = str(parts[5])

            rounded_start = round_start_to_bin(start, bin_size)
            #rounded_end = round_end_to_bin(end, bin_size)
            rounded_end = rounded_start + bin_size
            # Ensure start is never greater than end
            if rounded_start > rounded_end:
                rounded_start, rounded_end = rounded_end, rounded_start

            parts[1] = str(rounded_start)
            parts[2] = str(rounded_end)

            outfile.write('\t'.join(parts) + '\n')

def main():
    parser = argparse.ArgumentParser(description="Round BED file start and end positions to the nearest bin size.")
    parser.add_argument("--input_bed", help="Input BED file")
    parser.add_argument("--output_bed", help="Output BED file")
    parser.add_argument("--bin_size", type=int, help="Bin size to round positions to")
    args = parser.parse_args()

    process_bed_file(args.input_bed, args.output_bed, args.bin_size)

if __name__ == "__main__":
    main()

