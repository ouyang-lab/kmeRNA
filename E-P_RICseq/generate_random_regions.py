#!/usr/bin/env python3

# MIT License
# Copyright (c) 2025 Eric Nels Pederson, University of Massachusetts Amherst
# Licensed under the MIT License. See LICENSE file in the project root for full license information.

import argparse
import random

def load_bed(file_path):
    """Load BED file and store regions as a list of tuples."""
    regions = []
    with open(file_path, 'r') as f:
        for line in f:
            if line.strip():
                parts = line.strip().split('\t')
                chrom, start, end = parts[:3]
                regions.append((chrom, int(start), int(end)))
    return regions

def generate_samples_unbin(regions, n_samples, width, seed):
    """Generate random samples within the BED regions."""
    random.seed(seed)
    samples = []

    for _ in range(n_samples):
        chrom, start, end = random.choice(regions)
        if end - start < width:
            continue  # Skip if region is too small
        sample_start = random.randint(start, end - width)
        sample_end = sample_start + width
        samples.append((chrom, sample_start, sample_end))

    return samples

def generate_samples_rand(regions, n_samples, width, seed):
    """Generate random samples within the BED regions, binned to width-sized bins."""
    random.seed(seed)
    samples = []
    
    for _ in range(n_samples):
        chrom, start, end = random.choice(regions)
        if end - start < width:
            continue  # Skip if region is too small
        
        # Compute the number of bins within the region
        num_bins = (end - start) // width
        if num_bins == 0:
            continue
        
        # Select a random bin
        bin_index = random.randint(0, num_bins - 1)
        sample_start = start + bin_index * width
        sample_end = sample_start + width + 1
        samples.append((chrom, sample_start, sample_end))
    
    return samples

def generate_samples_with_dup(regions, n_samples, width, seed):
    """Generate random samples within the BED regions, rounded to the nearest width-sized interval."""
    random.seed(seed)
    samples = []
    for _ in range(n_samples):
        chrom, start, end = random.choice(regions)
        if end - start < width:
            continue  # Skip if region is too small
        sample_start = random.randint(start, end - width)
        sample_start = (sample_start // width) * width  # Round to nearest bin
        sample_end = sample_start + width
        sample_strand = random.choice(['+','-'])
        samples.append((chrom, sample_start, sample_end,sample_strand))
    return samples


def generate_samples(regions, n_samples, width, seed):
    """Generate unique random samples within the BED regions, rounded to the nearest width-sized interval."""
    random.seed(seed)
    samples = []
    seen = set()  # Keep track of unique (chrom, start, end)

    attempts = 0  # Prevent infinite loops in case of limited unique possibilities
    max_attempts = n_samples * 10  # Arbitrary limit

    while len(samples) < n_samples and attempts < max_attempts:
        chrom, start, end = random.choice(regions)
        if end - start < width:
            attempts += 1
            continue  # Skip if region is too small
        sample_start = random.randint(start, end - width)
        sample_start = (sample_start // width) * width  # Round to nearest bin
        sample_end = sample_start + width
        key = (chrom, sample_start, sample_end)
        if key in seen:
            attempts += 1
            continue  # Skip if already sampled
        sample_strand = random.choice(['+', '-'])
        samples.append((chrom, sample_start, sample_end, sample_strand))
        seen.add(key)
        attempts += 1

    return samples


def write_bed(samples, output_file):
    """Write generated samples to an output BED file."""
    idx=0
    with open(output_file, 'w') as f:
        for chrom, start, end, strand in samples:
            idx+=1
            f.write(f"{chrom}\t{start}\t{end}\tneg_region_{idx}\t.\t{strand}\n")

def main():
    parser = argparse.ArgumentParser(description="Generate random BED samples.")
    parser.add_argument("--bed", required=True, help="Input BED file")
    parser.add_argument("--n_samples", type=int, required=True, help="Number of samples to generate")
    parser.add_argument("--width", type=int, required=True, help="Width of each sample")
    parser.add_argument("--seed", type=int, required=True, help="Random seed for reproducibility")
    parser.add_argument("--output", required=True, help="Output BED file")

    args = parser.parse_args()

    regions = load_bed(args.bed)
    samples = generate_samples(regions, args.n_samples, args.width, args.seed)
    write_bed(samples, args.output)

if __name__ == "__main__":
    main()

