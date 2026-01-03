#!/usr/bin/env python3

# MIT License
# Copyright (c) 2026 Eric Nels Pederson, University of Massachusetts Amherst
# Licensed under the MIT License. See LICENSE file in the project root for full license information.

import argparse
import pandas as pd
import random
from collections import defaultdict

def load_data(file_path):
    df = pd.read_csv(file_path, sep='\t', header=None, names=['chrom1', 'chrom2', 'percent','count'])
    return df

def compute_chrom_weights(df):
    weights = defaultdict(float)
    for _, row in df.iterrows():
        weights[row['chrom1']] += row['percent']
        weights[row['chrom2']] += row['percent']
    return weights




def split_chromosomes_by_weight(chrom_weights, target_ratios=(0.6, 0.2, 0.2), seed=42):
    random.seed(seed)
    chroms = list(chrom_weights.items())
    random.shuffle(chroms)

    total = sum(weight for _, weight in chroms)
    targets = {
        'train': total * target_ratios[0],
        'test': total * target_ratios[1],
        'val': total * target_ratios[2],
    }

    buckets = {'train': set(), 'test': set(), 'val': set()}
    bucket_weights = {'train': 0, 'test': 0, 'val': 0}

    for chrom, weight in chroms:
        # Only consider buckets that are under their target
        available = {
            b: targets[b] - bucket_weights[b]
            for b in buckets
            if bucket_weights[b] + weight <= targets[b] * 1.05  # allow small overshoot
        }

        if not available:
            # If all buckets are full, assign to the one with least overshoot
            best_fit = min(buckets, key=lambda b: (bucket_weights[b] + weight) - targets[b])
        else:
            # Choose the bucket that is farthest under target (greedy fill)
            best_fit = max(available, key=available.get)

        buckets[best_fit].add(chrom)
        bucket_weights[best_fit] += weight

    return buckets['train'], buckets['test'], buckets['val']

def assign_pairs(df, train_set, test_set, val_set):
    categories = defaultdict(list)

    for _, row in df.iterrows():
        c1, c2 = row['chrom1'], row['chrom2']
        if c1 in train_set and c2 in train_set:
            categories['train'].append((c1, c2))
        elif c1 in test_set and c2 in test_set:
            categories['test'].append((c1, c2))
        elif c1 in val_set and c2 in val_set:
            categories['validate'].append((c1, c2))
        # else: ignore pair

    return categories

def write_output(categories, prefix):
    for split in ['train', 'test', 'validate']:
        with open(f'{prefix}_{split}.txt', 'w') as f:
            for pair in categories[split]:
                f.write(f'{pair[0]}\t{pair[1]}\n')

def main():
    parser = argparse.ArgumentParser(description="Optimally split chrom-pair list into train/test/validate sets.")
    parser.add_argument('--input_file', help="Path to input file (tab-separated with chrom1, chrom2, percent)")
    parser.add_argument('--output_prefix', default='output', help="Prefix for output files")
    parser.add_argument('--seed', type=int, default=98, help="Random seed for reproducibility")
    args = parser.parse_args()

    df = load_data(args.input_file)
    chrom_weights = compute_chrom_weights(df)

    train_set, test_set, val_set = split_chromosomes_by_weight(chrom_weights, seed=args.seed)
    categories = assign_pairs(df, train_set, test_set, val_set)
    write_output(categories, args.output_prefix)

if __name__ == "__main__":
    main()







