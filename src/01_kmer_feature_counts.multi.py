#!/usr/bin/env python3
# MIT License
# Copyright (c) 2026 Eric Nels Pederson, University of Massachusetts Amherst
# Licensed under the MIT License. See LICENSE file in the project root for full license information.

from itertools import product
import argparse
import pickle
import numpy as np
import csv
import gzip
import time
import subprocess
import tempfile
import os
import sys
import random
from concurrent.futures import ProcessPoolExecutor, as_completed

def parse_args():
    parser = argparse.ArgumentParser(description="RNA Sequence Embedding with per-kmer Z-scores + RNAduplex/RNAcofold energies")
    parser.add_argument("--input", type=str, help="Input .tsv or .csv file with ID, RNA1, RNA2")
    parser.add_argument("--output", type=str, help="Output pickle file for embedded weights")
    parser.add_argument("--out_format", type=str, default='pkl', choices=['pkl','csv','pkl,csv'])
    parser.add_argument("--max_k", type=int, help="Max k-mer length in bases")
    parser.add_argument("--min_k", default=1, type=int, help="Min k-mer length in bases")
    parser.add_argument("--type",type=str,default='Zscore', choices=['count','Zscore'])
    parser.add_argument("-n",default=1, type=int,help="Number of cores")
    parser.add_argument("--verbose", action='store_true')
    return parser.parse_args()

def read_file(filename):
    ids, seq1, seq2 = [], [], []
    open_func = gzip.open if filename.endswith('.gz') else open
    mode = 'rt' if filename.endswith('.gz') else 'r'
    with open_func(filename, mode, newline='\n') as file:
        reader = csv.reader(file, delimiter='\t')
        for row in reader:
            if 'id' in row:
                continue
            if len(row) >= 3:
                ids.append(str(row[0]))
                seq1.append(str(row[1]).replace('T', 'U'))
                seq2.append(str(row[2]).replace('T', 'U'))
            else:
                print(f"Ignoring row with unexpected number of columns: {row}")
                sys.stdout.flush()
    return ids, seq1, seq2

def initialize_kmers(min_k, max_k):
    bases = ['A', 'C', 'G', 'U']
    kmers = []
    for length in range(min_k, max_k + 1):
        kmers.extend([''.join(p) for p in product(bases, repeat=length)])
    return {kmer: 0 for kmer in kmers}

def count_kmers(seq, min_k, max_k, kmer_init):
    kmers = kmer_init.copy()
    seq_len = len(seq)
    total = 0

    for k in range(min_k, max_k + 1):
        for i in range(seq_len - k + 1):  # ← uses k directly
            kmer = seq[i:i+k]
            if kmer in kmers:
                kmers[kmer] += 1
                total += 1

    #if total > 0:
    #    for k in kmers:
    #        kmers[k] /= total
    max_count = max(kmers.values()) if kmers else 1
    if max_count > 0:
        for k in kmers:
            kmers[k] /= max_count

    return kmers


def process_pair(args):
    """Worker function for one (s1, s2) pair."""
    s1, s2, min_k, max_k, kmer_init = args

    freqs1 = count_kmers(s1, min_k, max_k, kmer_init)
    freqs2 = count_kmers(s2, min_k, max_k, kmer_init)

    duplex_energy = run_rnaduplex(s1, s2)
    cofold_energy = run_rnacofold(s1, s2)

    return freqs1, freqs2, duplex_energy, cofold_energy

def shuffle_embedded_output(embedded_out, seed=None):
    """
    Shuffle the order of rows (ID-feature pairs) in embedded_out.

    Args:
        embedded_out (list): flat list [id1, vec1, id2, vec2, ...]
        seed (int, optional): random seed for reproducibility

    Returns:
        list: same structure [id, vec, id, vec, ...] but shuffled order
    """
    if seed is not None:
        random.seed(seed)

    # Split into rows (ID, feature_vector)
    paired = [(embedded_out[i], embedded_out[i+1]) for i in range(0, len(embedded_out), 2)]

    # Shuffle the order of rows
    for idx in range(1,10):
        random.shuffle(paired)

    # Flatten back into original structure
    shuffled = []
    for id_, vec in paired:
        shuffled += [id_, vec]

    return shuffled

def compute_zscores(kmer_counts_list):
    """Compute per-kmer Z-scores across all samples."""
    all_kmers = sorted(list(kmer_counts_list[0].keys()))
    counts = np.array([[sample[k] for k in all_kmers] for sample in kmer_counts_list])
    means = np.mean(counts, axis=0)
    stds = np.std(counts, axis=0, ddof=1)
    stds[stds == 0] = 1e-8
    zscores = (counts - means) / stds
    return all_kmers, zscores

def run_rnaduplex(seq1, seq2):
    """Compute minimum free energy (MFE) using RNAduplex."""
    seq1 = seq1.upper().replace('T', 'U')
    seq2 = seq2.upper().replace('T', 'U')
    try:
        cmd = ["RNAduplex"]
        process = subprocess.run(cmd, input=f"{seq1}\n{seq2}\n",
                                 universal_newlines=True,
                                 stdout=subprocess.PIPE,   # replaces capture_output=True
                                 stderr=subprocess.PIPE,
                                 check=True)
        output = process.stdout.strip().split('\n')[-1]
        if '(' in output and ')' in output:
            energy = float(output.split('(')[-1].split(')')[0])
        else:
            energy = float(0)
        return energy

    except Exception as e:
        print("RNAduplex error:", e)
        return float(0)

def run_rnacofold(seq1, seq2):
    """Compute minimum free energy (MFE) using RNAcofold."""
    seq1 = seq1.upper().replace('T', 'U')
    seq2 = seq2.upper().replace('T', 'U')
    try:
        cmd = ["RNAcofold","--noPS"]
        process = subprocess.run(cmd, input=f"{seq1}&{seq2}\n",
                                 universal_newlines=True,
                                 stdout=subprocess.PIPE,   # replaces capture_output=True
                                 stderr=subprocess.PIPE,
                                 check=True)
        output = process.stdout.strip().split('\n')[-1]
        if '(' in output and ')' in output:
            energy = float(output.split('(')[-1].split(')')[0])
        else:
            energy = float(0)
        return energy
    except Exception as e:
        print("RNAcofold error:", e)
        return float(0)

def main():
    args = parse_args()
    start_time = time.time()

    ids, seq1, seq2 = read_file(args.input)
    if args.verbose:
        print(f"Read {len(ids)} sequence pairs")
        sys.stdout.flush()
    kmer_init = initialize_kmers(args.min_k, args.max_k)

    # --- Count kmers (normalized per sequence) ---
    #kmer_freq_list_s1 = []
    #kmer_freq_list_s2 = []
    #duplex_energies = []
    #cofold_energies = []
    # --- Parallel processing of each (s1, s2) pair ---
    kmer_freq_list_s1, kmer_freq_list_s2 = [], []
    duplex_energies, cofold_energies = [], []
    # Create argument list for all workers
    task_args = [(s1, s2, args.min_k, args.max_k, kmer_init) for s1, s2 in zip(seq1, seq2)]

    # Adjust number of workers if needed (default: # of CPUs)
    #with ProcessPoolExecutor(max_workers=args.n) as executor:
    #    futures = [executor.submit(process_pair, a) for a in task_args]
    #    for i, f in enumerate(as_completed(futures)):
    #        freqs1, freqs2, duplex_e, cofold_e = f.result()
    #        kmer_freq_list_s1.append(freqs1)
    #        kmer_freq_list_s2.append(freqs2)
    #        duplex_energies.append(duplex_e)
    #        cofold_energies.append(cofold_e)
    #        if args.verbose and (i+1) % 10 == 0:
    #            print(f"Processed {i+1}/{len(task_args)} samples")
    with ProcessPoolExecutor(max_workers=args.n) as executor:
        futures = [executor.submit(process_pair, a) for a in task_args]

        # Retrieve results in submission order
        for i, f in enumerate(futures):
            freqs1, freqs2, duplex_e, cofold_e = f.result()

            kmer_freq_list_s1.append(freqs1)
            kmer_freq_list_s2.append(freqs2)
            duplex_energies.append(duplex_e)
            cofold_energies.append(cofold_e)

            if args.verbose and (i + 1) % 10 == 0:
                print(f"Processed {i+1}/{len(task_args)} samples")

    if args.type == "Zscore":
        # --- Compute per-kmer Z-scores across all samples ---
        all_kmers, zscores_s1 = compute_zscores(kmer_freq_list_s1)
        _, zscores_s2 = compute_zscores(kmer_freq_list_s2)
    
        # --- Combine features: [Z-scores + RNAduplex + RNAcofold] ---
        embedded_out = []
        for idx, id_ in enumerate(ids):
            feature_vector = (
                zscores_s1[idx, :].tolist() +
                zscores_s2[idx, :].tolist() +
                [float(duplex_energies[idx]), float(cofold_energies[idx])]
            )
            embedded_out += [id_, feature_vector]
    elif args.type == "count":
        # k-mer frequencies (normalized counts by total # k-mers)
        #feats_s1 = kmer_freq_list_s1
        #feats_s2 = kmer_freq_list_s2
        
        #feats_s1 = [[float(x) for x in arr] for arr in kmer_freq_list_s1]
        #feats_s2 = [[float(x) for x in arr] for arr in kmer_freq_list_s2]
        all_kmers = sorted(list(kmer_freq_list_s1[0].keys()))
        feats_s1 = np.array([[sample[k] for k in all_kmers] for sample in kmer_freq_list_s1])
        feats_s2 = np.array([[sample[k] for k in all_kmers] for sample in kmer_freq_list_s2])
        
        embedded_out = []
        for idx, id_ in enumerate(ids):
            feature_vector = (
                list(feats_s1[idx]) +
                list(feats_s2[idx]) +
                [float(duplex_energies[idx]), float(cofold_energies[idx])]
            )
            embedded_out += [id_, feature_vector]
    else:
        raise ValueError(f"Unknown feature type: {args.type}")

    # --- Save outputs ---
    if 'pkl' in args.out_format:
        with open(args.output + '.pkl', 'wb') as f:
            pickle.dump(embedded_out, f)

    if 'csv' in args.out_format:
        with open(args.output + '.csv', 'w', newline='') as f:
            writer = csv.writer(f)
            header_s1 = [f"{k}_s1" for k in all_kmers]
            header_s2 = [f"{k}_s2" for k in all_kmers]
            writer.writerow(['ID'] + header_s1 + header_s2 + ['RNAduplex_energy', 'RNAcofold_energy'])
            for idx, id_ in enumerate(ids):
                writer.writerow([id_] + zscores_s1[idx, :].tolist() + zscores_s2[idx, :].tolist() +
                            [duplex_energies[idx], cofold_energies[idx]])
    if args.verbose:
        print(f"Completed in {time.time() - start_time:.2f}s")
        sys.stdout.flush()
if __name__ == "__main__":
    main()

