#!/usr/bin/env python3
# MIT License
# Copyright (c) 2025 Eric Nels Pederson, University of Massachusetts Amherst
# Licensed under the MIT License. See LICENSE file in the project root for full license information.
from itertools import product, permutations, combinations_with_replacement
import argparse
import pickle
import numpy as np
import csv
import gzip
import time

def parse_args():
    parser = argparse.ArgumentParser(description="RNA Sequence Embedding using normalized k-mer and contiguous k-mer pair counts")
    parser.add_argument("--input", type=str, help="Input .tsv file containing an identifier (col 1) RNA sequences 1 & 2 (col 2-3)")
    parser.add_argument("--output", type=str, help="Output pickle file for embedded weights")
    parser.add_argument("--out_format", type=str, help="Output pickle file for embedded weights",default='pkl',choices=['pkl','csv','pkl,csv'])
    #parser.add_argument("--label", type=int,help="0/1 label",choices=[0,1])
    parser.add_argument("--max_k", type=int, help="Max k-mer length in bases")
    parser.add_argument("--min_k", default=1,type=int, help="Min k-mer length in bases")
    parser.add_argument("--verbose", action='store_true')
    return parser.parse_args()


def read_file(filename):
    ids = []
    seq1 = []
    seq2 = []
    if filename.endswith('.gz'):
        open_func = lambda f: gzip.open(f, 'rt', newline='\n')
    else:
        open_func = lambda f: open(f, 'r', newline='\n')
    with open_func(filename) as file:
        reader = csv.reader(file, delimiter='\t')
        for row in reader:
            if len(row) == 3:
                ids.append(str(row[0]))
                seq1.append(str(row[1]))
                seq2.append(str(row[2]))
            elif len(row) > 3:
                ids.append(str(row[0]))
                seq1.append(str(row[1]))
                seq2.append(str(row[2]))
            else:
                print(f"Ignoring row with unexpected number of columns: {row}") 
    return ids, seq1, seq2 



def initialize_kmers(min_k,max_k):
    """Initialize all possible kmers with count 0"""
    bases = ['A', 'C', 'G', 'T']
    kmers = []
    for length in range(min_k, max_k + 1):
        kmers.extend([''.join(p) for p in product(bases, repeat=length)])
    return {kmer: 0 for kmer in kmers}


def count_kmers(kmer_init,seq,min_k, max_k):
    #kmers = initialize_kmers(k)  # Initialize all possible kmers with count 0
    kmers = kmer_init
    # Iterate through the sequence to count kmers
    for i in range(len(seq) - max_k + 1):
        for j in range(min_k, max_k + 1):
            kmer = seq[i:i+j]
            if kmer in kmers:
                kmers[kmer] += 1
    
    # Sort kmers alphabetically and create a list of tuples (kmer, count)
    sorted_kmers = sorted(kmers.items())
    return sorted_kmers





def main():
    args = parse_args()
    if args.input:
        ids, input_sequence1, input_sequence2 = read_file(args.input)
        num_samples=len(input_sequence1)
        kmer_count_vecs = {}
        begin = time.time()
        embedded_out = []
        kmer_init = initialize_kmers(args.min_k,args.max_k)
        for idx in range(len(input_sequence1)):
            start = time.time()
            kmer_init1 = kmer_init.copy()
            kmer_init2 = kmer_init.copy()
            #kmer_init1 = initialize_kmers(args.min_k,args.max_k)
            #kmer_init2 = initialize_kmers(args.min_k,args.max_k)
            kmer_count_vec = count_kmers(kmer_init1,input_sequence1[idx], args.min_k, args.max_k)
            kmer_count_vec += count_kmers(kmer_init2,input_sequence2[idx], args.min_k, args.max_k)
            kmer_count_vec_out = [t[1] for t in kmer_count_vec]
            if max(kmer_count_vec_out)>0:
                kmer_counts_out_norm = [x / max(kmer_count_vec_out) for x in kmer_count_vec_out] 
            else:
                kmer_counts_out_norm = [x / 1 for x in kmer_count_vec_out]
            embedded_out += [ids[idx],kmer_counts_out_norm]
            if args.out_format == 'csv' or args.out_format == 'pkl,csv':
                with gzip.open(args.output+'.csv.gz','at',encoding='utf-8') as f:
                    f.write(str([ids[idx],kmer_counts_out_norm])+'\n')
            end = time.time()
            if args.verbose:
                if idx==0:
                    kmer_len = len(kmer_init1)
                    print(f"{args.min_k}\t= min_k")
                    print(f"{args.max_k}\t= max_k")
                    print(f"{kmer_len}\tk-mers")
                    #print(f"{contact_len}\tk-mer pairs")
                    print(f"{2*kmer_len}\tTotal vector length")
                if (idx+1)%100==0:
                    x = (idx+1)/num_samples*100
                    print(f"{idx+1}/{num_samples} lines . . . {x:.2f} % Complete . . . {round((end-start),1)} min | {round((end-begin)/60,2)} min elapsed")
                if idx+1==(num_samples):
                    print(f"{idx+1}/{num_samples} lines . . . 100 % Complete . . . {round((end-begin)/60,2)} min elapsed")

            #embedded_out = []
        if args.input is None or args.output is None or args.max_k is None:
            print("Incomplete number of arguments given")
            exit()
        if args.out_format == 'pkl' or args.out_format == 'pkl,csv':
            with open(args.output+'.pkl', 'wb') as f:
                pickle.dump(embedded_out, f)


if __name__ == "__main__":
    main()


