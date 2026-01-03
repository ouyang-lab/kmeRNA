#!/bin/bash
# MIT License
# Copyright (c) 2026 Eric Nels Pederson, University of Massachusetts Amherst
# Licensed under the MIT License. See LICENSE file in the project root for full license information.

# Comparing AUROC by DeLong script:
# 1. If using SLURM, confirm the following settings:
#SBATCH -c 2 # Number of cores per task
#SBATCH --mem=2G  # Requested Memory
#SBATCH -p cpu  # Partition
#SBATCH -t 8:00:00  # Job time limit
#SBATCH -o /path/to/log_file

kmeRNA_dir=/path/to/package/directory/
work_dir=${kmeRNA_dir}/kmeRNA_eRNA-paRNA_RICseq_example_data/
script_dir=${kmeRNA_dir}/src

cell_type=HeLa
python3 ${script_dir}/aggregate_delong.py \                                                                                                       --inputA \
	/path/to/output/${cell_type}.enhancers.split_1.prob.tsv \
	/path/to/output/${cell_type}.enhancers.split_2.prob.tsv \
	/path/to/output/${cell_type}.enhancers.split_3.prob.tsv \
	--inputB \
	/path/to/output/${cell_type}.promoters.split_1.prob.tsv \
        /path/to/output/${cell_type}.promoters.split_2.prob.tsv \
        /path/to/output/${cell_type}.promoters.split_3.prob.tsv \
	|\
	tee /path/to/output/${cell_type}.enhancers_v_promoters.delong.txt





