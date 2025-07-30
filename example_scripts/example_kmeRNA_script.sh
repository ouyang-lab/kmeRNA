#!/bin/bash
# MIT License
# Copyright (c) 2025 Eric Nels Pederson, University of Massachusetts Amherst
# Licensed under the MIT License. See LICENSE file in the project root for full license information.

# Preparing script:
# 1. If using SLURM, confirm the following settings:
#SBATCH -c 32 # Number of cores per task
#SBATCH --mem=60G  # Requested Memory
#SBATCH -p cpu  # Partition
#SBATCH -t 8:00:00  # Job time limit
#SBATCH -o /path/to/log_file

# 2. Set seed, k-mer range and max number of epochs (for deep learning only)
SEED_VAL=42813
RANDOM=${SEED_VAL}
Kmin=1
K=5
max_epochs=20

# 3. Cell type and query region type (enhancers, promoters or DN.enhancer-promoters):
cell_type=HeLa
region_type=enhancers

# 4. adjust absolute paths
# Directory paths (requires setup)
kmeRNA_dir=/path/to/package/directory/
work_dir=${kmeRNA_dir}/kmeRNA_eRNA-paRNA_RICseq_example_data/
script_dir=${kmeRNA_dir}/src

# File paths (requires setup). Suffixes will be added to input and output data
train_input_dat=${work_dir}/split_train/embeddings/GSE190214_HeLa.Enhancer_Promoter.chimeric_fragment.out.150nt.split_train.${region_type}
test_input_dat=${work_dir}/split_test/embeddings/GSE190214_HeLa.Enhancer_Promoter.chimeric_fragment.out.150nt.split_test.${region_type}
validate_input_dat=${work_dir}/split_validate/embeddings/GSE190214_HeLa.Enhancer_Promoter.chimeric_fragment.out.150nt.split_validate.${region_type}

train_labels=${work_dir}/split_train/embeddings/GSE190214_HeLa.Enhancer_Promoter.chimeric_fragment.out.150nt.split_train.seed${SEED_VAL}.${region_type}.kmeRNA.labels.tsv
test_labels=${work_dir}/split_test/embeddings/GSE190214_HeLa.Enhancer_Promoter.chimeric_fragment.out.150nt.split_test.seed${SEED_VAL}.${region_type}.kmeRNA.labels.tsv
validate_labels=${work_dir}/split_validate/embeddings/GSE190214_HeLa.Enhancer_Promoter.chimeric_fragment.out.150nt.split_validate.seed${SEED_VAL}.${region_type}.kmeRNA.labels.tsv

train_output_dat=${work_dir}/split_train/results/GSE190214_HeLa.Enhancer_Promoter.chimeric_fragment.out.150nt.split_train.seed${SEED_VAL}.${region_type}
test_output_dat=${work_dir}/split_test/results/GSE190214_HeLa.Enhancer_Promoter.chimeric_fragment.out.150nt.split_test.seed${SEED_VAL}.${region_type}
validate_output_dat=${work_dir}/split_validate/GSE190214_HeLa.Enhancer_Promoter.chimeric_fragment.out.150nt.split_validate.seed${SEED_VAL}.${region_type}



# 5. Driver options: set values to 0 to skip and 1 to run
create_labels=1
embedding=1
deep_learning=1
deep_learning_shap=1
tree_learning=1
tree_shap=1

# use trained model to predict across cell types:
cross_cell_deep_learning=1
cross_cell_trees=1

# Create labels
if [ $create_labels -gt 0 ]
then
	sh ${script_dir}/00_extract_labels.sh \
		--input_data ${train_input_dat} \
		--labels ${train_labels}
	sh ${script_dir}/00_extract_labels.sh \
		--input_data ${test_input_dat} \
		--labels ${test_labels}
	sh ${script_dir}/00_extract_labels.sh \
		--input_data ${validate_input_dat} \
		--labels ${validate_labels}

fi

# Embedding
if [ $embedding -gt 0 ]
then
	python3 ${script_dir}/01_kmer_feature_counts.py	 \
		--input ${train_input_dat}.tsv.gz \
		--output ${train_input_dat}.k${Kmin}_k${K}.embed.out \
		--out_format pkl \
		--min_k ${Kmin} \
		--max_k ${K} \
		--verbose
	python3 ${script_dir}/01_kmer_feature_counts.py	 \
		--input ${test_input_dat}.tsv.gz \
		--output ${test_input_dat}.k${Kmin}_k${K}.embed.out \
		--out_format pkl \
		--min_k ${Kmin} \
		--max_k ${K} \
		--verbose
	python3 ${script_dir}/01_kmer_feature_counts.py	 \
		--input ${validate_input_dat}.tsv.gz \
		--output ${validate_input_dat}.k${Kmin}_k${K}.embed.out \
		--out_format pkl \
		--min_k ${Kmin} \
		--max_k ${K} \
		--verbose



fi

# Deep learning model training
if [ $deep_learning -gt 0 ]
then
	python3 ${script_dir}/02_train_test_kmer_keras.py \
		--task RNA-RNA \
                --train_data ${train_input_dat}.k${Kmin}_k${K}.embed.out.pkl \
                --train_labels $train_labels \
                --test_data ${validate_input_dat}.k${Kmin}_k${K}.embed.out.pkl \
                --test_labels ${validate_labels} \
                --cv \
                --model_out ${train_output_dat}.model.k${Kmin}_k${K} \
                --max_epochs ${max_epochs} \
                --min_k ${Kmin} \
                --max_k ${K} \
                --seed_val ${SEED_VAL} \
                --batch_size 32 \
                --drop1 3 \
                --drop2 2 \
                --drop3 1 \
                --drop4 1 \
                --output_train_prob ${train_output_dat}.k${Kmin}_k${K}.prob.tsv \
                --plot_path ${train_output_dat}.k${Kmin}_k${K}.train_stats

	python3 ${script_dir}/02_train_test_kmer_keras.py \
                --task RNA-RNA \
                --test_data ${test_input_dat}.k${Kmin}_k${K}.embed.out.pkl \
                --test_labels $test_labels \
                --model_in ${train_output_dat}.model.k${Kmin}_k${K} \
                --batch_size 4 \
                --output_test_prob ${test_output_dat}.pairs.RNA.k${Kmin}_k${K}.prob.tsv

        python3 ${script_dir}/get_roc_auc.py \
                --input ${test_output_dat}.pairs.RNA.k${Kmin}_k${K}.prob.tsv \
                --plot_roc ${test_output_dat}.pairs.RNA.k${Kmin}_k${K}.plot |\
                tee ${test_output_dat}.pairs.RNA.k${Kmin}_k${K}.perf.tsv
fi

# Deep learning model feature importance
if [ $deep_learning_shap -gt 0 ]
then
	model_in=${train_output_dat}.model.k${Kmin}_k${K}.keras
        python3 ${script_dir}/04_get_SHAP_diff_values.py \
                --model ${model_in} \
                --data ${train_input_dat}.k${Kmin}_k${K}.embed.out.pkl \
                --output ${train_output_dat}.k${Kmin}_k${K}.embed.out.SHAP_diff \
		--min_k ${K_min} \
		--max_k ${K}

fi

# Tree-based model training
if [ $tree_learning -gt 0 ]
then
	echo "Training ensemble models..."
        python3 ${script_dir}/03_train_test_kmer_ensemble.py \
                --train_data ${train_input_dat}.k${K_min}_k${K}.embed.out.pkl \
                --train_labels $train_labels \
                --model_out ${train_output_dat}.RNA.k${K_min}_k${K} \
                --seed_val ${SEED_VAL} \
                --output_train_prob ${train_output_dat}.k${K_min}_k${K}.prob.tsv
	echo "Loading trained model and testing set..."
        python3 ${script_dir}/03_train_test_kmer_ensemble.py \
                --model_in ${train_output_dat}.RNA.k${K_min}_k${K} \
                --seed_val ${SEED_VAL} \
                --test_data ${test_input_dat}.k${K_min}_k${K}.embed.out.pkl \
                --test_labels $test_labels \
                --output_test_prob ${test_output_dat}.k${K_min}_k${K}.prob.tsv
        echo "Loading trained model and validation set..."
        python3 ${script_dir}/03_train_test_kmer_ensemble.py \
                --model_in ${train_output_dat}.RNA.k${K_min}_k${K} \
                --seed_val ${SEED_VAL} \
                --test_data ${validate_input_dat}.k${K_min}_k${K}.embed.out.pkl \
                --test_labels $validate_labels \
                --output_test_prob ${validate_output_dat}.k${K_min}_k${K}.prob.tsv
        python3 ${script_dir}/get_roc_auc.ens.py \
                --input ${train_output_dat}.k${K_min}_k${K}.prob.tsv \
                --plot_roc ${train_output_dat}.k${K_min}_k${K}.plot |\
		tee ${train_output_dat}.k${K_min}_k${K}.perf.tsv
        python3 ${script_dir}/get_roc_auc.ens.py \
                --input ${test_output_dat}.k${K_min}_k${K}.prob.tsv \
                --plot_roc ${test_output_dat}.k${K_min}_k${K}.plot |\
		tee ${test_output_dat}.k${K_min}_k${K}.perf.tsv
	python3 ${script_dir}/get_roc_auc.ens.py \
                --input ${validate_output_dat}.k${K_min}_k${K}.prob.tsv \
                --plot_roc ${validate_output_dat}.k${K_min}_k${K}.plot |\
		tee ${validate_output_dat}.k${K_min}_k${K}.perf.tsv

fi
# Complete Tree-based model feature importance (SHAP)
if [ $tree_shap -gt 0 ]
then
	model_in=${train_output_dat}.RNA.k${Kmin}_k${K}.RF.joblib
        echo "Training data SHAP RandomForest model"
        python3 ${script_dir}/05_get_SHAP_values.py \
                --model ${model_in} \
                --tree \
                --data ${train_input_dat}.k${Kmin}_k${K}.embed.out.pkl \
                --output ${train_output_dat}.k${Kmin}_k${K}.embed.out.RF.SHAP_diff

        model_in=${train_output_dat}.RNA.k${Kmin}_k${K}.ET.joblib
        echo "Training data SHAP ExtraTrees model"
        python3 ${script_dir}/05_get_SHAP_values.py \
                --model ${model_in} \
                --tree \
                --data ${train_input_dat}.k${Kmin}_k${K}.embed.out.pkl \
                --output ${train_output_dat}.k${Kmin}_k${K}.embed.out.ET.SHAP_diff
        
	model_in=${train_output_dat}.RNA.k${Kmin}_k${K}.GB.joblib
        echo "Training data SHAP GradientBoosting model"
        python3 ${script_dir}/05_get_SHAP_values.py \
                --model ${model_in} \
                --tree \
                --data ${train_input_dat}.k${Kmin}_k${K}.embed.out.pkl \
                --output ${train_output_dat}.k${Kmin}_k${K}.embed.out.GB.SHAP_diff


fi

# Cross cell (deep learning)
if [[ $cross_cell_deep_learning -gt 0 ]]
then
        model_in=${train_output_dat}.model.k${Kmin}_k${K}.keras
	for cell_type2 in HeLa hNPC GM12878 IMR90 H1 HepG2 K562
        do
                if [[ $cell_type == $cell_type2 ]]
                then
                        continue
                fi
        	test_input_cross_cell=${work_dir}/split_test/embeddings/GSE190214_${cell_type}-${cell_type2}.Enhancer_Promoter.chimeric_fragment.out.150nt.split_test.${region_type}
        	test_labels_cross_cell=${work_dir}/split_test/embeddings/GSE190214_${cell_type}-${cell_type2}.Enhancer_Promoter.chimeric_fragment.out.150nt.split_test.seed${SEED_VAL}.${region_type}.labels.tsv
		test_output_cross_cell=${work_dir}/split_test/results/GSE190214_${cell_type}-${cell_type2}.Enhancer_Promoter.chimeric_fragment.out.150nt.split_test.seed${SEED_VAL}.${region_type}
		
		zcat ${test_input_cross_cell}.tsv.gz |\
                        awk '{OFS="\t"}{if ($1~/neg/) {print 0} else {print 1}}' >\
                ${test_labels_cross_cell}
	                python3 ${script_dir}/01_kmer_feature_counts.py \
                        	--input ${test_input_cross_cell}.tsv.gz \
                       		--output ${test_input_cross_cell}.k${Kmin}_k${K}.embed.out.pkl \
                        	--out_format pkl \
				--min_k ${Kmin} \
                        	--max_k ${K}

		echo $model_in | awk '{gsub(/\//," "); print $NF}'
                echo ${test_input_cross_cell}
                python3 ${script_dir}/02_train_test_kmer_keras_predict.py \
                        --task RNA-RNA \
                        --test_data ${test_input_cross_cell}.k${Kmin}_k${K}.embed.out.pkl \
                        --test_labels $test_labels_cross_cell \
                        --model_in ${model_in} \
                        --batch_size 4 \
                        --output_test_prob ${test_output_cross_cell}.k${Kmin}_k${K}.prob.tsv
		
		python3 ${script_dir}/get_roc_auc.py \
                	--input ${test_output_cross_cell}.k${Kmin}_k${K}.prob.tsv \
                	--plot_roc ${test_output_cross_cell}.k${Kmin}_k${K}.plot |\
                	tee ${test_output_cross_cell}.k${Kmin}_k${K}.perf.tsv
	done


fi

# Cross cell (Tree-based models)
if [[ $cross_cell_trees -gt 0 ]]
then
        model_in=${train_output_dat}.model.k${Kmin}_k${K}.keras
	for cell_type2 in HeLa hNPC GM12878 IMR90 H1 HepG2 K562
        do
                if [[ $cell_type == $cell_type2 ]]
                then
                        continue
                fi
        	test_input_cross_cell=${work_dir}/split_test/embeddings/GSE190214_${cell_type}-${cell_type2}.Enhancer_Promoter.chimeric_fragment.out.150nt.split_test.${region_type}
        	test_labels_cross_cell=${work_dir}/split_test/embeddings/GSE190214_${cell_type}-${cell_type2}.Enhancer_Promoter.chimeric_fragment.out.150nt.split_test.seed${SEED_VAL}.${region_type}.labels.tsv
		test_output_cross_cell=${work_dir}/split_test/results/GSE190214_${cell_type}-${cell_type2}.Enhancer_Promoter.chimeric_fragment.out.150nt.split_test.seed${SEED_VAL}.${region_type}
		
		zcat ${test_input_cross_cell}.tsv.gz |\
                        awk '{OFS="\t"}{if ($1~/neg/) {print 0} else {print 1}}' >\
                ${test_labels_cross_cell}
	                python3 ${script_dir}/01_kmer_feature_counts.py \
                        	--input ${test_input_cross_cell}.tsv.gz \
                       		--output ${test_input_cross_cell}.k${Kmin}_k${K}.embed.out.pkl \
                        	--out_format pkl \
				--min_k ${Kmin} \
                        	--max_k ${K}
			model_in=${train_output_dat}.RNA.k${Kmin}_k${K}
			echo $model_in | awk '{gsub(/\//," "); print $NF}'
        	        echo ${test_input_cross_cell}		
			python3 ${script_dir}/03_train_test_kmer_ensemble.py \
        	                --model_in ${model_in} \
        	                --seed_val ${SEED_VAL} \
        	                --test_data ${test_input_cross_cell}.RNA.k${Kmin}_k${K}.embed.out.pkl \
        	                --test_labels $test_labels_cross_cell \
        	                --output_test_prob ${test_output_cross_cell}.RNA.k${Kmin}_k${K}.ens.prob.tsv

			python3 ${script_dir}/get_roc_auc.ens.py \
        	        	--input ${test_output_cross_cell}.k${Kmin}_k${K}.ens.prob.tsv \
        	        	--plot_roc ${test_output_cross_cell}.k${Kmin}_k${K}.ens.plot |\
        	        	tee ${test_output_cross_cell}.k${Kmin}_k${K}.ens.perf.tsv
	done


fi



