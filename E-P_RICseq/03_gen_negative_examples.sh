#!/bin/bash

# MIT License
# Copyright (c) 2026 Eric Nels Pederson, University of Massachusetts Amherst
# Licensed under the MIT License. See LICENSE file in the project root for full license information.

genome_faidx=/path/to/genome/Homo_sapiens.GRCh37.75.dna_sm.toplevel.fa.fai
genome_fasta=/path/to/genome/Homo_sapiens.GRCh37.75.dna_sm.toplevel.fa
GTF=/path/to/genome/Homo_sapiens.GRCh37.75.sort.gtf
seed=42813
bin_size=150
out_dir=/path/to/work_dir/

for cell_type in K562 IMR90 HeLa hNPC HepG2 GM12878 H1;
do
	train_bedpe=GSE190214_${cell_type}.Enhancer_Promoter.chimeric_fragment.out.150nt.EP_anno.uniq.v3.split_train.bedpe
	test_bedpe=GSE190214_${cell_type}.Enhancer_Promoter.chimeric_fragment.out.150nt.EP_anno.uniq.v3.split_test.bedpe
	validate_bedpe=GSE190214_${cell_type}.Enhancer_Promoter.chimeric_fragment.out.150nt.EP_anno.uniq.v3.split_validate.bedpe
	
	## Generate split negative region files
	awk '{OFS="\t"}FNR==NR{a[$1]+=1;next}{if (a["chr"$1]>0) {print $0}}' \
		<(awk '{print $1"\n"$2}' GSE190214_${cell_type}.Enhancer_Promoter.chimeric_fragment.out.150nt.EP_anno.uniq.v3.split_train.txt | sort | uniq) \
		<(zcat genome/Homo_sapiens.GRCh37.75.noncoding.regions.no_EP.bed.gz) >\
		Homo_sapiens.GRCh37.75.noncoding.regions.${cell_type}.v3.split_train.bed
	awk '{OFS="\t"}FNR==NR{a[$1]+=1;next}{if (a["chr"$1]>0) {print $0}}' \
		<(awk '{print $1"\n"$2}' GSE190214_${cell_type}.Enhancer_Promoter.chimeric_fragment.out.150nt.EP_anno.uniq.v3.split_test.txt | sort | uniq) \
		<(zcat genome/Homo_sapiens.GRCh37.75.noncoding.regions.no_EP.bed.gz) >\
		Homo_sapiens.GRCh37.75.noncoding.regions.${cell_type}.v3.split_test.bed
	awk '{OFS="\t"}FNR==NR{a[$1]+=1;next}{if (a["chr"$1]>0) {print $0}}' \
		<(awk '{print $1"\n"$2}' GSE190214_${cell_type}.Enhancer_Promoter.chimeric_fragment.out.150nt.EP_anno.uniq.v3.split_validate.txt | sort | uniq) \
		<(zcat genome/Homo_sapiens.GRCh37.75.noncoding.regions.no_EP.bed.gz) >\
		Homo_sapiens.GRCh37.75.noncoding.regions.${cell_type}.v3.split_validate.bed

	#exit
	cat GSE190214_${cell_type}.Enhancer_Promoter.chimeric_fragment.out.150nt.EP_anno.uniq.v3.split_train.bedpe |\
		awk '{OFS="\t"}{gsub(/chr/,"");printf "%s\t%.0f\t%.0f\t%s\t%s\t%s\n%s\t%.0f\t%.0f\t%s\t%s\t%s\n", $1,$2,$3,$7"=1",".",$9,$4,$5,$6,$7"=2",".",$10}' >\
		GSE190214_${cell_type}.Enhancer_Promoter.chimeric_fragment.out.150nt.EP_anno.uniq.v3.split_train.tmp.bed

	bedtools getfasta \
		-fi $genome_fasta \
		-bed GSE190214_${cell_type}.Enhancer_Promoter.chimeric_fragment.out.150nt.EP_anno.uniq.v3.split_train.tmp.bed \
		-s \
		-nameOnly \
		-tab |\
		sed 's/(-)//g' | sed 's/(+)//g' |\
		gzip >\
		GSE190214_${cell_type}.Enhancer_Promoter.chimeric_fragment.out.150nt.EP_anno.uniq.v3.split_train.fa.tab.gz
	#exit
	rm -f GSE190214_${cell_type}.Enhancer_Promoter.chimeric_fragment.out.150nt.EP_anno.uniq.v3.split_train.tmp.bed
	zcat GSE190214_${cell_type}.Enhancer_Promoter.chimeric_fragment.out.150nt.EP_anno.uniq.v3.split_train.fa.tab.gz |\
	awk '{OFS="\t"; gsub(/=/," ")}$2==1{id1[$1]=$3;getline;id2[$1]=$2}END{for (key in id1) {if (id1[key]!=""&&id2[key"=2"]!="") {print key,id1[key],id2[key"=2"]}}}' |\
		awk '{OFS="\t"}$2!~/N/&&$3!~/N/' |\
		gzip >\
		GSE190214_${cell_type}.Enhancer_Promoter.chimeric_fragment.out.150nt.EP_anno.uniq.v3.split_train.fa.pairs.tab.gz
#exit
	cat GSE190214_${cell_type}.Enhancer_Promoter.chimeric_fragment.out.150nt.EP_anno.uniq.v3.split_test.bedpe |\
		awk '{OFS="\t"}{gsub(/chr/,"");printf "%s\t%.0f\t%.0f\t%s\t%s\t%s\n%s\t%.0f\t%.0f\t%s\t%s\t%s\n", $1,$2,$3,$7"=1",".",$9,$4,$5,$6,$7"=2",".",$10}' >\
		GSE190214_${cell_type}.Enhancer_Promoter.chimeric_fragment.out.150nt.EP_anno.uniq.v3.split_test.tmp.bed

	bedtools getfasta \
		-fi $genome_fasta \
		-bed GSE190214_${cell_type}.Enhancer_Promoter.chimeric_fragment.out.150nt.EP_anno.uniq.v3.split_test.tmp.bed \
		-s \
		-nameOnly \
		-tab |\
		sed 's/(-)//g' | sed 's/(+)//g' |\
		gzip >\
		GSE190214_${cell_type}.Enhancer_Promoter.chimeric_fragment.out.150nt.EP_anno.uniq.v3.split_test.fa.tab.gz
	#exit
	rm -f GSE190214_${cell_type}.Enhancer_Promoter.chimeric_fragment.out.150nt.EP_anno.uniq.v3.split_test.tmp.bed
	zcat GSE190214_${cell_type}.Enhancer_Promoter.chimeric_fragment.out.150nt.EP_anno.uniq.v3.split_test.fa.tab.gz |\
	awk '{OFS="\t"; gsub(/=/," ")}$2==1{id1[$1]=$3;getline;id2[$1]=$2}END{for (key in id1) {if (id1[key]!=""&&id2[key"=2"]!="") {print key,id1[key],id2[key"=2"]}}}' |\
		awk '{OFS="\t"}$2!~/N/&&$3!~/N/' |\
		gzip >\
		GSE190214_${cell_type}.Enhancer_Promoter.chimeric_fragment.out.150nt.EP_anno.uniq.v3.split_test.fa.pairs.tab.gz
	cat GSE190214_${cell_type}.Enhancer_Promoter.chimeric_fragment.out.150nt.EP_anno.uniq.v3.split_validate.bedpe |\
		awk '{OFS="\t"}{gsub(/chr/,"");printf "%s\t%.0f\t%.0f\t%s\t%s\t%s\n%s\t%.0f\t%.0f\t%s\t%s\t%s\n", $1,$2,$3,$7"=1",".",$9,$4,$5,$6,$7"=2",".",$10}' >\
		GSE190214_${cell_type}.Enhancer_Promoter.chimeric_fragment.out.150nt.EP_anno.uniq.v3.split_validate.tmp.bed

	bedtools getfasta \
		-fi $genome_fasta \
		-bed GSE190214_${cell_type}.Enhancer_Promoter.chimeric_fragment.out.150nt.EP_anno.uniq.v3.split_validate.tmp.bed \
		-s \
		-nameOnly \
		-tab |\
		sed 's/(-)//g' | sed 's/(+)//g' |\
		gzip >\
		GSE190214_${cell_type}.Enhancer_Promoter.chimeric_fragment.out.150nt.EP_anno.uniq.v3.split_validate.fa.tab.gz
	#exit
	rm -f GSE190214_${cell_type}.Enhancer_Promoter.chimeric_fragment.out.150nt.EP_anno.uniq.v3.split_validate.tmp.bed
	zcat GSE190214_${cell_type}.Enhancer_Promoter.chimeric_fragment.out.150nt.EP_anno.uniq.v3.split_validate.fa.tab.gz |\
	awk '{OFS="\t"; gsub(/=/," ")}$2==1{id1[$1]=$3;getline;id2[$1]=$2}END{for (key in id1) {if (id1[key]!=""&&id2[key"=2"]!="") {print key,id1[key],id2[key"=2"]}}}' |\
		awk '{OFS="\t"}$2!~/N/&&$3!~/N/' |\
		gzip >\
		GSE190214_${cell_type}.Enhancer_Promoter.chimeric_fragment.out.150nt.EP_anno.uniq.v3.split_validate.fa.pairs.tab.gz

	## Generate unique and random regions
# Split Train	
	python3 generate_random_regions.py \
        	--bed Homo_sapiens.GRCh37.75.noncoding.regions.${cell_type}.v3.split_train.bed \
		--n_samples $(zcat GSE190214_${cell_type}.Enhancer_Promoter.chimeric_fragment.out.150nt.EP_anno.uniq.v3.split_train.fa.pairs.tab.gz | wc -l) \
		--width ${bin_size} \
        	--seed ${seed} \
        	--output ./negative_dataset/Homo_sapiens.GRCh37.75.noncoding.regions.${cell_type}.seed${seed}.${bin_size}nt.v3.split_train.bed
	sort -k1,1 -k2,2n ./negative_dataset/Homo_sapiens.GRCh37.75.noncoding.regions.${cell_type}.seed${seed}.${bin_size}nt.v3.split_train.bed >\
		./negative_dataset/Homo_sapiens.GRCh37.75.noncoding.regions.${cell_type}.seed${seed}.${bin_size}nt.v3.split_train.sort.bed
	bedtools getfasta \
                -fi $genome_fasta \
                -bed ./negative_dataset/Homo_sapiens.GRCh37.75.noncoding.regions.${cell_type}.seed${seed}.${bin_size}nt.v3.split_train.sort.bed \
                -s \
                -nameOnly \
                -tab |\
                sed 's/(-)//g' | sed 's/(+)//g' |\
                gzip >\
		./negative_dataset/Homo_sapiens.GRCh37.75.noncoding.regions.${cell_type}.seed${seed}.${bin_size}nt.v3.split_train.sort.fa.tab.gz
	# Enhancers constant
	paste -d"\t" \
		<(zcat GSE190214_${cell_type}.Enhancer_Promoter.chimeric_fragment.out.150nt.EP_anno.uniq.v3.split_train.fa.pairs.tab.gz) \
		<(zcat ./negative_dataset/Homo_sapiens.GRCh37.75.noncoding.regions.${cell_type}.seed${seed}.${bin_size}nt.v3.split_train.sort.fa.tab.gz) |\
		awk '{OFS="\t"}{print $1"_neg",$2,$5}' |\
		gzip >\
		./negative_dataset/GSE190214_${cell_type}.Enhancer_Promoter.chimeric_fragment.out.150nt.EP_anno.uniq.noncoding_seed${seed}.v3.split_train.enhancers.fa.pairs.tab.gz
	paste -d"\n" \
		<(zcat GSE190214_${cell_type}.Enhancer_Promoter.chimeric_fragment.out.150nt.EP_anno.uniq.v3.split_train.fa.pairs.tab.gz) \
		<(zcat ./negative_dataset/GSE190214_${cell_type}.Enhancer_Promoter.chimeric_fragment.out.150nt.EP_anno.uniq.noncoding_seed${seed}.v3.split_train.enhancers.fa.pairs.tab.gz) |\
		gzip >\
		${out_dir}/split_train/GSE190214_${cell_type}.Enhancer_Promoter.chimeric_fragment.out.150nt.EP_anno.uniq.noncoding_seed${seed}.v3.split_train.enhancers.fa.pairs.tab.gz
        
	zcat ${out_dir}/split_train/GSE190214_${cell_type}.Enhancer_Promoter.chimeric_fragment.out.150nt.EP_anno.uniq.noncoding_seed${seed}.v3.split_train.enhancers.fa.pairs.tab.gz |\
                awk '{OFS="\t"}{gsub(/U/,"T",$2);gsub(/U/,"T",$3); print $1,$2,$3}' | gzip >\
			${out_dir}/split_train/GSE190214_${cell_type}.Enhancer_Promoter.chimeric_fragment.out.150nt.EP_anno.uniq.noncoding_seed${seed}.v3.split_train.enhancers.fa.pairs.RNA.out.tsv.gz &

	# Promoters constant
	paste -d"\t" \
		<(zcat GSE190214_${cell_type}.Enhancer_Promoter.chimeric_fragment.out.150nt.EP_anno.uniq.v3.split_train.fa.pairs.tab.gz) \
		<(zcat ./negative_dataset/Homo_sapiens.GRCh37.75.noncoding.regions.${cell_type}.seed${seed}.${bin_size}nt.v3.split_train.sort.fa.tab.gz) |\
		awk '{OFS="\t"}{print $1"_neg",$3,$5}' |\
		gzip >\
		./negative_dataset/GSE190214_${cell_type}.Enhancer_Promoter.chimeric_fragment.out.150nt.EP_anno.uniq.noncoding_seed${seed}.v3.split_train.promoters.fa.pairs.tab.gz
	paste -d"\n" \
		<(zcat GSE190214_${cell_type}.Enhancer_Promoter.chimeric_fragment.out.150nt.EP_anno.uniq.v3.split_train.fa.pairs.tab.gz | awk '{OFS="\t"}{print $1,$3,$2}') \
		<(zcat ./negative_dataset/GSE190214_${cell_type}.Enhancer_Promoter.chimeric_fragment.out.150nt.EP_anno.uniq.noncoding_seed${seed}.v3.split_train.promoters.fa.pairs.tab.gz) |\
		gzip >\
		${out_dir}/split_train/GSE190214_${cell_type}.Enhancer_Promoter.chimeric_fragment.out.150nt.EP_anno.uniq.noncoding_seed${seed}.v3.split_train.promoters.fa.pairs.tab.gz

		zcat ${out_dir}/split_train/GSE190214_${cell_type}.Enhancer_Promoter.chimeric_fragment.out.150nt.EP_anno.uniq.noncoding_seed${seed}.v3.split_train.promoters.fa.pairs.tab.gz |\
			awk '{OFS="\t"}{gsub(/U/,"T",$2);gsub(/U/,"T",$3); print $1,$2,$3}' | gzip >\
			${out_dir}/split_train/GSE190214_${cell_type}.Enhancer_Promoter.chimeric_fragment.out.150nt.EP_anno.uniq.noncoding_seed${seed}.v3.split_train.promoters.fa.pairs.RNA.out.tsv.gz &
wait
       	
#exit
# Split Test
	python3 generate_random_regions.py \
        	--bed Homo_sapiens.GRCh37.75.noncoding.regions.${cell_type}.v3.split_test.bed \
        	--n_samples $(cat GSE190214_${cell_type}.Enhancer_Promoter.chimeric_fragment.out.150nt.EP_anno.uniq.v3.split_test.bedpe | wc -l) \
		--width ${bin_size} \
        	--seed ${seed} \
        	--output ./negative_dataset/Homo_sapiens.GRCh37.75.noncoding.regions.${cell_type}.seed${seed}.${bin_size}nt.v3.split_test.bed
	sort -k1,1 -k2,2n ./negative_dataset/Homo_sapiens.GRCh37.75.noncoding.regions.${cell_type}.seed${seed}.${bin_size}nt.v3.split_test.bed >\
		./negative_dataset/Homo_sapiens.GRCh37.75.noncoding.regions.${cell_type}.seed${seed}.${bin_size}nt.v3.split_test.sort.bed
	bedtools getfasta \
                -fi $genome_fasta \
                -bed ./negative_dataset/Homo_sapiens.GRCh37.75.noncoding.regions.${cell_type}.seed${seed}.${bin_size}nt.v3.split_test.sort.bed \
                -s \
                -nameOnly \
                -tab |\
                sed 's/(-)//g' | sed 's/(+)//g' |\
                gzip >\
		./negative_dataset/Homo_sapiens.GRCh37.75.noncoding.regions.${cell_type}.seed${seed}.${bin_size}nt.v3.split_test.sort.fa.tab.gz

	# Enhancers constant
	paste -d"\t" \
		<(zcat GSE190214_${cell_type}.Enhancer_Promoter.chimeric_fragment.out.150nt.EP_anno.uniq.v3.split_test.fa.pairs.tab.gz) \
		<(zcat ./negative_dataset/Homo_sapiens.GRCh37.75.noncoding.regions.${cell_type}.seed${seed}.${bin_size}nt.v3.split_test.sort.fa.tab.gz) |\
		awk '{OFS="\t"}{print $1"_neg",$2,$5}' |\
		gzip >\
		./negative_dataset/GSE190214_${cell_type}.Enhancer_Promoter.chimeric_fragment.out.150nt.EP_anno.uniq.noncoding_seed${seed}.v3.split_test.enhancers.fa.pairs.tab.gz
	paste -d"\n" \
		<(zcat GSE190214_${cell_type}.Enhancer_Promoter.chimeric_fragment.out.150nt.EP_anno.uniq.v3.split_test.fa.pairs.tab.gz) \
		<(zcat ./negative_dataset/GSE190214_${cell_type}.Enhancer_Promoter.chimeric_fragment.out.150nt.EP_anno.uniq.noncoding_seed${seed}.v3.split_test.enhancers.fa.pairs.tab.gz) |\
		gzip >\
		${out_dir}/split_test/GSE190214_${cell_type}.Enhancer_Promoter.chimeric_fragment.out.150nt.EP_anno.uniq.noncoding_seed${seed}.v3.split_test.enhancers.fa.pairs.tab.gz
	
	zcat ${out_dir}/split_test/GSE190214_${cell_type}.Enhancer_Promoter.chimeric_fragment.out.150nt.EP_anno.uniq.noncoding_seed${seed}.v3.split_test.enhancers.fa.pairs.tab.gz |\
		awk '{OFS="\t"}{gsub(/U/,"T",$2);gsub(/U/,"T",$3); print $1,$2,$3}' | gzip >\
			${out_dir}/split_test/GSE190214_${cell_type}.Enhancer_Promoter.chimeric_fragment.out.150nt.EP_anno.uniq.noncoding_seed${seed}.v3.split_test.enhancers.fa.pairs.RNA.out.tsv.gz &


	# Promoters constant
	paste -d"\t" \
		<(zcat GSE190214_${cell_type}.Enhancer_Promoter.chimeric_fragment.out.150nt.EP_anno.uniq.v3.split_test.fa.pairs.tab.gz) \
		<(zcat ./negative_dataset/Homo_sapiens.GRCh37.75.noncoding.regions.${cell_type}.seed${seed}.${bin_size}nt.v3.split_test.sort.fa.tab.gz) |\
		awk '{OFS="\t"}{print $1"_neg",$3,$5}' |\
		gzip >\
		./negative_dataset/GSE190214_${cell_type}.Enhancer_Promoter.chimeric_fragment.out.150nt.EP_anno.uniq.noncoding_seed${seed}.v3.split_test.promoters.fa.pairs.tab.gz
	paste -d"\n" \
		<(zcat GSE190214_${cell_type}.Enhancer_Promoter.chimeric_fragment.out.150nt.EP_anno.uniq.v3.split_test.fa.pairs.tab.gz | awk '{OFS="\t"}{print $1,$3,$2}') \
		<(zcat ./negative_dataset/GSE190214_${cell_type}.Enhancer_Promoter.chimeric_fragment.out.150nt.EP_anno.uniq.noncoding_seed${seed}.v3.split_test.promoters.fa.pairs.tab.gz) |\
		gzip >\
		${out_dir}/split_test/GSE190214_${cell_type}.Enhancer_Promoter.chimeric_fragment.out.150nt.EP_anno.uniq.noncoding_seed${seed}.v3.split_test.promoters.fa.pairs.tab.gz
		
	zcat ${out_dir}/split_test/GSE190214_${cell_type}.Enhancer_Promoter.chimeric_fragment.out.150nt.EP_anno.uniq.noncoding_seed${seed}.v3.split_test.promoters.fa.pairs.tab.gz |\
		awk '{OFS="\t"}{gsub(/U/,"T",$2);gsub(/U/,"T",$3); print $1,$2,$3}' | gzip >\
			${out_dir}/split_test/GSE190214_${cell_type}.Enhancer_Promoter.chimeric_fragment.out.150nt.EP_anno.uniq.noncoding_seed${seed}.v3.split_test.promoters.fa.pairs.RNA.out.tsv.gz &

wait
	# Split Validate
	python3 generate_random_regions.py \
        	--bed Homo_sapiens.GRCh37.75.noncoding.regions.${cell_type}.v3.split_validate.bed \
        	--n_samples $(cat GSE190214_${cell_type}.Enhancer_Promoter.chimeric_fragment.out.150nt.EP_anno.uniq.v3.split_validate.bedpe | wc -l) \
		--width ${bin_size} \
        	--seed ${seed} \
        	--output ./negative_dataset/Homo_sapiens.GRCh37.75.noncoding.regions.${cell_type}.seed${seed}.${bin_size}nt.v3.split_validate.bed
	sort -k1,1 -k2,2n ./negative_dataset/Homo_sapiens.GRCh37.75.noncoding.regions.${cell_type}.seed${seed}.${bin_size}nt.v3.split_validate.bed >\
		./negative_dataset/Homo_sapiens.GRCh37.75.noncoding.regions.${cell_type}.seed${seed}.${bin_size}nt.v3.split_validate.sort.bed
	bedtools getfasta \
                -fi $genome_fasta \
                -bed ./negative_dataset/Homo_sapiens.GRCh37.75.noncoding.regions.${cell_type}.seed${seed}.${bin_size}nt.v3.split_validate.sort.bed \
                -s \
                -nameOnly \
                -tab |\
                sed 's/(-)//g' | sed 's/(+)//g' |\
                gzip >\
		./negative_dataset/Homo_sapiens.GRCh37.75.noncoding.regions.${cell_type}.seed${seed}.${bin_size}nt.v3.split_validate.sort.fa.tab.gz
	# Enhancers constant
	paste -d"\t" \
		<(zcat GSE190214_${cell_type}.Enhancer_Promoter.chimeric_fragment.out.150nt.EP_anno.uniq.v3.split_validate.fa.pairs.tab.gz) \
		<(zcat ./negative_dataset/Homo_sapiens.GRCh37.75.noncoding.regions.${cell_type}.seed${seed}.${bin_size}nt.v3.split_validate.sort.fa.tab.gz) |\
		awk '{OFS="\t"}{print $1"_neg",$2,$5}' |\
		gzip >\
		./negative_dataset/GSE190214_${cell_type}.Enhancer_Promoter.chimeric_fragment.out.150nt.EP_anno.uniq.noncoding_seed${seed}.v3.split_validate.enhancers.fa.pairs.tab.gz
	paste -d"\n" \
		<(zcat GSE190214_${cell_type}.Enhancer_Promoter.chimeric_fragment.out.150nt.EP_anno.uniq.v3.split_validate.fa.pairs.tab.gz) \
		<(zcat ./negative_dataset/GSE190214_${cell_type}.Enhancer_Promoter.chimeric_fragment.out.150nt.EP_anno.uniq.noncoding_seed${seed}.v3.split_validate.enhancers.fa.pairs.tab.gz) |\
		gzip >\
		${out_dir}/split_validate/GSE190214_${cell_type}.Enhancer_Promoter.chimeric_fragment.out.150nt.EP_anno.uniq.noncoding_seed${seed}.v3.split_validate.enhancers.fa.pairs.tab.gz

	zcat ${out_dir}/split_validate/GSE190214_${cell_type}.Enhancer_Promoter.chimeric_fragment.out.150nt.EP_anno.uniq.noncoding_seed${seed}.v3.split_validate.enhancers.fa.pairs.tab.gz |\
		awk '{OFS="\t"}{gsub(/U/,"T",$2);gsub(/U/,"T",$3); print $1,$2,$3}' | gzip >\
			${out_dir}/split_validate/GSE190214_${cell_type}.Enhancer_Promoter.chimeric_fragment.out.150nt.EP_anno.uniq.noncoding_seed${seed}.v3.split_validate.enhancers.fa.pairs.RNA.out.tsv.gz &



	# Promoters constant
	paste -d"\t" \
		<(zcat GSE190214_${cell_type}.Enhancer_Promoter.chimeric_fragment.out.150nt.EP_anno.uniq.v3.split_validate.fa.pairs.tab.gz) \
		<(zcat ./negative_dataset/Homo_sapiens.GRCh37.75.noncoding.regions.${cell_type}.seed${seed}.${bin_size}nt.v3.split_validate.sort.fa.tab.gz) |\
		awk '{OFS="\t"}{print $1"_neg",$3,$5}' |\
		gzip >\
		./negative_dataset/GSE190214_${cell_type}.Enhancer_Promoter.chimeric_fragment.out.150nt.EP_anno.uniq.noncoding_seed${seed}.v3.split_validate.promoters.fa.pairs.tab.gz
	paste -d"\n" \
		<(zcat GSE190214_${cell_type}.Enhancer_Promoter.chimeric_fragment.out.150nt.EP_anno.uniq.v3.split_validate.fa.pairs.tab.gz | awk '{OFS="\t"}{print $1,$3,$2}') \
		<(zcat ./negative_dataset/GSE190214_${cell_type}.Enhancer_Promoter.chimeric_fragment.out.150nt.EP_anno.uniq.noncoding_seed${seed}.v3.split_validate.promoters.fa.pairs.tab.gz) |\
		gzip >\
		${out_dir}/split_validate/GSE190214_${cell_type}.Enhancer_Promoter.chimeric_fragment.out.150nt.EP_anno.uniq.noncoding_seed${seed}.v3.split_validate.promoters.fa.pairs.tab.gz
	
	zcat ${out_dir}/split_validate/GSE190214_${cell_type}.Enhancer_Promoter.chimeric_fragment.out.150nt.EP_anno.uniq.noncoding_seed${seed}.v3.split_validate.promoters.fa.pairs.tab.gz |\
		awk '{OFS="\t"}{gsub(/U/,"T",$2);gsub(/U/,"T",$3); print $1,$2,$3}' | gzip >\
			${out_dir}/split_validate/GSE190214_${cell_type}.Enhancer_Promoter.chimeric_fragment.out.150nt.EP_anno.uniq.noncoding_seed${seed}.v3.split_validate.promoters.fa.pairs.RNA.out.tsv.gz &


wait
	#exit
	#exit
done


