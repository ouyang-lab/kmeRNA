#!/bin/bash

# MIT License
# Copyright (c) 2025 Eric Nels Pederson, University of Massachusetts Amherst
# Licensed under the MIT License. See LICENSE file in the project root for full license information.

genome_faidx=/path/to/genome/Homo_sapiens.GRCh37.75.dna_sm.toplevel.fa.fai
genome_fasta=/path/to/genome/Homo_sapiens.GRCh37.75.dna_sm.toplevel.fa
GTF=/path/to/genome/Homo_sapiens.GRCh37.75.sort.out.gtf
seed=42813
bin_size=150
out_dir=/path/to/work_dir/
final_dir=${out_dir}/cross_cell/
mkdir -p $final_dir

for cell_type in HeLa hNPC GM12878
do
for cell_type2 in HeLa hNPC GM12878 IMR90 H1 HepG2 K562
do
if [[ $cell_type == $cell_type2 ]]
then
	continue
fi

echo "Trained on: "${cell_type}
echo "Testing on: "${cell_type2}
#echo "test"
awk '{OFS="\t"}FNR==NR{a[$1]+=1;next}{if (a[$1]>0&&a[$4]>0) {print $0}}' \
	<(awk '{print $1"\n"$2}' GSE190214_${cell_type}.Enhancer_Promoter.chimeric_fragment.out.150nt.EP_anno.uniq.v3.split_test.txt | sort  | uniq) \
	GSE190214_${cell_type2}.Enhancer_Promoter.chimeric_fragment.out.150nt.EP_anno.uniq.bedpe >\
	GSE190214_${cell_type}-${cell_type2}.Enhancer_Promoter.chimeric_fragment.out.150nt.EP_anno.uniq.v3.split_test.bedpe
	
	awk '{OFS="\t"}FNR==NR{a[$1]+=1;next}{if (a["chr"$1]>0) {print $0}}' \
		<(awk '{print $1"\n"$2}' GSE190214_${cell_type}.Enhancer_Promoter.chimeric_fragment.out.150nt.EP_anno.uniq.v3.split_test.txt | sort | uniq) \
		Homo_sapiens.GRCh37.75.noncoding.regions.bed >\
		Homo_sapiens.GRCh37.75.noncoding.regions.${cell_type}-${cell_type2}.v3.split_test.bed
	cat GSE190214_${cell_type}-${cell_type2}.Enhancer_Promoter.chimeric_fragment.out.150nt.EP_anno.uniq.v3.split_test.bedpe |\
		awk '{OFS="\t"}{gsub(/chr/,"");printf "%s\t%.0f\t%.0f\t%s\t%s\t%s\n%s\t%.0f\t%.0f\t%s\t%s\t%s\n", $1,$2,$3,$7"=1",".",$9,$4,$5,$6,$7"=2",".",$10}' >\
		GSE190214_${cell_type}-${cell_type2}.Enhancer_Promoter.chimeric_fragment.out.150nt.EP_anno.uniq.v3.split_test.tmp.bed
	#exit
	bedtools getfasta \
		-fi $genome_fasta \
		-bed GSE190214_${cell_type}-${cell_type2}.Enhancer_Promoter.chimeric_fragment.out.150nt.EP_anno.uniq.v3.split_test.tmp.bed \
		-s \
		-nameOnly \
		-tab |\
		sed 's/(-)//g' | sed 's/(+)//g' |\
		gzip >\
		GSE190214_${cell_type}-${cell_type2}.Enhancer_Promoter.chimeric_fragment.out.150nt.EP_anno.uniq.v3.split_test.fa.tab.gz
	#exit
	rm -f GSE190214_${cell_type}-${cell_type2}.Enhancer_Promoter.chimeric_fragment.out.150nt.EP_anno.uniq.v3.split_test.tmp.bed
	zcat GSE190214_${cell_type}-${cell_type2}.Enhancer_Promoter.chimeric_fragment.out.150nt.EP_anno.uniq.v3.split_test.fa.tab.gz |\
	awk '{OFS="\t"; gsub(/=/," ")}$2==1{id1[$1]=$3;getline;id2[$1]=$2}END{for (key in id1) {if (id1[key]!=""&&id2[key"=2"]!="") {print key,id1[key],id2[key"=2"]}}}' |\
		awk '{OFS="\t"}$2!~/N/&&$3!~/N/' |\
		gzip >\
		GSE190214_${cell_type}-${cell_type2}.Enhancer_Promoter.chimeric_fragment.out.150nt.EP_anno.uniq.v3.split_test.fa.pairs.tab.gz
	#exit
	# Split Test
	python3 generate_random_regions.py \
        	--bed Homo_sapiens.GRCh37.75.noncoding.regions.${cell_type}-${cell_type2}.v3.split_test.bed \
        	--n_samples $(cat GSE190214_${cell_type}-${cell_type2}.Enhancer_Promoter.chimeric_fragment.out.150nt.EP_anno.uniq.v3.split_test.bedpe | wc -l) \
		--width ${bin_size} \
        	--seed ${seed} \
        	--output ./negative_dataset/Homo_sapiens.GRCh37.75.noncoding.regions.${cell_type}-${cell_type2}.seed${seed}.${bin_size}nt.v3.split_test.bed
	sort -k1,1 -k2,2n ./negative_dataset/Homo_sapiens.GRCh37.75.noncoding.regions.${cell_type}-${cell_type2}.seed${seed}.${bin_size}nt.v3.split_test.bed >\
		./negative_dataset/Homo_sapiens.GRCh37.75.noncoding.regions.${cell_type}-${cell_type2}.seed${seed}.${bin_size}nt.v3.split_test.sort.bed
	#exit
	bedtools getfasta \
                -fi $genome_fasta \
                -bed ./negative_dataset/Homo_sapiens.GRCh37.75.noncoding.regions.${cell_type}-${cell_type2}.seed${seed}.${bin_size}nt.v3.split_test.sort.bed \
                -s \
                -nameOnly \
                -tab |\
                sed 's/(-)//g' | sed 's/(+)//g' |\
                gzip >\
		./negative_dataset/Homo_sapiens.GRCh37.75.noncoding.regions.${cell_type}-${cell_type2}.seed${seed}.${bin_size}nt.v3.split_test.sort.fa.tab.gz
#exit
	# Enhancers constant
	paste -d"\t" \
		<(zcat GSE190214_${cell_type}-${cell_type2}.Enhancer_Promoter.chimeric_fragment.out.150nt.EP_anno.uniq.v3.split_test.fa.pairs.tab.gz) \
		<(zcat ./negative_dataset/Homo_sapiens.GRCh37.75.noncoding.regions.${cell_type}-${cell_type2}.seed${seed}.${bin_size}nt.v3.split_test.sort.fa.tab.gz) |\
		awk '{OFS="\t"}{print $1"_neg",$2,$5}' |\
		gzip >\
		./negative_dataset/GSE190214_${cell_type}-${cell_type2}.Enhancer_Promoter.chimeric_fragment.out.150nt.EP_anno.uniq.noncoding_seed${seed}.v3.split_test.enhancers.fa.pairs.tab.gz
	paste -d"\n" \
		<(zcat GSE190214_${cell_type}-${cell_type2}.Enhancer_Promoter.chimeric_fragment.out.150nt.EP_anno.uniq.v3.split_test.fa.pairs.tab.gz) \
		<(zcat ./negative_dataset/GSE190214_${cell_type}-${cell_type2}.Enhancer_Promoter.chimeric_fragment.out.150nt.EP_anno.uniq.noncoding_seed${seed}.v3.split_test.enhancers.fa.pairs.tab.gz) |\
		gzip >\
		${out_dir}/split_test/GSE190214_${cell_type}-${cell_type2}.Enhancer_Promoter.chimeric_fragment.out.150nt.EP_anno.uniq.noncoding_seed${seed}.v3.split_test.enhancers.fa.pairs.tab.gz
	#exit	
	zcat ${out_dir}/split_test/GSE190214_${cell_type}-${cell_type2}.Enhancer_Promoter.chimeric_fragment.out.150nt.EP_anno.uniq.noncoding_seed${seed}.v3.split_test.enhancers.fa.pairs.tab.gz |\
		awk '{OFS="\t"}{gsub(/U/,"T",$2);gsub(/U/,"T",$3); print $1,$2,$3}' | gzip >\
			${out_dir}/split_test/GSE190214_${cell_type}-${cell_type2}.Enhancer_Promoter.chimeric_fragment.out.150nt.EP_anno.uniq.noncoding_seed${seed}.v3.split_test.enhancers.fa.pairs.RNA.out.tsv.gz &

	# Promoters constant
	paste -d"\t" \
		<(zcat GSE190214_${cell_type}-${cell_type2}.Enhancer_Promoter.chimeric_fragment.out.150nt.EP_anno.uniq.v3.split_test.fa.pairs.tab.gz) \
		<(zcat ./negative_dataset/Homo_sapiens.GRCh37.75.noncoding.regions.${cell_type}-${cell_type2}.seed${seed}.${bin_size}nt.v3.split_test.sort.fa.tab.gz) |\
		awk '{OFS="\t"}{print $1"_neg",$3,$5}' |\
		gzip >\
		./negative_dataset/GSE190214_${cell_type}-${cell_type2}.Enhancer_Promoter.chimeric_fragment.out.150nt.EP_anno.uniq.noncoding_seed${seed}.v3.split_test.promoters.fa.pairs.tab.gz
	paste -d"\n" \
		<(zcat GSE190214_${cell_type}-${cell_type2}.Enhancer_Promoter.chimeric_fragment.out.150nt.EP_anno.uniq.v3.split_test.fa.pairs.tab.gz | awk '{OFS="\t"}{print $1,$3,$2}') \
		<(zcat ./negative_dataset/GSE190214_${cell_type}-${cell_type2}.Enhancer_Promoter.chimeric_fragment.out.150nt.EP_anno.uniq.noncoding_seed${seed}.v3.split_test.promoters.fa.pairs.tab.gz) |\
		gzip >\
		${out_dir}/split_test/GSE190214_${cell_type}-${cell_type2}.Enhancer_Promoter.chimeric_fragment.out.150nt.EP_anno.uniq.noncoding_seed${seed}.v3.split_test.promoters.fa.pairs.tab.gz
		
	zcat ${out_dir}/split_test/GSE190214_${cell_type}-${cell_type2}.Enhancer_Promoter.chimeric_fragment.out.150nt.EP_anno.uniq.noncoding_seed${seed}.v3.split_test.promoters.fa.pairs.tab.gz |\
		awk '{OFS="\t"}{gsub(/U/,"T",$2);gsub(/U/,"T",$3); print $1,$2,$3}' | gzip >\
			${out_dir}/split_test/GSE190214_${cell_type}-${cell_type2}.Enhancer_Promoter.chimeric_fragment.out.150nt.EP_anno.uniq.noncoding_seed${seed}.v3.split_test.promoters.fa.pairs.RNA.out.tsv.gz &

wait
#exit
mv ${out_dir}/split_test/GSE190214_${cell_type}-${cell_type2}.Enhancer_Promoter.chimeric_fragment.out.150nt.EP_anno.uniq.noncoding_seed${seed}.v3.split_test.enhancers.fa.pairs.RNA.out.tsv.gz $final_dir
mv ${out_dir}/split_test/GSE190214_${cell_type}-${cell_type2}.Enhancer_Promoter.chimeric_fragment.out.150nt.EP_anno.uniq.noncoding_seed${seed}.v3.split_test.promoters.fa.pairs.RNA.out.tsv.gz $final_dir
done
done



