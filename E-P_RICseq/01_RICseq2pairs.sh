#!/bin/bash
# MIT License
# Copyright (c) 2026 Eric Nels Pederson, University of Massachusetts Amherst
# Licensed under the MIT License. See LICENSE file in the project root for full license information.

genome_fasta=/path/to/genome/Homo_sapiens.GRCh37.75.dna_sm.toplevel.fa
out_bin_size=150

for cell in K562 IMR90 HeLa hNPC HepG2 GM12878 H1; 
do 
	awk '{OFS="\t"}NR>1{gsub(/,/," ");gsub(/:/," ");gsub(/-/," ");print $2,$3,$4,$6,$7,$8,$1";"$5,"."}' ${cell}_enhancer-promoter.csv | awk 'NF==8' |\
		sed 's/\r//g' >\
		${cell}_enhancer-promoter.bedpe
	awk '{OFS="\t"}{gsub(/chr/,"");gsub(/;/,"\t"); print $1,$2,$3,$7"\n"$4,$5,$6,$8}' \
		${cell}_enhancer-promoter.bedpe |\
		sort -k1,1 -k2,2g | uniq >\
		${cell}_enhancer-promoter.regions.bed

	./bed_bin.py \
		--input_bed GSE190214_${cell}.Enhancer_Promoter.chimeric_fragment.out.bed \
		--output_bed GSE190214_${cell}.Enhancer_Promoter.chimeric_fragment.out.${out_bin_size}nt.bed \
		--bin_size ${out_bin_size}
	bedtools intersect \
		-a <(sort -k1,1 -k2,2g GSE190214_${cell}.Enhancer_Promoter.chimeric_fragment.out.${out_bin_size}nt.bed) \
		-b ${cell}_enhancer-promoter.regions.bed \
		-wa -wb |\
		awk '{OFS="\t"}{print $1,$2,$3,$4,$5,$6,$10}' |\
		sort -k1,1 -k2,2g | uniq | sort -k4,4 >\
		GSE190214_${cell}.Enhancer_Promoter.chimeric_fragment.out.${out_bin_size}nt.EP_anno.bed
	awk '{OFS="\t"}/=1/{region="chr"$1"\t"$2"\t"$3;id=substr($4,1,length($4)-2);strand=$6;comment=$7; getline; print region,"chr"$1,$2,$3,id,".",strand,$6,comment,$7}' \
		GSE190214_${cell}.Enhancer_Promoter.chimeric_fragment.out.${out_bin_size}nt.EP_anno.bed |\
		awk '{OFS="\t"}{if (index($12,"MPT")>0) {print $0} else if (index($11,"MPT")>0) {print $4,$5,$6,$1,$2,$3,$7,$8,$10,$9,$12,$11}}' >\
		GSE190214_${cell}.Enhancer_Promoter.chimeric_fragment.out.${out_bin_size}nt.EP_anno.bedpe
	awk '{OFS="\t"}{$7=".";print $0}' \
		GSE190214_${cell}.Enhancer_Promoter.chimeric_fragment.out.${out_bin_size}nt.EP_anno.bedpe |\
		sort -k1,1 -k2,2g | uniq -c | awk '{OFS="\t"}{$8="bin_pair_"NR;print $2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$1}' >\
		GSE190214_${cell}.Enhancer_Promoter.chimeric_fragment.out.${out_bin_size}nt.EP_anno.uniq.bedpe

done
