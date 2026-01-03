#!/bin/bash

# MIT License
# Copyright (c) 2026 Eric Nels Pederson, University of Massachusetts Amherst
# Licensed under the MIT License. See LICENSE file in the project root for full license information.

# K562 IMR90 HeLa hNPC HepG2 GM12878 H1
cell_type=HeLa
awk '{OFS="\t"}{gsub(/chr/,"");if ($1<$4) {print $1,$4} else {print $4,$1}}' \
	GSE190214_${cell_type}.Enhancer_Promoter.chimeric_fragment.out.150nt.EP_anno.uniq.bedpe |\
	sort  | uniq -c | sort -k1,1gr | awk '{OFS="\t"}{print "chr"$2,"chr"$3,$1/95919*100,$1}' >\
	GSE190214_${cell_type}.Enhancer_Promoter.chimeric_fragment.out.150nt.EP_anno.uniq.pct.tsv

python3 ./split_dat.py \
	--input_file GSE190214_${cell_type}.Enhancer_Promoter.chimeric_fragment.out.150nt.EP_anno.uniq.pct.tsv \
	--output_prefix GSE190214_${cell_type}.Enhancer_Promoter.chimeric_fragment.out.150nt.EP_anno.uniq.v3.split \
	--seed 182
echo ${cell_type}
tot=$(grep -f <(cat GSE190214_${cell_type}.Enhancer_Promoter.chimeric_fragment.out.150nt.EP_anno.uniq.v3.split_train.txt \
	GSE190214_${cell_type}.Enhancer_Promoter.chimeric_fragment.out.150nt.EP_anno.uniq.v3.split_test.txt \
	GSE190214_${cell_type}.Enhancer_Promoter.chimeric_fragment.out.150nt.EP_anno.uniq.v3.split_validate.txt) GSE190214_${cell_type}.Enhancer_Promoter.chimeric_fragment.out.150nt.EP_anno.uniq.pct.tsv | awk '{sum+=$4}END{print sum}')
echo "train"
grep -f GSE190214_${cell_type}.Enhancer_Promoter.chimeric_fragment.out.150nt.EP_anno.uniq.v3.split_train.txt \
	GSE190214_${cell_type}.Enhancer_Promoter.chimeric_fragment.out.150nt.EP_anno.uniq.pct.tsv |\
	awk -v tot=$tot '{sum2+=$4}END{print sum2/tot*100"%"}'
awk '{OFS="\t"}FNR==NR{a[$1"\t"$2]+=1;next}{if (a[$1"\t"$4]>0||a[$4"\t"$1]>0) {print $0}}' \
	GSE190214_${cell_type}.Enhancer_Promoter.chimeric_fragment.out.150nt.EP_anno.uniq.v3.split_train.txt \
	GSE190214_${cell_type}.Enhancer_Promoter.chimeric_fragment.out.150nt.EP_anno.uniq.bedpe >\
	GSE190214_${cell_type}.Enhancer_Promoter.chimeric_fragment.out.150nt.EP_anno.uniq.v3.split_train.bedpe
echo "test"
grep -f GSE190214_${cell_type}.Enhancer_Promoter.chimeric_fragment.out.150nt.EP_anno.uniq.v3.split_test.txt \
	GSE190214_${cell_type}.Enhancer_Promoter.chimeric_fragment.out.150nt.EP_anno.uniq.pct.tsv |\
	awk -v tot=$tot '{sum2+=$4}END{print sum2/tot*100"%"}'
awk '{OFS="\t"}FNR==NR{a[$1"\t"$2]+=1;next}{if (a[$1"\t"$4]>0||a[$4"\t"$1]>0) {print $0}}' \
	GSE190214_${cell_type}.Enhancer_Promoter.chimeric_fragment.out.150nt.EP_anno.uniq.v3.split_test.txt \
	GSE190214_${cell_type}.Enhancer_Promoter.chimeric_fragment.out.150nt.EP_anno.uniq.bedpe >\
	GSE190214_${cell_type}.Enhancer_Promoter.chimeric_fragment.out.150nt.EP_anno.uniq.v3.split_test.bedpe

echo "validate"
grep -f GSE190214_${cell_type}.Enhancer_Promoter.chimeric_fragment.out.150nt.EP_anno.uniq.v3.split_validate.txt \
	GSE190214_${cell_type}.Enhancer_Promoter.chimeric_fragment.out.150nt.EP_anno.uniq.pct.tsv |\
	awk -v tot=$tot '{sum2+=$4}END{print sum2/tot*100"%"}'
awk '{OFS="\t"}FNR==NR{a[$1"\t"$2]+=1;next}{if (a[$1"\t"$4]>0||a[$4"\t"$1]>0) {print $0}}' \
	GSE190214_${cell_type}.Enhancer_Promoter.chimeric_fragment.out.150nt.EP_anno.uniq.v3.split_validate.txt \
	GSE190214_${cell_type}.Enhancer_Promoter.chimeric_fragment.out.150nt.EP_anno.uniq.bedpe >\
	GSE190214_${cell_type}.Enhancer_Promoter.chimeric_fragment.out.150nt.EP_anno.uniq.v3.split_validate.bedpe

#exit

cell_type=hNPC
awk '{OFS="\t"}{gsub(/chr/,"");if ($1<$4) {print $1,$4} else {print $4,$1}}' \
	GSE190214_${cell_type}.Enhancer_Promoter.chimeric_fragment.out.150nt.EP_anno.uniq.bedpe |\
	sort  | uniq -c | sort -k1,1gr | awk '{OFS="\t"}{print "chr"$2,"chr"$3,$1/95919*100,$1}' >\
	GSE190214_${cell_type}.Enhancer_Promoter.chimeric_fragment.out.150nt.EP_anno.uniq.pct.tsv

python3 ./split_dat.py \
	--input_file GSE190214_${cell_type}.Enhancer_Promoter.chimeric_fragment.out.150nt.EP_anno.uniq.pct.tsv \
	--output_prefix GSE190214_${cell_type}.Enhancer_Promoter.chimeric_fragment.out.150nt.EP_anno.uniq.v3.split \
	--seed 363
echo ${cell_type}
tot=$(grep -f <(cat GSE190214_${cell_type}.Enhancer_Promoter.chimeric_fragment.out.150nt.EP_anno.uniq.v3.split_train.txt \
	GSE190214_${cell_type}.Enhancer_Promoter.chimeric_fragment.out.150nt.EP_anno.uniq.v3.split_test.txt \
	GSE190214_${cell_type}.Enhancer_Promoter.chimeric_fragment.out.150nt.EP_anno.uniq.v3.split_validate.txt) GSE190214_${cell_type}.Enhancer_Promoter.chimeric_fragment.out.150nt.EP_anno.uniq.pct.tsv | awk '{sum+=$4}END{print sum}')
echo "train"
grep -f GSE190214_${cell_type}.Enhancer_Promoter.chimeric_fragment.out.150nt.EP_anno.uniq.v3.split_train.txt \
	GSE190214_${cell_type}.Enhancer_Promoter.chimeric_fragment.out.150nt.EP_anno.uniq.pct.tsv |\
	awk -v tot=$tot '{sum2+=$4}END{print sum2/tot*100"%"}'
awk '{OFS="\t"}FNR==NR{a[$1"\t"$2]+=1;next}{if (a[$1"\t"$4]>0||a[$4"\t"$1]>0) {print $0}}' \
	GSE190214_${cell_type}.Enhancer_Promoter.chimeric_fragment.out.150nt.EP_anno.uniq.v3.split_train.txt \
	GSE190214_${cell_type}.Enhancer_Promoter.chimeric_fragment.out.150nt.EP_anno.uniq.bedpe >\
	GSE190214_${cell_type}.Enhancer_Promoter.chimeric_fragment.out.150nt.EP_anno.uniq.v3.split_train.bedpe
echo "test"
grep -f GSE190214_${cell_type}.Enhancer_Promoter.chimeric_fragment.out.150nt.EP_anno.uniq.v3.split_test.txt \
	GSE190214_${cell_type}.Enhancer_Promoter.chimeric_fragment.out.150nt.EP_anno.uniq.pct.tsv |\
	awk -v tot=$tot '{sum2+=$4}END{print sum2/tot*100"%"}'
awk '{OFS="\t"}FNR==NR{a[$1"\t"$2]+=1;next}{if (a[$1"\t"$4]>0||a[$4"\t"$1]>0) {print $0}}' \
	GSE190214_${cell_type}.Enhancer_Promoter.chimeric_fragment.out.150nt.EP_anno.uniq.v3.split_test.txt \
	GSE190214_${cell_type}.Enhancer_Promoter.chimeric_fragment.out.150nt.EP_anno.uniq.bedpe >\
	GSE190214_${cell_type}.Enhancer_Promoter.chimeric_fragment.out.150nt.EP_anno.uniq.v3.split_test.bedpe

echo "validate"
grep -f GSE190214_${cell_type}.Enhancer_Promoter.chimeric_fragment.out.150nt.EP_anno.uniq.v3.split_validate.txt \
	GSE190214_${cell_type}.Enhancer_Promoter.chimeric_fragment.out.150nt.EP_anno.uniq.pct.tsv |\
	awk -v tot=$tot '{sum2+=$4}END{print sum2/tot*100"%"}'
awk '{OFS="\t"}FNR==NR{a[$1"\t"$2]+=1;next}{if (a[$1"\t"$4]>0||a[$4"\t"$1]>0) {print $0}}' \
	GSE190214_${cell_type}.Enhancer_Promoter.chimeric_fragment.out.150nt.EP_anno.uniq.v3.split_validate.txt \
	GSE190214_${cell_type}.Enhancer_Promoter.chimeric_fragment.out.150nt.EP_anno.uniq.bedpe >\
	GSE190214_${cell_type}.Enhancer_Promoter.chimeric_fragment.out.150nt.EP_anno.uniq.v3.split_validate.bedpe
#exit

cell_type=IMR90
awk '{OFS="\t"}{gsub(/chr/,"");if ($1<$4) {print $1,$4} else {print $4,$1}}' \
	GSE190214_${cell_type}.Enhancer_Promoter.chimeric_fragment.out.150nt.EP_anno.uniq.bedpe |\
	sort  | uniq -c | sort -k1,1gr | awk '{OFS="\t"}{print "chr"$2,"chr"$3,$1/95919*100,$1}' >\
	GSE190214_${cell_type}.Enhancer_Promoter.chimeric_fragment.out.150nt.EP_anno.uniq.pct.tsv

python3 ./split_dat.py \
	--input_file GSE190214_${cell_type}.Enhancer_Promoter.chimeric_fragment.out.150nt.EP_anno.uniq.pct.tsv \
	--output_prefix GSE190214_${cell_type}.Enhancer_Promoter.chimeric_fragment.out.150nt.EP_anno.uniq.v3.split \
	--seed 614
echo ${cell_type}
tot=$(grep -f <(cat GSE190214_${cell_type}.Enhancer_Promoter.chimeric_fragment.out.150nt.EP_anno.uniq.v3.split_train.txt \
	GSE190214_${cell_type}.Enhancer_Promoter.chimeric_fragment.out.150nt.EP_anno.uniq.v3.split_test.txt \
	GSE190214_${cell_type}.Enhancer_Promoter.chimeric_fragment.out.150nt.EP_anno.uniq.v3.split_validate.txt) GSE190214_${cell_type}.Enhancer_Promoter.chimeric_fragment.out.150nt.EP_anno.uniq.pct.tsv | awk '{sum+=$4}END{print sum}')
echo "train"
grep -f GSE190214_${cell_type}.Enhancer_Promoter.chimeric_fragment.out.150nt.EP_anno.uniq.v3.split_train.txt \
	GSE190214_${cell_type}.Enhancer_Promoter.chimeric_fragment.out.150nt.EP_anno.uniq.pct.tsv |\
	awk -v tot=$tot '{sum2+=$4}END{print sum2/tot*100"%"}'
awk '{OFS="\t"}FNR==NR{a[$1"\t"$2]+=1;next}{if (a[$1"\t"$4]>0||a[$4"\t"$1]>0) {print $0}}' \
	GSE190214_${cell_type}.Enhancer_Promoter.chimeric_fragment.out.150nt.EP_anno.uniq.v3.split_train.txt \
	GSE190214_${cell_type}.Enhancer_Promoter.chimeric_fragment.out.150nt.EP_anno.uniq.bedpe >\
	GSE190214_${cell_type}.Enhancer_Promoter.chimeric_fragment.out.150nt.EP_anno.uniq.v3.split_train.bedpe
echo "test"
grep -f GSE190214_${cell_type}.Enhancer_Promoter.chimeric_fragment.out.150nt.EP_anno.uniq.v3.split_test.txt \
	GSE190214_${cell_type}.Enhancer_Promoter.chimeric_fragment.out.150nt.EP_anno.uniq.pct.tsv |\
	awk -v tot=$tot '{sum2+=$4}END{print sum2/tot*100"%"}'
awk '{OFS="\t"}FNR==NR{a[$1"\t"$2]+=1;next}{if (a[$1"\t"$4]>0||a[$4"\t"$1]>0) {print $0}}' \
	GSE190214_${cell_type}.Enhancer_Promoter.chimeric_fragment.out.150nt.EP_anno.uniq.v3.split_test.txt \
	GSE190214_${cell_type}.Enhancer_Promoter.chimeric_fragment.out.150nt.EP_anno.uniq.bedpe >\
	GSE190214_${cell_type}.Enhancer_Promoter.chimeric_fragment.out.150nt.EP_anno.uniq.v3.split_test.bedpe

echo "validate"
grep -f GSE190214_${cell_type}.Enhancer_Promoter.chimeric_fragment.out.150nt.EP_anno.uniq.v3.split_validate.txt \
	GSE190214_${cell_type}.Enhancer_Promoter.chimeric_fragment.out.150nt.EP_anno.uniq.pct.tsv |\
	awk -v tot=$tot '{sum2+=$4}END{print sum2/tot*100"%"}'
awk '{OFS="\t"}FNR==NR{a[$1"\t"$2]+=1;next}{if (a[$1"\t"$4]>0||a[$4"\t"$1]>0) {print $0}}' \
	GSE190214_${cell_type}.Enhancer_Promoter.chimeric_fragment.out.150nt.EP_anno.uniq.v3.split_validate.txt \
	GSE190214_${cell_type}.Enhancer_Promoter.chimeric_fragment.out.150nt.EP_anno.uniq.bedpe >\
	GSE190214_${cell_type}.Enhancer_Promoter.chimeric_fragment.out.150nt.EP_anno.uniq.v3.split_validate.bedpe


#exit

cell_type=HepG2
awk '{OFS="\t"}{gsub(/chr/,"");if ($1<$4) {print $1,$4} else {print $4,$1}}' \
	GSE190214_${cell_type}.Enhancer_Promoter.chimeric_fragment.out.150nt.EP_anno.uniq.bedpe |\
	sort  | uniq -c | sort -k1,1gr | awk '{OFS="\t"}{print "chr"$2,"chr"$3,$1/95919*100,$1}' >\
	GSE190214_${cell_type}.Enhancer_Promoter.chimeric_fragment.out.150nt.EP_anno.uniq.pct.tsv

python3 ./split_dat.py \
	--input_file GSE190214_${cell_type}.Enhancer_Promoter.chimeric_fragment.out.150nt.EP_anno.uniq.pct.tsv \
	--output_prefix GSE190214_${cell_type}.Enhancer_Promoter.chimeric_fragment.out.150nt.EP_anno.uniq.v3.split \
	--seed 290
echo ${cell_type}
tot=$(grep -f <(cat GSE190214_${cell_type}.Enhancer_Promoter.chimeric_fragment.out.150nt.EP_anno.uniq.v3.split_train.txt \
	GSE190214_${cell_type}.Enhancer_Promoter.chimeric_fragment.out.150nt.EP_anno.uniq.v3.split_test.txt \
	GSE190214_${cell_type}.Enhancer_Promoter.chimeric_fragment.out.150nt.EP_anno.uniq.v3.split_validate.txt) GSE190214_${cell_type}.Enhancer_Promoter.chimeric_fragment.out.150nt.EP_anno.uniq.pct.tsv | awk '{sum+=$4}END{print sum}')
echo "train"
grep -f GSE190214_${cell_type}.Enhancer_Promoter.chimeric_fragment.out.150nt.EP_anno.uniq.v3.split_train.txt \
	GSE190214_${cell_type}.Enhancer_Promoter.chimeric_fragment.out.150nt.EP_anno.uniq.pct.tsv |\
	awk -v tot=$tot '{sum2+=$4}END{print sum2/tot*100"%"}'
awk '{OFS="\t"}FNR==NR{a[$1"\t"$2]+=1;next}{if (a[$1"\t"$4]>0||a[$4"\t"$1]>0) {print $0}}' \
	GSE190214_${cell_type}.Enhancer_Promoter.chimeric_fragment.out.150nt.EP_anno.uniq.v3.split_train.txt \
	GSE190214_${cell_type}.Enhancer_Promoter.chimeric_fragment.out.150nt.EP_anno.uniq.bedpe >\
	GSE190214_${cell_type}.Enhancer_Promoter.chimeric_fragment.out.150nt.EP_anno.uniq.v3.split_train.bedpe
echo "test"
grep -f GSE190214_${cell_type}.Enhancer_Promoter.chimeric_fragment.out.150nt.EP_anno.uniq.v3.split_test.txt \
	GSE190214_${cell_type}.Enhancer_Promoter.chimeric_fragment.out.150nt.EP_anno.uniq.pct.tsv |\
	awk -v tot=$tot '{sum2+=$4}END{print sum2/tot*100"%"}'
awk '{OFS="\t"}FNR==NR{a[$1"\t"$2]+=1;next}{if (a[$1"\t"$4]>0||a[$4"\t"$1]>0) {print $0}}' \
	GSE190214_${cell_type}.Enhancer_Promoter.chimeric_fragment.out.150nt.EP_anno.uniq.v3.split_test.txt \
	GSE190214_${cell_type}.Enhancer_Promoter.chimeric_fragment.out.150nt.EP_anno.uniq.bedpe >\
	GSE190214_${cell_type}.Enhancer_Promoter.chimeric_fragment.out.150nt.EP_anno.uniq.v3.split_test.bedpe

echo "validate"
grep -f GSE190214_${cell_type}.Enhancer_Promoter.chimeric_fragment.out.150nt.EP_anno.uniq.v3.split_validate.txt \
	GSE190214_${cell_type}.Enhancer_Promoter.chimeric_fragment.out.150nt.EP_anno.uniq.pct.tsv |\
	awk -v tot=$tot '{sum2+=$4}END{print sum2/tot*100"%"}'
awk '{OFS="\t"}FNR==NR{a[$1"\t"$2]+=1;next}{if (a[$1"\t"$4]>0||a[$4"\t"$1]>0) {print $0}}' \
	GSE190214_${cell_type}.Enhancer_Promoter.chimeric_fragment.out.150nt.EP_anno.uniq.v3.split_validate.txt \
	GSE190214_${cell_type}.Enhancer_Promoter.chimeric_fragment.out.150nt.EP_anno.uniq.bedpe >\
	GSE190214_${cell_type}.Enhancer_Promoter.chimeric_fragment.out.150nt.EP_anno.uniq.v3.split_validate.bedpe
#exit


cell_type=GM12878
awk '{OFS="\t"}{gsub(/chr/,"");if ($1<$4) {print $1,$4} else {print $4,$1}}' \
	GSE190214_${cell_type}.Enhancer_Promoter.chimeric_fragment.out.150nt.EP_anno.uniq.bedpe |\
	sort  | uniq -c | sort -k1,1gr | awk '{OFS="\t"}{print "chr"$2,"chr"$3,$1/95919*100,$1}' >\
	GSE190214_${cell_type}.Enhancer_Promoter.chimeric_fragment.out.150nt.EP_anno.uniq.pct.tsv

python3 ./split_dat.py \
	--input_file GSE190214_${cell_type}.Enhancer_Promoter.chimeric_fragment.out.150nt.EP_anno.uniq.pct.tsv \
	--output_prefix GSE190214_${cell_type}.Enhancer_Promoter.chimeric_fragment.out.150nt.EP_anno.uniq.v3.split \
	--seed 2827
echo ${cell_type}
tot=$(grep -f <(cat GSE190214_${cell_type}.Enhancer_Promoter.chimeric_fragment.out.150nt.EP_anno.uniq.v3.split_train.txt \
	GSE190214_${cell_type}.Enhancer_Promoter.chimeric_fragment.out.150nt.EP_anno.uniq.v3.split_test.txt \
	GSE190214_${cell_type}.Enhancer_Promoter.chimeric_fragment.out.150nt.EP_anno.uniq.v3.split_validate.txt) GSE190214_${cell_type}.Enhancer_Promoter.chimeric_fragment.out.150nt.EP_anno.uniq.pct.tsv | awk '{sum+=$4}END{print sum}')
echo "train"
grep -f GSE190214_${cell_type}.Enhancer_Promoter.chimeric_fragment.out.150nt.EP_anno.uniq.v3.split_train.txt \
	GSE190214_${cell_type}.Enhancer_Promoter.chimeric_fragment.out.150nt.EP_anno.uniq.pct.tsv |\
	awk -v tot=$tot '{sum2+=$4}END{print sum2/tot*100"%"}'
awk '{OFS="\t"}FNR==NR{a[$1"\t"$2]+=1;next}{if (a[$1"\t"$4]>0||a[$4"\t"$1]>0) {print $0}}' \
	GSE190214_${cell_type}.Enhancer_Promoter.chimeric_fragment.out.150nt.EP_anno.uniq.v3.split_train.txt \
	GSE190214_${cell_type}.Enhancer_Promoter.chimeric_fragment.out.150nt.EP_anno.uniq.bedpe >\
	GSE190214_${cell_type}.Enhancer_Promoter.chimeric_fragment.out.150nt.EP_anno.uniq.v3.split_train.bedpe
echo "test"
grep -f GSE190214_${cell_type}.Enhancer_Promoter.chimeric_fragment.out.150nt.EP_anno.uniq.v3.split_test.txt \
	GSE190214_${cell_type}.Enhancer_Promoter.chimeric_fragment.out.150nt.EP_anno.uniq.pct.tsv |\
	awk -v tot=$tot '{sum2+=$4}END{print sum2/tot*100"%"}'
awk '{OFS="\t"}FNR==NR{a[$1"\t"$2]+=1;next}{if (a[$1"\t"$4]>0||a[$4"\t"$1]>0) {print $0}}' \
	GSE190214_${cell_type}.Enhancer_Promoter.chimeric_fragment.out.150nt.EP_anno.uniq.v3.split_test.txt \
	GSE190214_${cell_type}.Enhancer_Promoter.chimeric_fragment.out.150nt.EP_anno.uniq.bedpe >\
	GSE190214_${cell_type}.Enhancer_Promoter.chimeric_fragment.out.150nt.EP_anno.uniq.v3.split_test.bedpe

echo "validate"
grep -f GSE190214_${cell_type}.Enhancer_Promoter.chimeric_fragment.out.150nt.EP_anno.uniq.v3.split_validate.txt \
	GSE190214_${cell_type}.Enhancer_Promoter.chimeric_fragment.out.150nt.EP_anno.uniq.pct.tsv |\
	awk -v tot=$tot '{sum2+=$4}END{print sum2/tot*100"%"}'
awk '{OFS="\t"}FNR==NR{a[$1"\t"$2]+=1;next}{if (a[$1"\t"$4]>0||a[$4"\t"$1]>0) {print $0}}' \
	GSE190214_${cell_type}.Enhancer_Promoter.chimeric_fragment.out.150nt.EP_anno.uniq.v3.split_validate.txt \
	GSE190214_${cell_type}.Enhancer_Promoter.chimeric_fragment.out.150nt.EP_anno.uniq.bedpe >\
	GSE190214_${cell_type}.Enhancer_Promoter.chimeric_fragment.out.150nt.EP_anno.uniq.v3.split_validate.bedpe


#exit

cell_type=H1
awk '{OFS="\t"}{gsub(/chr/,"");if ($1<$4) {print $1,$4} else {print $4,$1}}' \
	GSE190214_${cell_type}.Enhancer_Promoter.chimeric_fragment.out.150nt.EP_anno.uniq.bedpe |\
	sort  | uniq -c | sort -k1,1gr | awk '{OFS="\t"}{print "chr"$2,"chr"$3,$1/95919*100,$1}' >\
	GSE190214_${cell_type}.Enhancer_Promoter.chimeric_fragment.out.150nt.EP_anno.uniq.pct.tsv

python3 ./split_dat.py \
	--input_file GSE190214_${cell_type}.Enhancer_Promoter.chimeric_fragment.out.150nt.EP_anno.uniq.pct.tsv \
	--output_prefix GSE190214_${cell_type}.Enhancer_Promoter.chimeric_fragment.out.150nt.EP_anno.uniq.v3.split \
	--seed 549
echo ${cell_type}
tot=$(grep -f <(cat GSE190214_${cell_type}.Enhancer_Promoter.chimeric_fragment.out.150nt.EP_anno.uniq.v3.split_train.txt \
	GSE190214_${cell_type}.Enhancer_Promoter.chimeric_fragment.out.150nt.EP_anno.uniq.v3.split_test.txt \
	GSE190214_${cell_type}.Enhancer_Promoter.chimeric_fragment.out.150nt.EP_anno.uniq.v3.split_validate.txt) GSE190214_${cell_type}.Enhancer_Promoter.chimeric_fragment.out.150nt.EP_anno.uniq.pct.tsv | awk '{sum+=$4}END{print sum}')
echo "train"
grep -f GSE190214_${cell_type}.Enhancer_Promoter.chimeric_fragment.out.150nt.EP_anno.uniq.v3.split_train.txt \
	GSE190214_${cell_type}.Enhancer_Promoter.chimeric_fragment.out.150nt.EP_anno.uniq.pct.tsv |\
	awk -v tot=$tot '{sum2+=$4}END{print sum2/tot*100"%"}'
awk '{OFS="\t"}FNR==NR{a[$1"\t"$2]+=1;next}{if (a[$1"\t"$4]>0||a[$4"\t"$1]>0) {print $0}}' \
	GSE190214_${cell_type}.Enhancer_Promoter.chimeric_fragment.out.150nt.EP_anno.uniq.v3.split_train.txt \
	GSE190214_${cell_type}.Enhancer_Promoter.chimeric_fragment.out.150nt.EP_anno.uniq.bedpe >\
	GSE190214_${cell_type}.Enhancer_Promoter.chimeric_fragment.out.150nt.EP_anno.uniq.v3.split_train.bedpe
echo "test"
grep -f GSE190214_${cell_type}.Enhancer_Promoter.chimeric_fragment.out.150nt.EP_anno.uniq.v3.split_test.txt \
	GSE190214_${cell_type}.Enhancer_Promoter.chimeric_fragment.out.150nt.EP_anno.uniq.pct.tsv |\
	awk -v tot=$tot '{sum2+=$4}END{print sum2/tot*100"%"}'
awk '{OFS="\t"}FNR==NR{a[$1"\t"$2]+=1;next}{if (a[$1"\t"$4]>0||a[$4"\t"$1]>0) {print $0}}' \
	GSE190214_${cell_type}.Enhancer_Promoter.chimeric_fragment.out.150nt.EP_anno.uniq.v3.split_test.txt \
	GSE190214_${cell_type}.Enhancer_Promoter.chimeric_fragment.out.150nt.EP_anno.uniq.bedpe >\
	GSE190214_${cell_type}.Enhancer_Promoter.chimeric_fragment.out.150nt.EP_anno.uniq.v3.split_test.bedpe

echo "validate"
grep -f GSE190214_${cell_type}.Enhancer_Promoter.chimeric_fragment.out.150nt.EP_anno.uniq.v3.split_validate.txt \
	GSE190214_${cell_type}.Enhancer_Promoter.chimeric_fragment.out.150nt.EP_anno.uniq.pct.tsv |\
	awk -v tot=$tot '{sum2+=$4}END{print sum2/tot*100"%"}'
awk '{OFS="\t"}FNR==NR{a[$1"\t"$2]+=1;next}{if (a[$1"\t"$4]>0||a[$4"\t"$1]>0) {print $0}}' \
	GSE190214_${cell_type}.Enhancer_Promoter.chimeric_fragment.out.150nt.EP_anno.uniq.v3.split_validate.txt \
	GSE190214_${cell_type}.Enhancer_Promoter.chimeric_fragment.out.150nt.EP_anno.uniq.bedpe >\
	GSE190214_${cell_type}.Enhancer_Promoter.chimeric_fragment.out.150nt.EP_anno.uniq.v3.split_validate.bedpe
#exit

cell_type=K562
awk '{OFS="\t"}{gsub(/chr/,"");if ($1<$4) {print $1,$4} else {print $4,$1}}' \
	GSE190214_${cell_type}.Enhancer_Promoter.chimeric_fragment.out.150nt.EP_anno.uniq.bedpe |\
	sort  | uniq -c | sort -k1,1gr | awk '{OFS="\t"}{print "chr"$2,"chr"$3,$1/95919*100,$1}' >\
	GSE190214_${cell_type}.Enhancer_Promoter.chimeric_fragment.out.150nt.EP_anno.uniq.pct.tsv

python3 ./split_dat.py \
	--input_file GSE190214_${cell_type}.Enhancer_Promoter.chimeric_fragment.out.150nt.EP_anno.uniq.pct.tsv \
	--output_prefix GSE190214_${cell_type}.Enhancer_Promoter.chimeric_fragment.out.150nt.EP_anno.uniq.v3.split \
	--seed 392
echo ${cell_type}
tot=$(grep -f <(cat GSE190214_${cell_type}.Enhancer_Promoter.chimeric_fragment.out.150nt.EP_anno.uniq.v3.split_train.txt \
	GSE190214_${cell_type}.Enhancer_Promoter.chimeric_fragment.out.150nt.EP_anno.uniq.v3.split_test.txt \
	GSE190214_${cell_type}.Enhancer_Promoter.chimeric_fragment.out.150nt.EP_anno.uniq.v3.split_validate.txt) GSE190214_${cell_type}.Enhancer_Promoter.chimeric_fragment.out.150nt.EP_anno.uniq.pct.tsv | awk '{sum+=$4}END{print sum}')
echo "train"
grep -f GSE190214_${cell_type}.Enhancer_Promoter.chimeric_fragment.out.150nt.EP_anno.uniq.v3.split_train.txt \
	GSE190214_${cell_type}.Enhancer_Promoter.chimeric_fragment.out.150nt.EP_anno.uniq.pct.tsv |\
	awk -v tot=$tot '{sum2+=$4}END{print sum2/tot*100"%"}'
awk '{OFS="\t"}FNR==NR{a[$1"\t"$2]+=1;next}{if (a[$1"\t"$4]>0||a[$4"\t"$1]>0) {print $0}}' \
	GSE190214_${cell_type}.Enhancer_Promoter.chimeric_fragment.out.150nt.EP_anno.uniq.v3.split_train.txt \
	GSE190214_${cell_type}.Enhancer_Promoter.chimeric_fragment.out.150nt.EP_anno.uniq.bedpe >\
	GSE190214_${cell_type}.Enhancer_Promoter.chimeric_fragment.out.150nt.EP_anno.uniq.v3.split_train.bedpe
echo "test"
grep -f GSE190214_${cell_type}.Enhancer_Promoter.chimeric_fragment.out.150nt.EP_anno.uniq.v3.split_test.txt \
	GSE190214_${cell_type}.Enhancer_Promoter.chimeric_fragment.out.150nt.EP_anno.uniq.pct.tsv |\
	awk -v tot=$tot '{sum2+=$4}END{print sum2/tot*100"%"}'
awk '{OFS="\t"}FNR==NR{a[$1"\t"$2]+=1;next}{if (a[$1"\t"$4]>0||a[$4"\t"$1]>0) {print $0}}' \
	GSE190214_${cell_type}.Enhancer_Promoter.chimeric_fragment.out.150nt.EP_anno.uniq.v3.split_test.txt \
	GSE190214_${cell_type}.Enhancer_Promoter.chimeric_fragment.out.150nt.EP_anno.uniq.bedpe >\
	GSE190214_${cell_type}.Enhancer_Promoter.chimeric_fragment.out.150nt.EP_anno.uniq.v3.split_test.bedpe

echo "validate"
grep -f GSE190214_${cell_type}.Enhancer_Promoter.chimeric_fragment.out.150nt.EP_anno.uniq.v3.split_validate.txt \
	GSE190214_${cell_type}.Enhancer_Promoter.chimeric_fragment.out.150nt.EP_anno.uniq.pct.tsv |\
	awk -v tot=$tot '{sum2+=$4}END{print sum2/tot*100"%"}'
awk '{OFS="\t"}FNR==NR{a[$1"\t"$2]+=1;next}{if (a[$1"\t"$4]>0||a[$4"\t"$1]>0) {print $0}}' \
	GSE190214_${cell_type}.Enhancer_Promoter.chimeric_fragment.out.150nt.EP_anno.uniq.v3.split_validate.txt \
	GSE190214_${cell_type}.Enhancer_Promoter.chimeric_fragment.out.150nt.EP_anno.uniq.bedpe >\
	GSE190214_${cell_type}.Enhancer_Promoter.chimeric_fragment.out.150nt.EP_anno.uniq.v3.split_validate.bedpe
exit




