#!/bin/bash

#mv GSE190214_HeLa.Enhancer_Promoter.chimeric_fragment.out.150nt.EP_anno.uniq.noncoding_seed42813.split_train.DN.enhancers.fa.pairs.RNA.out.tsv.gz \
#	GSE190214_HeLa.Enhancer_Promoter.chimeric_fragment.out.150nt.split_train.DN.enhancer-promoters.tsv.gz
#mv GSE190214_HeLa.Enhancer_Promoter.chimeric_fragment.out.150nt.EP_anno.uniq.noncoding_seed42813.split_train.enhancers.fa.pairs.RNA.out.tsv.gz \
#	GSE190214_HeLa.Enhancer_Promoter.chimeric_fragment.out.150nt.split_train.enhancers.tsv.gz
#mv GSE190214_HeLa.Enhancer_Promoter.chimeric_fragment.out.150nt.EP_anno.uniq.noncoding_seed42813.split_train.promoters.fa.pairs.RNA.out.tsv.gz \
#	GSE190214_HeLa.Enhancer_Promoter.chimeric_fragment.out.150nt.split_train.promoters.tsv.gz

#mv GSE190214_HeLa.Enhancer_Promoter.chimeric_fragment.out.150nt.EP_anno.uniq.noncoding_seed42813.split_test.DN.enhancers.fa.pairs.RNA.out.tsv.gz \
#	GSE190214_HeLa.Enhancer_Promoter.chimeric_fragment.out.150nt.split_test.DN.enhancer-promoters.tsv.gz
#mv GSE190214_HeLa.Enhancer_Promoter.chimeric_fragment.out.150nt.EP_anno.uniq.noncoding_seed42813.split_test.enhancers.fa.pairs.RNA.out.tsv.gz \
#	GSE190214_HeLa.Enhancer_Promoter.chimeric_fragment.out.150nt.split_test.enhancers.tsv.gz
#mv GSE190214_HeLa.Enhancer_Promoter.chimeric_fragment.out.150nt.EP_anno.uniq.noncoding_seed42813.split_test.promoters.fa.pairs.RNA.out.tsv.gz \
#	GSE190214_HeLa.Enhancer_Promoter.chimeric_fragment.out.150nt.split_test.promoters.tsv.gz

for cell_type2 in hNPC GM12878 H1 IMR90 HepG2 K562
do
	mv GSE190214_HeLa-${cell_type2}.Enhancer_Promoter.chimeric_fragment.out.150nt.EP_anno.uniq.noncoding_seed42813.split_test.DN.enhancers.fa.pairs.RNA.out.tsv.gz \
		GSE190214_HeLa-${cell_type2}.Enhancer_Promoter.chimeric_fragment.out.150nt.split_test.DN.enhancer-promoters.tsv.gz
	mv GSE190214_HeLa-${cell_type2}.Enhancer_Promoter.chimeric_fragment.out.150nt.EP_anno.uniq.noncoding_seed42813.split_test.enhancers.fa.pairs.RNA.out.tsv.gz \
		GSE190214_HeLa-${cell_type2}.Enhancer_Promoter.chimeric_fragment.out.150nt.split_test.enhancers.tsv.gz
 	mv GSE190214_HeLa-${cell_type2}.Enhancer_Promoter.chimeric_fragment.out.150nt.EP_anno.uniq.noncoding_seed42813.split_test.promoters.fa.pairs.RNA.out.tsv.gz \
		GSE190214_HeLa-${cell_type2}.Enhancer_Promoter.chimeric_fragment.out.150nt.split_test.promoters.tsv.gz
 

done




