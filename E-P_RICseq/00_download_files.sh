#!/bin/bash
# MIT License
# Copyright (c) 2025 Eric Nels Pederson, University of Massachusetts Amherst
# Licensed under the MIT License. See LICENSE file in the project root for full license information.

# Network files for locations of enhancer and promoter regions
wget https://ftp.ncbi.nlm.nih.gov/geo/series/GSE190nnn/GSE190214/suppl/GSE190214%5FEP%5Fand%5FPP%2EGM12878%2Enetwork%2Exlsx
wget https://ftp.ncbi.nlm.nih.gov/geo/series/GSE190nnn/GSE190214/suppl/GSE190214%5FEP%5Fand%5FPP%2EH1%2Enetwork%2Exlsx
wget https://ftp.ncbi.nlm.nih.gov/geo/series/GSE190nnn/GSE190214/suppl/GSE190214%5FEP%5Fand%5FPP%2EHeLa%2Enetwork%2Exlsx
wget https://ftp.ncbi.nlm.nih.gov/geo/series/GSE190nnn/GSE190214/suppl/GSE190214%5FEP%5Fand%5FPP%2EHepG2%2Enetwork%2Exlsx
wget https://ftp.ncbi.nlm.nih.gov/geo/series/GSE190nnn/GSE190214/suppl/GSE190214%5FEP%5Fand%5FPP%2EIMR90%2Enetwork%2Exlsx
wget https://ftp.ncbi.nlm.nih.gov/geo/series/GSE190nnn/GSE190214/suppl/GSE190214%5FEP%5Fand%5FPP%2EK562%2Enetwork%2Exlsx
wget https://ftp.ncbi.nlm.nih.gov/geo/series/GSE190nnn/GSE190214/suppl/GSE190214%5FEP%5Fand%5FPP%2EhNPC%2Enetwork%2Exlsx

# Bed files with data for kmeRNA
wget https://ftp.ncbi.nlm.nih.gov/geo/series/GSE190nnn/GSE190214/suppl/GSE190214%5FGM12878%2EEnhancer%5FPromoter%2Echimeric%5Ffragment%2Ebed%2Egz
wget https://ftp.ncbi.nlm.nih.gov/geo/series/GSE190nnn/GSE190214/suppl/GSE190214%5FH1%2EEnhancer%5FPromoter%2Echimeric%5Ffragment%2Ebed%2Egz
wget https://ftp.ncbi.nlm.nih.gov/geo/series/GSE190nnn/GSE190214/suppl/GSE190214%5FHeLa%2EEnhancer%5FPromoter%2Echimeric%5Ffragment%2Ebed%2Egz
wget https://ftp.ncbi.nlm.nih.gov/geo/series/GSE190nnn/GSE190214/suppl/GSE190214%5FHepG2%2EEnhancer%5FPromoter%2Echimeric%5Ffragment%2Ebed%2Egz
wget https://ftp.ncbi.nlm.nih.gov/geo/series/GSE190nnn/GSE190214/suppl/GSE190214%5FIMR90%2EEnhancer%5FPromoter%2Echimeric%5Ffragment%2Ebed%2Egz
wget https://ftp.ncbi.nlm.nih.gov/geo/series/GSE190nnn/GSE190214/suppl/GSE190214%5FK562%2EEnhancer%5FPromoter%2Echimeric%5Ffragment%2Ebed%2Egz
wget https://ftp.ncbi.nlm.nih.gov/geo/series/GSE190nnn/GSE190214/suppl/GSE190214%5FhNPC%2EEnhancer%5FPromoter%2Echimeric%5Ffragment%2Ebed%2Egz
# Genome files

cd genome
wget http://ftp.ensembl.org/pub/release-75/fasta/homo_sapiens/dna/Homo_sapiens.GRCh37.75.dna_sm.toplevel.fa.gz
gunzip Homo_sapiens.GRCh37.75.dna_sm.toplevel.fa.gz
samtools faidx Homo_sapiens.GRCh37.75.dna_sm.toplevel.fa

wget http://ftp.ensembl.org/pub/release-75/gtf/homo_sapiens/Homo_sapiens.GRCh37.75.gtf.gz
gunzip Homo_sapiens.GRCh37.75.gtf.gz


