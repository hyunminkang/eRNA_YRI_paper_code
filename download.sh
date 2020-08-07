#!/bin/sh
cd TRE_identification/data
wget ftp://share.sph.umich.edu/public/GSE110638_example_data/procap.pl.sorted.bg.gz
wget ftp://share.sph.umich.edu/public/GSE110638_example_data/procap.mn.sorted.bg.gz
cd ../..
cd variant_alignment/genotype/
wget ftp://share.sph.umich.edu/public/GSE110638_example_data/genotype.imputed.vcf.gz
cd ../..
