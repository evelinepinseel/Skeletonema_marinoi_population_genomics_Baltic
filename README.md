# Skeletonema marinoi population genomics Baltic Sea

This repository contains scripts used to analyse SNP data of the diatom *Skeletonema marinoi*, obtained via Pool-seq.

1. **filter-pool-seq-sync.R** filters a sync file (see PoPoolation2 for more information on this file format) using specified parameters for minimum and maximum coverage, minimum allele count, and minimum allele frequency. The script also removes multiallelic SNPs (i.e., >2 alleles) and calculates allele frequencies.

2. **sync2VCF.R**: converts a sync file into a VCF file with a format that is compatible for SNP annotation with snpEff.

Both scripts can be run using the terminal.

**Citation:** please cite this GitHub repository.
