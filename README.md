# Skeletonema marinoi population genomics Baltic Sea

This repository contains two scripts used to analyse SNP data of the diatom Skeletonema marinoi.

1. **filter-pool-seq-sync.R** filters a sync file (see PoPoolation2 for more information on this format) using specified parameters for minimum and maximum coverage, minimum allele count, and minimum allele frequency. The script also calculates allele frequencies.

2. **sync2VCF.R**: converts a sync file into a VCF file with a format that is compatible for SNP annotation with snpEff.
