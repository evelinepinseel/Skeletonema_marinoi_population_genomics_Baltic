# Skeletonema marinoi population genomics Baltic Sea

This repository contains scripts used to analyse SNP data of the diatom *Skeletonema marinoi*, obtained via Pool-seq.

1. **filter-pool-seq-sync.R** filters a sync file (see PoPoolation2 for more information on this file format) using specified parameters for minimum and maximum coverage, minimum allele count, and minimum allele frequency. The script also removes multiallelic SNPs (i.e., >2 alleles) and calculates allele frequencies.

2. **sync2VCF.R**: converts a sync file into a VCF file with a format that is compatible for SNP annotation with snpEff.

Both scripts can be run using the terminal.

**Citation:** please cite: Pinseel E., Ruck E.C., Nakov T., Jonsson P., Kourtchenko O., Kremp A., Pinder M.I.M., Roberts W.R., Sjöqvist C., Töpel M., Godhe A, Hahn M.H., Alverson A.J.  (2023). ​Local adaptation of a marine diatom is governed by genome-wide changes in diverse metabolic processes. BioRxiv. doi: 10.1101/2023.09.22.559080.
