# CapuchinsDNAm

**Repository for the project on the non-invasive measurement of DNA methylation in capuchin monkeys.**


This repository contains the scripts used to analyze DNA methylation profiles measured from blood, fecal, and urine samples in two species of capuchin monkeys.

Note that all steps were performed on the Arizona State University SOL higher-performance computing clusters. We have used R version 4.4.0 in the RStudio interface.
The first processing steps of raw fastq files were performed using a bash shell, with required software and binaries in the PATH. Versions for softwares used in the shell
environment are specified in the corresponding scripts. We did our best to keep R packages up-to-date until the end of the analysis in May 2025.

Below we detail the steps to reproduce the results presented in the associated manuscript. When relevant, we highlight the need to define the paths to external files.

1. Bioinformatic_processing.sh
Download raw fastq files and indicate paths the data /path/to/fastq. A manifest file listing all library names should be stored in the same folder /path/to/fastq/manifest.txt. Update path in the script.
Download the C.imitator whole genome Cebus_imitator.Cebus_imitator-1.0.dna.toplevel.fa. Update path in the script.

3. Assembling_bsseq.R

4. TissueMarkers_fecal.R
Download data from HumanMethylationAtlas Loyfer et al.2024 Supplementary Table S4C. List of 50286 cell type-specific unmethylated markers (top 1000, hg19). Update path in the script.
Download capuchin metadata in Table S1. Update path in the script (valid for all scripts loading the metadata).
Download Cebus_imitator.Cebus_imitator-1.0.113.gtf. Update path in the script.

5. Multidimensional_analysis_DNAm.R

6. DNAm_clock_mixedsamplesource.R
Update outpath to export results of elastic net regressions.

7. 
