# CapuchinsDNAm

**Repository for the project on the non-invasive measurement of DNA methylation in capuchin monkeys.**


This repository contains the scripts used to analyze DNA methylation profiles measured from blood, fecal, and urine samples in two species of capuchin monkeys.

Note that all steps were performed on the Arizona State University SOL higher-performance computing clusters. We have used R version 4.4.0 in the RStudio interface.
The first processing steps of raw fastq files were performed using a bash shell, with required software and binaries in the PATH. Versions for softwares used in the shell
environment are specified in the corresponding scripts. We did our best to keep R packages up-to-date until the end of the analysis in May 2025.

Below we detail the steps to reproduce the results presented in the associated manuscript. When relevant, we highlight the need to define the paths to external files. Note that some scripts require to source the script Functions_intersect_bsseq.R in order to use the functions.

**1. Bioinformatic_processing.sh**
Download raw fastq files and indicate paths the data /path/to/fastq. A manifest file listing all library names should be stored in the same folder /path/to/fastq/manifest.txt. Update path.
Download the C.imitator whole genome Cebus_imitator.Cebus_imitator-1.0.dna.toplevel.fa. Update path.

**2. Assembling_bsseq.R**
Update output path for the bsseq formatted data.

**3. TissueMarkers_fecal.R**
Download data from HumanMethylationAtlas Loyfer et al.2024 Supplementary Table S4C. List of 50286 cell type-specific unmethylated markers (top 1000, hg19). Update path.
Update path to the bsseq formatted data (valid for all scripts loading the bsseq data).
Download capuchin metadata in Table S1. Update path (valid for all scripts loading the metadata).
Download Cebus_imitator.Cebus_imitator-1.0.113.gtf. Update path.

**4. Multidimensional_analysis_DNAm.R**

**5. Glmnet_sample_source_classifier.R**

**6. DNAm_clock_mixedsamplesource.R**
Update output/input paths to export/load formatted percent methylation matrices.
Update output/input paths to export/load results of elastic net regressions.

**7. DNAm_clock_Cebusfecal.R**
Update output/input paths to export/load formatted percent methylation matrices.
Update output/input paths to export/load results of elastic net regressions.

**8. Age_DMSs.R**
Update output/input paths to export/load the results of the differential methylation analysis.
Update path to the Cebus_imitator.Cebus_imitator-1.0.113.gtf.
Download the panmammalian data Lu et al. 2023 Supplementary Table S6.1_TopMetaEWASinAllTissues. Update path.

**9. Capuchin_Mammalian_EWAS.R**
Update path to the panmammalian data Lu et al. 2023 Supplementary Table S6.1_TopMetaEWASinAllTissues.
Update path to the Cebus_imitator.Cebus_imitator-1.0.113.gtf.
Update paths to the results of the differential methylation analysis.

**10. Make_DMRs.R**
Update paths to the results of the differential methylation analysis.
Update output paths for the age-associated DMRs and the set of background DMRs.

**11. TFBS_analysis.R**
Update put paths to the age-associated DMRs and the set of background DMRs in the script.
Update path to the C.imitator whole genome Cebus_imitator.Cebus_imitator-1.0.dna.toplevel.fa.
