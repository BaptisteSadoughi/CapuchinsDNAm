# CapuchinsDNAm

**Repository for the project on the non-invasive measurement of DNA methylation in capuchin monkeys.**


This repository contains the scripts used to analyze DNA methylation profiles measured from blood, fecal, and urine samples in two species of capuchin monkeys.

Note that all steps were performed on the Arizona State University SOL higher-performance computing cluster. We have used R version 4.4.0 in the RStudio interface.
The first processing steps of raw fastq files were performed using a bash shell, with required software and binaries in the PATH. Versions for softwares used in the shell
environment are specified in the corresponding scripts. We did our best to keep R packages up-to-date until the end of the analysis in May 2025.

Below we detail the steps to reproduce the results presented in the associated manuscript. When relevant, we highlight the need to define the paths to external files. Note that some scripts require to source the script Functions_intersect_bsseq.R in order to use the functions.

**1. Bioinformatic_processing.sh**<br />
Download raw fastq files and indicate the path to the data /path/to/fastq. A manifest file listing all library names should be stored in the same folder /path/to/fastq/manifest.txt. Update path.<br />
Download the C.imitator whole genome Cebus_imitator.Cebus_imitator-1.0.dna.toplevel.fa. Update path.<br />

**2. Assembling_bsseq.R**<br />
Update output path for the bsseq formatted data.<br />

**3. TissueMarkers_fecal.R**<br />
Download data from HumanMethylationAtlas Loyfer et al.2024 Supplementary Table S4C. List of 50286 cell type-specific unmethylated markers (top 1000, hg19). Update path.<br />
Update path to the bsseq formatted data (valid for all scripts loading the bsseq data).<br />
Download capuchin metadata in Table S1. Update path (valid for all scripts loading the metadata).<br />
Download Cebus_imitator.Cebus_imitator-1.0.113.gtf. Update path.<br />

**4. Multidimensional_analysis_DNAm.R**<br />

**5. Glmnet_sample_source_classifier.R**<br />

**6. DNAm_clock_mixedsamplesource.R**<br />
Update output/input paths to export/load formatted percent methylation matrices.<br />
Update output/input paths to export/load results of elastic net regressions.<br />

**7. DNAm_clock_Cebusfecal.R**<br />
Update output/input paths to export/load formatted percent methylation matrices.<br />
Update output/input paths to export/load results of elastic net regressions.<br />

**8. Age_DMSs.R**<br />
Update output/input paths to export/load the results of the differential methylation analysis.<br />
Update path to the Cebus_imitator.Cebus_imitator-1.0.113.gtf.<br />
Download the panmammalian data Lu et al. 2023 Supplementary Table S6.1_TopMetaEWASinAllTissues. Update path.<br />

**9. Capuchin_Mammalian_EWAS.R**<br />
Update path to the panmammalian data Lu et al. 2023 Supplementary Table S6.1_TopMetaEWASinAllTissues.<br />
Update path to the Cebus_imitator.Cebus_imitator-1.0.113.gtf.<br />
Update paths to the results of the differential methylation analysis.<br />

**10. Make_DMRs.R**<br />
Update paths to the results of the differential methylation analysis.<br />
Update output paths for the age-associated DMRs and the set of background DMRs.<br />

**11. TFBS_analysis.R**<br />
Update put paths to the age-associated DMRs and the set of background DMRs in the script.<br />
Update path to the C.imitator whole genome Cebus_imitator.Cebus_imitator-1.0.dna.toplevel.fa.<br />
