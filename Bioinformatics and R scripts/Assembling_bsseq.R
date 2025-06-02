#!/usr/bin/env /packages/apps/spack/18/opt/spack/gcc-11.2.0/r-4.2.2-kpl/bin/Rscript

# sbatch --cpus-per-task=1 --mem=180G -p general -q public -t 0-04:00:00 /path/to/script/Assembling_bsseq.R

library(Biostrings)
library(bsseq)
library(rlang)

#The samples sequences are spread across several locations. Careful to retrieve all available data.
files_to_load <- list.files(path = "/path/to/final_cov", full.names = TRUE, pattern = ".cov.gz")

bismarkBSseq <- read.bismark(files = files_to_load,
                             strandCollapse = TRUE,
                             verbose = TRUE,
                             BACKEND = "HDF5Array",
                             replace=TRUE,
                             rmZeroCov = FALSE,
                             dir="/path/to/bsseq/hdf5")

saveRDS(bismarkBSseq,file="/path/to/bsseq/bismarkBSseq_capuchin_feces.rds")

q(save="no")
