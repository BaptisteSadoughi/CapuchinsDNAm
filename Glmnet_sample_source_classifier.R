####### DNAm in capuchin research project
## This script uses a glmnet multinomial classifier to test the ability to
## differentiate between fecal, urine and blood samples based on methylation profiles at a reduced set of CpGs.

rm(list = ls())

library_list <- c("tidyverse","stringr","comethyl","GenomicRanges","bsseq","missMDA","glmnet")
lapply(library_list, require, character.only=TRUE)

##### Data Prep
#####
{
  # Load Twist capture
  data_ttms <- readRDS("/path/to/bsseq/bismarkBSseq_capuchin_feces.rds")
  
  # Simplify sample names
  column_names <- colnames(data_ttms)
  # change first 121 columns
  column_names[1:121] <- gsub(pattern =".CpG_report.merged_CpG_evidence.cov.gz", "",
                              str_split_i(column_names[1:121],"/",6))
  
  # change remaining columns
  column_names[122:length(column_names)] <- gsub(pattern =".CpG_report.merged_CpG_evidence.cov.gz", "",
                                                 str_split_i(column_names[122:length(column_names)],"/",7))
  
  # Assign the new column names to the dataframe
  colnames(data_ttms) <- column_names
  
  # Load metadata #Table S1
  metadata=read.csv("/path/to/metadata/metadata_capuchins.csv")
  
  #### DROP samples with low mapping efficiency or conversion CHH >2 or duplicated or never mapped or found to have bad coverage post filtering
  samples_to_exclude <- c("11432-NS-0095", "11432-NS-0116", "11432-NS-0069", "11432-NS-0091", "11432-NS-0083", "11432-NS-0088",
                          "11432-NS-0024", "11432-NS-0105", "11432-NS-0118", "11432-NS-0096", "11432-NS-0102", "11432-NS-0077",
                          "11432-NS-0100", "11432-NS-0099", "11432-NS-0072", "11432-NS-0092", "11432-NS-0106", "11432-NS-0064",
                          "11432-NS-0070", "11432-NS-0115", "11432-NS-0075", "11119-AJ-0172", "11119-AJ-0171", "11119-AJ-0168",
                          "11119-AJ-0165", "11119-AJ-0173", "11119-AJ-0167", "11119-AJ-0164", "11119-AJ-0170", "11119-AJ-0174",
                          "11119-AJ-0175")
  
  metadata_clean <- metadata %>% filter(!sample_name %in% samples_to_exclude)
  
  # For filtering we separate fecal and blood samples
  
  # Blood
  data_ttms_blood = data_ttms[,sampleNames(data_ttms) %in%  metadata_clean$sample_name[metadata_clean$sample_type=="blood"]]
  
  metadata_blood <- metadata_clean[metadata_clean$sample_name %in% colnames(data_ttms_blood),]
  metadata_blood <- metadata_blood[match(colnames(data_ttms_blood),metadata_blood$sample_name),]
    
  # Fecal
  data_ttms_fecal = data_ttms[,sampleNames(data_ttms) %in%  metadata_clean$sample_name[metadata_clean$sample_type=="fecal"]]
  
  metadata_fecal <- metadata_clean[metadata_clean$sample_name %in% colnames(data_ttms_fecal),]
  metadata_fecal <- metadata_fecal[match(colnames(data_ttms_fecal),metadata_fecal$sample_name),]
  
  # Urine
  data_ttms_urine = data_ttms[,sampleNames(data_ttms) %in%  metadata_clean$sample_name[metadata_clean$sample_type=="urine"]]
  
  metadata_urine <- metadata_clean[metadata_clean$sample_name %in% colnames(data_ttms_urine),]
  metadata_urine <- metadata_urine[match(colnames(data_ttms_urine),metadata_urine$sample_name),]
}
identical(colnames(data_ttms_urine), metadata_urine$sample_name)
identical(colnames(data_ttms_fecal), metadata_fecal$sample_name)
identical(colnames(data_ttms_blood), metadata_blood$sample_name)


##########################################################################
#                                                                        #
#### Section 2: Multidimensional analysis of methylation profiles        #
#                                                                        #
##########################################################################

# Filter datasets independently
data_ttms_blood_cov5p75 <- filterCpGs(data_ttms_blood, cov = 5, perSample = 0.75,
                                      verbose = FALSE, save = FALSE, file = NULL)


data_ttms_fecal_cebus = data_ttms_fecal[,sampleNames(data_ttms_fecal) %in%  metadata_fecal$sample_name[metadata_fecal$species=="Cebus"]]
data_ttms_fecal_cebus_cov5p75 <- filterCpGs(data_ttms_fecal_cebus, cov = 5, perSample = 0.75,
                                            verbose = FALSE, save = FALSE, file = NULL)


data_ttms_fecal_sapajus = data_ttms_fecal[,sampleNames(data_ttms_fecal) %in%  metadata_fecal$sample_name[metadata_fecal$species=="Sapajus"]]
data_ttms_fecal_sapajus_cov5p75 <- filterCpGs(data_ttms_fecal_sapajus, cov = 5, perSample = 0.75,
                                              verbose = FALSE, save = FALSE, file = NULL)


data_ttms_urine_cov5p75 <- filterCpGs(data_ttms_urine, cov = 5, perSample = 0.75,
                                      verbose = FALSE, save = FALSE, file = NULL)

# Intersect datasets
source("Functions_intersect_bsseq.R") #it may be necessary to adapt the function to pass the desired number of bsseq objects

# Apply function to the bsseq objects
result_intersect_fecal_species <- extract_common_cpgs_four(
  bsseq_1 = data_ttms_fecal_cebus_cov5p75, 
  bsseq_2 = data_ttms_fecal_sapajus_cov5p75, 
  bsseq_3 = data_ttms_blood_cov5p75,
  bsseq_4 = data_ttms_urine_cov5p75,
  bsseq_common_1_name = "fecal_cebus_common", 
  bsseq_common_2_name = "fecal_sapajus_common",
  bsseq_common_3_name = "blood_sapajus_common",
  bsseq_common_4_name = "urine_sapajus_common"
)


data_ttms_fecal_cebus_cov5p75_common <- result_intersect_fecal_species$fecal_cebus_common
data_ttms_fecal_sapajus_cov5p75_common <- result_intersect_fecal_species$fecal_sapajus_common
data_ttms_blood_sapajus_cov5p75_common <- result_intersect_fecal_species$blood_sapajus_common
data_ttms_urine_sapajus_cov5p75_common <- result_intersect_fecal_species$urine_sapajus_common

# add percent methylation matrices
data_ttms_fecal_cebus_cov5p75_common@assays@data@listData[["pmeth"]]=getCoverage(data_ttms_fecal_cebus_cov5p75_common, type = "M")/getCoverage(data_ttms_fecal_cebus_cov5p75_common, type = "Cov")
data_ttms_fecal_sapajus_cov5p75_common@assays@data@listData[["pmeth"]]=getCoverage(data_ttms_fecal_sapajus_cov5p75_common, type = "M")/getCoverage(data_ttms_fecal_sapajus_cov5p75_common, type = "Cov")
data_ttms_blood_sapajus_cov5p75_common@assays@data@listData[["pmeth"]]=getCoverage(data_ttms_blood_sapajus_cov5p75_common, type = "M")/getCoverage(data_ttms_blood_sapajus_cov5p75_common, type = "Cov")
data_ttms_urine_sapajus_cov5p75_common@assays@data@listData[["pmeth"]]=getCoverage(data_ttms_urine_sapajus_cov5p75_common, type = "M")/getCoverage(data_ttms_urine_sapajus_cov5p75_common, type = "Cov")


pmeth_fecal_cebus_cov5p75_common <- data_ttms_fecal_cebus_cov5p75_common@assays@data@listData$pmeth
pmeth_fecal_cebus_cov5p75_common <- as.matrix(pmeth_fecal_cebus_cov5p75_common)

pmeth_fecal_sapajus_cov5p75_common <- data_ttms_fecal_sapajus_cov5p75_common@assays@data@listData$pmeth
pmeth_fecal_sapajus_cov5p75_common <- as.matrix(pmeth_fecal_sapajus_cov5p75_common)

pmeth_blood_sapajus_cov5p75_common <- data_ttms_blood_sapajus_cov5p75_common@assays@data@listData$pmeth
pmeth_blood_sapajus_cov5p75_common <- as.matrix(pmeth_blood_sapajus_cov5p75_common)

pmeth_urine_sapajus_cov5p75_common <- data_ttms_urine_sapajus_cov5p75_common@assays@data@listData$pmeth
pmeth_urine_sapajus_cov5p75_common <- as.matrix(pmeth_urine_sapajus_cov5p75_common)


# Use the cpg coordinate for the methylation matrices
cpg_loc_all_cov5p75_common <- granges(data_ttms_fecal_cebus_cov5p75_common)
cpg_loc_all_cov5p75_common <- paste(seqnames(cpg_loc_all_cov5p75_common), start(cpg_loc_all_cov5p75_common), sep = "_")

rownames(pmeth_fecal_sapajus_cov5p75_common)<- cpg_loc_all_cov5p75_common
rownames(pmeth_fecal_cebus_cov5p75_common)<- cpg_loc_all_cov5p75_common
rownames(pmeth_blood_sapajus_cov5p75_common)<- cpg_loc_all_cov5p75_common
rownames(pmeth_urine_sapajus_cov5p75_common)<- cpg_loc_all_cov5p75_common

######### Impute missing methylation in the three datasets independently

# First remove site with low variance as they are not processed by missMDA

site_var_fecalsapa <- apply(pmeth_fecal_sapajus_cov5p75_common, 1, function(x) sd(x,na.rm = TRUE))
pmeth_fecal_sapajus_cov5p75_common_var <- pmeth_fecal_sapajus_cov5p75_common[site_var_fecalsapa > 0.05,]

site_var_fecalceb <- apply(pmeth_fecal_cebus_cov5p75_common, 1, function(x) sd(x,na.rm = TRUE))
pmeth_fecal_cebus_cov5p75_common_var <-  pmeth_fecal_cebus_cov5p75_common[site_var_fecalceb > 0.05,]

site_var_bloodsapa <- apply(pmeth_blood_sapajus_cov5p75_common, 1, function(x) sd(x,na.rm = TRUE))
pmeth_blood_sapajus_cov5p75_common_var <-  pmeth_blood_sapajus_cov5p75_common[site_var_bloodsapa > 0.05,]

site_var_urinesapa <- apply(pmeth_urine_sapajus_cov5p75_common, 1, function(x) sd(x,na.rm = TRUE))
pmeth_urine_sapajus_cov5p75_common_var <-  pmeth_urine_sapajus_cov5p75_common[site_var_urinesapa > 0.05,]

# ncp_best = estim_ncpPCA(t(pmeth_fecal_both_cov5p75_common_var[1:5000,]), scale = TRUE, ncp.min=0, ncp.max=5, method.cv = "Kfold")$ncp #ncp_best=2
pmeth_fecal_sapajus_cov5p75_common_var_imp=t(imputePCA(t(pmeth_fecal_sapajus_cov5p75_common_var), scale = TRUE, ncp = 2)$completeObs)
pmeth_fecal_cebus_cov5p75_common_var_imp=t(imputePCA(t(pmeth_fecal_cebus_cov5p75_common_var), scale = TRUE, ncp = 2)$completeObs)
pmeth_blood_sapajus_cov5p75_common_var_imp=t(imputePCA(t(pmeth_blood_sapajus_cov5p75_common_var), scale = TRUE, ncp = 2)$completeObs)
pmeth_urine_sapajus_cov5p75_common_var_imp=t(imputePCA(t(pmeth_urine_sapajus_cov5p75_common_var), scale = TRUE, ncp = 2)$completeObs)

## Given that we lost different sites in each matrix, we filter once again for common sites

final_site_set <- intersect(rownames(pmeth_fecal_sapajus_cov5p75_common_var_imp),
                            intersect(rownames(pmeth_fecal_cebus_cov5p75_common_var_imp),
                                      intersect(rownames(pmeth_blood_sapajus_cov5p75_common_var_imp),rownames(pmeth_urine_sapajus_cov5p75_common_var_imp))
                                      )
                            )

pmeth_blood_sapajus_cov5p75_common_var_imp_f <- pmeth_blood_sapajus_cov5p75_common_var_imp[rownames(pmeth_blood_sapajus_cov5p75_common_var_imp) %in% final_site_set,]
pmeth_fecal_sapajus_cov5p75_common_var_imp_f <- pmeth_fecal_sapajus_cov5p75_common_var_imp[rownames(pmeth_fecal_sapajus_cov5p75_common_var_imp) %in% final_site_set,]
pmeth_fecal_cebus_cov5p75_common_var_imp_f <- pmeth_fecal_cebus_cov5p75_common_var_imp[rownames(pmeth_fecal_cebus_cov5p75_common_var_imp) %in% final_site_set,]
pmeth_urine_sapajus_cov5p75_common_var_imp_f <- pmeth_urine_sapajus_cov5p75_common_var_imp[rownames(pmeth_urine_sapajus_cov5p75_common_var_imp) %in% final_site_set,]

identical(rownames(pmeth_fecal_cebus_cov5p75_common_var_imp_f), rownames(pmeth_fecal_sapajus_cov5p75_common_var_imp_f))
identical(rownames(pmeth_fecal_cebus_cov5p75_common_var_imp_f), rownames(pmeth_urine_sapajus_cov5p75_common_var_imp_f))
identical(rownames(pmeth_fecal_cebus_cov5p75_common_var_imp_f), rownames(pmeth_blood_sapajus_cov5p75_common_var_imp_f))

##########################################################################
#                                                                        #
#### Section 3: MULTINOMIAL GLMNET                                       #
#                                                                        #
##########################################################################

# Assemble back the three matrices
pmeth_all_cov5p75_common <- cbind(pmeth_fecal_cebus_cov5p75_common_var_imp_f,cbind(pmeth_fecal_sapajus_cov5p75_common_var_imp_f,cbind(pmeth_urine_sapajus_cov5p75_common_var_imp_f,pmeth_blood_sapajus_cov5p75_common_var_imp_f)))

# Simplify metadata for glmnet
meta <- metadata_clean %>% dplyr::select(sample_name,monkey_name,sample_type)

meta <- meta[match(colnames(pmeth_all_cov5p75_common), meta$sample_name),]

identical(meta$sample_name, colnames(pmeth_all_cov5p75_common))

# check that no residual NA remains in the data
any(is.na(pmeth_all_cov5p75_common))


#### FUNCTION
run_glmnet_multinomial <- function(SAMP, pmeth_matrix, meta) {
  
  # Ensure SAMP is an integer index
  SAMP <- as.integer(SAMP)
  
  # Extract training and test sets
  train <- pmeth_matrix[, -SAMP]
  test <- pmeth_matrix[, SAMP]
  
  # Extract corresponding metadata
  train_sampletype <- meta$sample_type[-SAMP]
  test_sampletype <- meta$sample_type[SAMP]
  test_monkeyname <- meta$monkey_name[SAMP]
  test_samplename <- meta$sample_name[SAMP]
  
  # Fit the glmnet model
  mod.glmnet <- cv.glmnet(x = t(train), y = train_sampletype, family = "multinomial", type.multinomial = "grouped")
  
  # Predict the response for the test data
  predicted_matrix <- predict(mod.glmnet, newx = t(test), s = "lambda.min", type = "response")
  
  # Convert the matrix to a data frame and name the columns
  pred <- as.data.frame(predicted_matrix)
  colnames(pred) <- c("blood", "fecal", "urine")
  
  # Find the best guess (the column name with the maximum value)
  best_guess <- colnames(pred)[which.max(pred)]
  
  # Create a result data frame
  result_samp <- data.frame(
    sample_name = test_samplename,
    sample_type = test_sampletype,
    monkey_name = test_monkeyname,
    blood_proba = as.numeric(pred[, "blood"]),
    fecal_proba = as.numeric(pred[, "fecal"]),
    urine_proba = as.numeric(pred[, "urine"]),
    mislabel = ifelse(test_sampletype != best_guess, 1, 0),
    stringsAsFactors = FALSE
  )
  
  return(result_samp)
}


# Apply to the methylation matrix
glmnet_LOOC_list <- mclapply(1:nrow(meta),
                             run_glmnet_multinomial,
                             pmeth_matrix = pmeth_all_cov5p75_common,
                             meta = meta,
                             mc.cores = 8)

glmnet_LOOC_df <- do.call(rbind,glmnet_LOOC_list)

############# END OF THE CODE #################
