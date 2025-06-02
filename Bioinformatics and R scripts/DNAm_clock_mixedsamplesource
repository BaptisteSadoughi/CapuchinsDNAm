####### DNAm in capuchin research project
## Build DNAm clock from blood, feces, and urine samples.

### DISCLAIMER: the glmnet algorithm is this script was run using 16 cores in parallel with 60GB on the SOL supercomputer. If needed, the function call can easily be extracted and submitted as a job.

########### Prep data for clock analysis

rm(list = ls())

library_list <- c("tidyverse","stringr","comethyl","GenomicRanges","bsseq", "glmnet","missMDA","bestNormalize")
lapply(library_list, require, character.only=TRUE)

# Function to compute clock performance
corr_and_mae <- function(df){
  coef <- round(cor.test(df$age,df$predicted_age, method = "pearson")$estimate,2)
  MAE <- round(median(abs(df$predicted_age - df$age)), digits = 2)
  return(c(coef, MAE= MAE))}

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
#### Section 2: Filtering and combination of all sample types            #
#                                                                        #
##########################################################################

# Filter datasets independently
data_ttms_blood_cov5p75 <- filterCpGs(data_ttms_blood, cov = 5, perSample = 0.75,
                                      verbose = FALSE, save = FALSE, file = NULL)


data_ttms_fecal_cebus = data_ttms_fecal[,sampleNames(data_ttms_fecal) %in%  metadata_fecal$sample_name[metadata_fecal$species=="Cebus" & metadata_fecal$estimated_dob == "n"]]
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

############# Impute missing methylation in the three datasets independently

# First remove site with low variance as they are not processed by missMDA

site_var_fecalsapa <- apply(pmeth_fecal_sapajus_cov5p75_common, 1, function(x) mean(x,na.rm = TRUE))
pmeth_fecal_sapajus_cov5p75_common_var <- pmeth_fecal_sapajus_cov5p75_common[site_var_fecalsapa > 0.1 & site_var_fecalsapa <0.9,]

site_var_fecalceb <- apply(pmeth_fecal_cebus_cov5p75_common, 1, function(x) mean(x,na.rm = TRUE))
pmeth_fecal_cebus_cov5p75_common_var <-  pmeth_fecal_cebus_cov5p75_common[site_var_fecalceb > 0.1 & site_var_fecalceb <0.9,]

site_var_bloodsapa <- apply(pmeth_blood_sapajus_cov5p75_common, 1, function(x) mean(x,na.rm = TRUE))
pmeth_blood_sapajus_cov5p75_common_var <-  pmeth_blood_sapajus_cov5p75_common[site_var_bloodsapa > 0.1 & site_var_bloodsapa <0.9,]

site_var_urinesapa <- apply(pmeth_urine_sapajus_cov5p75_common, 1, function(x) mean(x,na.rm = TRUE))
pmeth_urine_sapajus_cov5p75_common_var <-  pmeth_urine_sapajus_cov5p75_common[site_var_urinesapa > 0.1 & site_var_urinesapa <0.9,]

# ncp_best = estim_ncpPCA(t(pmeth_fecal_sapajus_cov5p75_common_var[1:5000,]), scale = TRUE, ncp.min=0, ncp.max=5, method.cv = "Kfold")$ncp #ncp_best=2
set.seed(4570)
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

# Assemble back the three matrices
pmeth_all_cov5p75_common <- cbind(pmeth_fecal_cebus_cov5p75_common_var_imp_f,cbind(pmeth_fecal_sapajus_cov5p75_common_var_imp_f,cbind(pmeth_urine_sapajus_cov5p75_common_var_imp_f,pmeth_blood_sapajus_cov5p75_common_var_imp_f)))

# Simplify metadata for glmnet
meta <- metadata_clean %>% dplyr::select(sample_name,monkey_name,sample_type,age)

meta <- meta[match(colnames(pmeth_all_cov5p75_common), meta$sample_name),]

identical(meta$sample_name, colnames(pmeth_all_cov5p75_common))

# saveRDS(pmeth_all_cov5p75_common, "/path/to/data_for_clocks_allsamples.rds")

###########################################################################################
#
# Section 3: SINGLE CpGs GLMNET CLOCK WITH/WITHOUT YEO JONHSON AND WITH/WITHOUT AGE TRANSFORMATION
#
# The code is supposed to be run 4 times with all combinations of YJ and age transformation possible.
# Results are exported after each run.
###########################################################################################

##### The data exported will be loaded in turn to generate different versions of the clock following different pre-processing steps.

pmeth_all_cov5p75_common <- readRDS("/path/to/data_for_clocks_allsamples.rds")

meta <- metadata_clean %>% dplyr::select(sample_name,monkey_name,sample_type,age)

meta <- meta[match(colnames(pmeth_all_cov5p75_common), meta$sample_name),]

identical(meta$sample_name, colnames(pmeth_all_cov5p75_common))

# Normalization
yeo_johnson_transformation <- function(x){  
  result <- yeojohnson(x)
  resultt <- result$x.t
  return(resultt)
}

# Mind that the function is returning a transposed matrix
pmeth_all_cov5p75_common_Norm <- apply(pmeth_all_cov5p75_common,1, yeo_johnson_transformation)

identical(rownames(pmeth_all_cov5p75_common_Norm), meta$sample_name)

# Transform age according to age at sexual maturity
# based on Horvath 2013 https://genomebiology.biomedcentral.com/articles/10.1186/gb-2013-14-10-r115#MOESM2
Fage = Vectorize(function(x){
  if (is.na(x) | is.na(maturity)) {return(NA)}
  k <- 1
  y <- 0
  if (x < maturity) {y = log(x+k)-log(maturity+k)}
  else {y = (x-maturity)/(maturity+k)}
  return(y)
})

## Inverse log linear transformation
F.inverse= Vectorize(function(y) {
  if (is.na(y) | is.na(maturity)) {return(NA)}
  
  k <- 1
  x <- 0
  if (y < 0) {x = (maturity+k)*exp(y)-k}
  else {x = (maturity+k)*y+maturity}
  return(x)
})

# Set age at maturity for the transformation
maturity=5

monkey_ids <- unique(meta$monkey_name)

#### FUNCTION Define the leave one out glmnet regression function
leave_one_out_glmnet <- function(mnky_id, pmeth) {
  
  if(ncol(pmeth)<nrow(pmeth))stop("Fewer sites than samples ! Does the matrix need to be transposed?")
  
  monkeys_samples_indices <- which(meta$monkey_name == mnky_id)
  sample_names <- meta$sample_name[monkeys_samples_indices]
  
  # Indices for train samples
  train_samples_indices <- which(!(rownames(pmeth) %in% sample_names))
  
  # Data partition
  train <- pmeth[train_samples_indices, ]
  test <- pmeth[monkeys_samples_indices, ]
  
  # Age vectors
  trainage <- meta$age[train_samples_indices]
  testage <- meta$age[monkeys_samples_indices]
  
  # Initialize a list to store results for different alpha values
  individual_results <- list()
  
  for (alph in seq(0, 1, by = 0.1)) {
    
    # Elastic net model
    if(exists("Fage", envir = .GlobalEnv) & exists("F.inverse", envir = .GlobalEnv)){
      # if the age transformation is loaded in the environment use it
      model <- cv.glmnet(train, Fage(trainage), type.measure = "mse", nfolds = 10, alpha = alph, standardize = FALSE) # with age transformation
    } else {
      # otherwise set the model without
      model <- cv.glmnet(train, trainage, type.measure = "mse", nfolds = 10, alpha = alph, standardize = FALSE)
    }
    
    # Predict age using the test sample from parameters that minimized MSE during internal CV
    if(exists("Fage", envir = .GlobalEnv) & exists("F.inverse", envir = .GlobalEnv)){
      predicted <- as.numeric(F.inverse(predict(model, newx = test, type = "response", s = "lambda.min"))) # with age transformation
    } else {
      predicted <- as.numeric(predict(model, newx = test, type = "response", s = "lambda.min"))
    }
    
    predicted_SAMP <- data.frame(sample_name = sample_names, age = testage, predicted_age = predicted)
    
    individual_results[[paste("alpha", alph)]] <- predicted_SAMP
  }
  
  return(individual_results)
}

set.seed(5500)

# Parallelize the leave one out glmnet regression over individuals ### CAREFUL !!! ALWAYS CHECK THAT SAMPLES=ROWS, when necessary transpose the matrix in the call below (e.g. pmeth = t(pmeth_all_cov5p75_common)).
predicted_out_list <- mclapply(monkey_ids, leave_one_out_glmnet, pmeth = pmeth_all_cov5p75_common_Norm, mc.cores = 16)

# Name the list elements by monkey IDs
names(predicted_out_list) <- monkey_ids

# Function to combine results for a given alpha across all monkeys
combine_results_for_alpha <- function(alpha_val) {
  alpha_key <- paste("alpha", alpha_val)
  combined_results <- do.call(rbind, lapply(predicted_out_list, `[[`, alpha_key))
  return(combined_results)
}

# Combine results for each alpha value
combined_results_list <- lapply(seq(0, 1, by = 0.1), combine_results_for_alpha)
names(combined_results_list) <- c(paste("alpha",seq(0, 1, by = 0.1)))

# saveRDS(combined_results_list,"/path/to/All_samples_Clockprediction_IDENTIFIER_PROCESSING.rds")

# --> The best model is uses Yeo-Johnson normalization and Transformation of age according to age at sexual maturity 

#### Load back results ####

# combined_results_list <- readRDS("/path/to/All_samples_Clockprediction_IDENTIFIER_PROCESSING.rds")

## Calculate alpha minimizing MSE
mse_values_clock = lapply(combined_results_list, function(x){
  mean((x$predicted_age - x$age)^2)
})
min(unlist(mse_values_clock))
View(mse_values_clock)

# predicted_out_min <- combined_results_list[["alpha 0.4"]] #classic
# predicted_out_min <- combined_results_list[["alpha 0.5"]] #age transfo
# predicted_out_min <- combined_results_list[["alpha 0.3"]] #Norm
predicted_out_min <- combined_results_list[["alpha 0.8"]] #NormAgetransfo

# Apply to the list of prediction to calculate summary statistics     
summ.stats_clock <- corr_and_mae(predicted_out_min)
summ.stats_clock

# add metadata
inspect_dev = merge(predicted_out_min, metadata_clean %>% dplyr::select(sample_name, batch, monkey_name, species, sex, sample_type)) #clock predictions

# change name of pearson for plot
names(summ.stats_clock) <- c("r","MAE")

# Fig 3A
ggplot(inspect_dev, aes(x=age, y=predicted_age))+
  geom_point(aes(col=sample_type, shape=paste(sex,species," ")),size=1.8) +
  scale_color_manual(values=c("darkred","black","orange"))+
  scale_shape_manual(values=c(16,1,17,24))+
  xlim(min(cbind(min(inspect_dev$age),min(inspect_dev$predicted_age))),max(cbind(max(inspect_dev$age),max(inspect_dev$predicted_age))))+
  ylim(min(cbind(min(inspect_dev$age),min(inspect_dev$predicted_age))),max(cbind(max(inspect_dev$age),max(inspect_dev$predicted_age))))+
  labs(x="chronological age", y="predicted age", col="sample source", shape="sex, species")+
  annotate("text",
           x = min(inspect_dev$age) + 33, y = 12,
           vjust=1, size=4.5,
           label = paste0("N=",nrow(inspect_dev),"\n",
                         paste(names(summ.stats_clock),summ.stats_clock, sep = "=",collapse = "\n")))+
  geom_abline(linetype="dashed")+
  stat_smooth(method = "lm", span = 1.1, se = FALSE, color="black")+
  theme_bw()+
  theme(axis.title = element_text(size = 14, color = "black"),
        axis.text = element_text(size = 14, color = "black"),
        legend.text = element_text(size = 12, color = "black"),
        legend.title = element_text(size = 12, color = "black"),
        legend.position = c(0.005, 1),
        legend.justification = c("left","top"),
        legend.background = element_blank(),
        legend.key = element_rect(fill=NA, color=NA),
        legend.spacing = unit(0, "lines"))

###########################################################################################
#
# Optional: CLOCK WITH PRESELECTION OR FILTERING OUT SITES ASSOCIATED WITH THIRD COVARIATES
#
###########################################################################################

## The procedure is exactly the same as shown above, with the addition of two alternative steps.

### OPTION 1: Filter CpG sites which show low-moderate correlation with age ## Careful works on CpGs as rows
corr_vals_strong <- apply(pmeth_all_cov5p75_common, 1, function(y) {
  na_omit_indices <- which(!is.na(y))
  filtered_y <- y[na_omit_indices]
  filtered_age <- meta$age[na_omit_indices]

  test_result <- cor.test(filtered_y, filtered_age, method = "spearman")
  if (!is.na(test_result$estimate) && abs(test_result$estimate) > 0.2) {
    return(test_result$estimate)
  } else {
    return(NA)
  }
})

# Filter rows ## Careful I apply to columns because I am using imputed data from missMDA
pmeth_all_cov5p75_common_filtered_strong <- pmeth_all_cov5p75_common[!is.na(corr_vals_strong),]

## From there on, Yeo-Johnson normalization and transformation of age may or may not be applied, and the glmnet function used.

#### OPTION 2: Remove sites associated with species or sample type
fit_sapajus <- read.table("/path/to/pqlseq_sapajus_sampletype_blood.txt", sep = "\t", header = TRUE) 
fit_fecal_both <- read.table("/path/to/pqlseq_fecal_both_species.txt", sep = "\t", header = TRUE) 
fit_cebus_sex <- read.table("/path/to/pqlseq_cebus_sex.txt", sep = "\t", header = TRUE) 

fit_sapajus_include <- fit_sapajus %>% filter(qval > 0.05) %>% pull(site)
fit_fecal_both_include <- fit_fecal_both %>% filter(qval > 0.05) %>% pull(site)
fit_cebus_sex_include <- fit_cebus_sex %>% filter(qval > 0.05) %>% pull(site)

site_included <- intersect(intersect(fit_sapajus_include, fit_fecal_both_include),fit_cebus_sex_include)

# Remove sites exhibiting a significant association with one of the covariate tested
pmeth_all_cov5p75_common_nobias <- pmeth_all_cov5p75_common[site_included,]

## From there on, Yeo-Johnson normalization and transformation of age may or may not be applied, and the glmnet function used.

###########################################################################################
#
# Section 5: Plotting the MAE and Corr coefficients for all clocks
#
###########################################################################################

# Load all the clock predictions
predicted_out_1 <- readRDS()
predicted_out_2 <- readRDS()
# ....

# Extract for each dataset, the prediction made with alpha minimizing MSE
predicted_out_1_min <- predicted_out_1[["alpha 0"]] #example
# ....

# Compute model performance
summ.stats_predicted_out_1 <- corr_and_mae(predicted_out_1)
# ....

# Gather model performance
result_elastic <- as.data.frame(rbind(summ.stats_predicted_out_1, #......))

result_elastic$model <- rownames(result_elastic)

result_elastic <- result_elastic %>%
  arrange(desc(-MAE)) %>% 
  mutate(model = factor(model, levels = unique(model)))
  
bu_pu_colors <- RColorBrewer::brewer.pal(6, "BuPu")

Cor_plot <- ggplot(result_elastic, aes(x=model,y=cor))+
  geom_bar(stat = "identity", fill = bu_pu_colors[5])+
  ylim(0,1)+
  labs(y="Pearson's correlation",x="")+
  theme_bw()+
  theme(axis.title = element_text(size=16, color="black"),
        axis.text.y =  element_text(size=16, color="black"),
        axis.text.x = element_text(size=14, color="black", angle = 55, vjust = 1, hjust=1))

MAE_plot <- ggplot(result_elastic, aes(x=model,y=MAE))+
  geom_bar(stat = "identity", fill = bu_pu_colors[3])+
  labs(x="", y="MAE (years)")+
  ylim(0,4)+
  theme_bw()+
  theme(axis.title = element_text(size=16, color="black"),
        axis.text.y =  element_text(size=16, color="black"),
        axis.text.x = element_text(size=14, color="black", angle=55, vjust = 1, hjust=1))

result_plots <- ggpubr::ggarrange(MAE_plot, Cor_plot, nrow=1)
result_plots

##### After determining which model performs best, express the error on repeated individuals 

inspect_dev = merge(predicted_out_1,metadata_clean %>% dplyr::select(sample_name, batch, monkey_name,species, sex, sample_type))

# Express the error
inspect_dev$age_error = abs(inspect_dev$predicted_age - inspect_dev$age)

# Error on repeated individuals
inspect_dev %>%
  group_by(monkey_name) %>% 
  summarize(n=length(unique(age))) %>% filter(n>1) %>% pull(monkey_name) -> repeated_names

# sd in years between repeated sampled individuals
inspect_dev %>%
  filter(monkey_name %in% repeated_names) %>% 
  group_by(monkey_name) %>% 
  summarize(sd_err= sd(age_error)) %>% 
  summarize(mean_err =mean(sd_err), sdr=sd(sd_err))

############# END #################################################
