####### DNAm in capuchin research project
## Build DNAm clock from fecal samples.

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

data_ttms_fecal_cebus = data_ttms_fecal[,sampleNames(data_ttms_fecal) %in%  metadata_fecal$sample_name[metadata_fecal$species=="Cebus"]]
data_ttms_fecal_cebus_cov5p75 <- filterCpGs(data_ttms_fecal_cebus, cov = 5, perSample = 0.75,
                                            verbose = FALSE, save = FALSE, file = NULL)

data_ttms_fecal_cebus_cov5p75@assays@data@listData[["pmeth"]]=getCoverage(data_ttms_fecal_cebus_cov5p75, type = "M")/getCoverage(data_ttms_fecal_cebus_cov5p75, type = "Cov")

pmeth_fecal_cebus_cov5p75 <- data_ttms_fecal_cebus_cov5p75@assays@data@listData$pmeth
pmeth_fecal_cebus_cov5p75 <- as.matrix(pmeth_fecal_cebus_cov5p75)

# Use the cpg coordinate for the methylation matrices
cpg_loc_cov5p75 <- granges(data_ttms_fecal_cebus_cov5p75)
cpg_loc_cov5p75 <- paste(seqnames(cpg_loc_cov5p75), start(cpg_loc_cov5p75), sep = "_")

rownames(pmeth_fecal_cebus_cov5p75)<- cpg_loc_cov5p75

site_var_ <- apply(pmeth_fecal_cebus_cov5p75, 1, function(x) mean(x,na.rm = TRUE))
pmeth_fecal_cebus_cov5p75_var <- pmeth_fecal_cebus_cov5p75[site_var_ > 0.1 & site_var_ <0.9,]

# Impute data
set.seed(126589)
pmeth_fecal_cebus_cov5p75_var_imp=t(imputePCA(t(pmeth_fecal_cebus_cov5p75_var), scale = TRUE, ncp = 2)$completeObs)

# Simplify metadata for glmnet
meta <- metadata_clean %>% dplyr::select(sample_name,monkey_name,sample_type,age,estimated_dob)

meta <- meta[match(colnames(pmeth_fecal_cebus_cov5p75_var_imp), meta$sample_name),]

identical(meta$sample_name, colnames(pmeth_fecal_cebus_cov5p75_var_imp))

# saveRDS(pmeth_fecal_cebus_cov5p75_var_imp, "/path/to/data_for_clocks_cebusfecal.rds")

###########################################################################################
#
# Section 3: SINGLE CpGs GLMNET CLOCK WITH YEO JONHSON AND AGE TRANSFORMATION
#
###########################################################################################

##### The data exported will be loaded in turn to generate different versions of the clock following different pre-processing steps.

pmeth_fecal_cebus_cov5p75_var_imp <- readRDS("/path/to/data_for_clocks_cebusfecal.rds")

meta <- metadata_clean %>% dplyr::select(sample_name,monkey_name,sample_type,age,estimated_dob)

meta <- meta[match(colnames(pmeth_fecal_cebus_cov5p75_var_imp), meta$sample_name),]

identical(meta$sample_name, colnames(pmeth_fecal_cebus_cov5p75_var_imp))

#### FUNCTION Define the leave one out glmnet regression function
leave_one_out_glmnet <- function(mnky_id, pmeth, meta) {
  
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


##### OPTIONAL: data normalization

yeo_johnson_transformation <- function(x){  
  result <- yeojohnson(x)
  resultt <- result$x.t
  return(resultt)
}

# Mind that the function is returning a transposed matrix
pmeth_fecal_cebus_cov5p75_var_imp_Norm <- apply(pmeth_fecal_cebus_cov5p75_var_imp,1, yeo_johnson_transformation)

identical(rownames(pmeth_fecal_cebus_cov5p75_var_imp_Norm), meta$sample_name)


##### OPTIONAL: Age Transformation

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


### Remove individuals with uncertain age
meta_known <- meta %>% filter(estimated_dob == "n")

# Filter non normalized data
# pmeth_fecal_cebus_cov5p75_var_imp_known <- pmeth_fecal_cebus_cov5p75_var_imp[,colnames(pmeth_fecal_cebus_cov5p75_var_imp) %in% meta_known$sample_name]

# Filter normalized data
pmeth_fecal_cebus_cov5p75_var_imp_known <- pmeth_fecal_cebus_cov5p75_var_imp_Norm[rownames(pmeth_fecal_cebus_cov5p75_var_imp_Norm) %in% meta_known$sample_name,]

# identical(colnames(pmeth_fecal_cebus_cov5p75_var_imp_known), meta_known$sample_name)
identical(rownames(pmeth_fecal_cebus_cov5p75_var_imp_known), meta_known$sample_name)

## Define list of individuals for LOOV
monkey_ids <- unique(meta_known$monkey_name)

set.seed(4560)

# Parallelize the leave one out glmnet regression over individuals ### CAREFUL !!! ALWAYS CHECK THAT SAMPLES=ROWS, when necessary transpose the matrix in the call below (e.g. pmeth = t(pmeth_all_cov5p75_common)).
predicted_out_list <- mclapply(monkey_ids, leave_one_out_glmnet, pmeth = pmeth_fecal_cebus_cov5p75_var_imp_known, meta=meta_known, mc.cores = 16)

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

# saveRDS(combined_results_list,"/path/to/Cebusfecal_Clockprediction_IDENTIFIER_PROCESSING.rds"")

## Calculate alpha minimizing MSE
mse_values_clock = lapply(combined_results_list, function(x){
  mean((x$predicted_age - x$age)^2)
})
min(unlist(mse_values_clock))
View(mse_values_clock)

# Extract data at alpha minimizing MSE
predicted_out_min <- combined_results_list[["alpha 0.9"]] #example

# Apply to the list of prediction to calculate summary statistics     
summ.stats_clock <- corr_and_mae(predicted_out_min)
summ.stats_clock

# add metadata
inspect_dev = merge(predicted_out_min, metadata_clean %>% dplyr::select(sample_name, batch, monkey_name, species, sex, sample_type))

# Fig 3B
ggplot(inspect_dev, aes(x=age, y=predicted_age))+
  geom_point(aes(shape=sex),size=2.0, color="#9370DB") +
  xlim(min(cbind(min(inspect_dev$age),min(inspect_dev$predicted_age))),max(cbind(max(inspect_dev$age),max(inspect_dev$predicted_age))))+
  ylim(min(cbind(min(inspect_dev$age),min(inspect_dev$predicted_age))),max(cbind(max(inspect_dev$age),max(inspect_dev$predicted_age))))+
  labs(x="chronological age", y="predicted age")+
  annotate("text",
           x = min(inspect_dev$age) + 5, y = max(inspect_dev$predicted_age) - 2,
           vjust=1,
           label = paste(names(summ.stats_clock), summ.stats_clock, sep = "=",collapse = "\n"))+
  geom_abline(linetype="dashed")+
  stat_smooth(method = "lm", span = 1.1, se = FALSE, color="black")+
  theme_bw()+
  theme(axis.title = element_text(size = 14, color = "black"),
        axis.text = element_text(size = 14, color = "black"))

###########################################################################################
#
# Optional: CLOCK WITH PRESELECTION
#
###########################################################################################

## The procedure is exactly the same as shown above, with the addition of two alternative steps.

### OPTION: Filter CpG sites which show low-moderate correlation with age ## Careful works on CpGs as rows
corr_vals_strong <- apply(pmeth_fecal_cebus_cov5p75_var_imp, 1, function(y) {
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
pmeth_fecal_cebus_cov5p75_var_imp_strong <- pmeth_fecal_cebus_cov5p75_var_imp[!is.na(corr_vals_strong),]

identical(colnames(pmeth_fecal_cebus_cov5p75_var_imp_strong), meta$sample_name)

## From there on, Yeo-Johnson normalization and transformation of age may or may not be applied, and the glmnet function used.

###########################################################################################
#
# Section 4: Plotting the MAE and Corr coefficients for all clocks
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

###########################################################################################
#
# Section 5: Predicting age in out samples using the best clock pre-processing steps
#
###########################################################################################

##### Cebus individuals with estimated age 

identical(meta$sample_name, colnames(pmeth_fecal_cebus_cov5p75_var_imp))

# filter for CpG sites which show low-moderate correlation with age ## Careful works on CpGs as rows
corr_vals_strong <- apply(pmeth_fecal_cebus_cov5p75_var_imp, 1, function(y) {
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
pmeth_fecal_cebus_cov5p75_var_imp_strong <- pmeth_fecal_cebus_cov5p75_var_imp[!is.na(corr_vals_strong),]

identical(colnames(pmeth_fecal_cebus_cov5p75_var_imp_strong), meta$sample_name)

test_samples_indices <- which(meta$estimated_dob == "y")
test_sample_names <- meta$sample_name[test_samples_indices]
  
# Indices for train samples
train_samples_indices <- which(!(colnames(pmeth_fecal_cebus_cov5p75_var_imp_strong) %in% test_sample_names))
  
# Data partition
train <- pmeth_fecal_cebus_cov5p75_var_imp_strong[,train_samples_indices]
test <- pmeth_fecal_cebus_cov5p75_var_imp_strong[,test_samples_indices]
  
# Age vectors
trainage <- meta$age[train_samples_indices]
testage <- meta$age[test_samples_indices]

set.seed(8800)
model <- cv.glmnet(t(train), trainage, type.measure = "mse", nfolds = 10, alpha=1, standardize = FALSE)

predicted <- as.numeric(predict(model, newx = t(test), type = "response", s = "lambda.min"))
    
predicted_ages <- data.frame(sample_name = test_sample_names, age = testage, predicted_age = predicted)

# Apply to the list of prediction to calculate summary statistics     
summ.stats_clock <- corr_and_mae(predicted_ages)
summ.stats_clock

############ END OF THE CODE ######################
