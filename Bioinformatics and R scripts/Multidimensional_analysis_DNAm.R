####### DNAm in capuchin research project
## This script explores all samples together on PCA to assess
## the impact of batch effects, sex, species, and age on
## overall methylation profiles.

rm(list = ls())

library_list <- c("tidyverse","stringr","comethyl","GenomicRanges","bsseq","ggpubr","emmeans")
lapply(library_list, require, character.only=TRUE)

setwd("/home/sbaptis7/Capuchins_feces_DNA/")


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


data_ttms_fecal_cebus = data_ttms_fecal[,sampleNames(data_ttms_fecal) %in%
                                          metadata_fecal$sample_name[metadata_fecal$species=="Cebus" & metadata_fecal$estimated_dob=="n"]]
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

# Assemble back the three matrices
pmeth_all_cov5p75_common <- cbind(pmeth_blood_sapajus_cov5p75_common,cbind(pmeth_fecal_cebus_cov5p75_common,cbind(pmeth_fecal_sapajus_cov5p75_common,pmeth_urine_sapajus_cov5p75_common)))

# Select fully covered sites for the PCA
pmeth_all_cov5p75_common_var_full <- pmeth_all_cov5p75_common[complete.cases(pmeth_all_cov5p75_common),]

# Run PCA on fecal samples
pca_all <- prcomp(t(pmeth_all_cov5p75_common_var_full), scale. = TRUE)

# Extract principal component scores
pca_all_scores <- as.data.frame(pca_all$x)

# Variance explained by the PCs
variance_explained_all <- pca_all$sdev^2
prop_var_explained_all <- variance_explained_all / sum(variance_explained_all)
prop_var_explained_all
  
metadata_pca <- rbind(metadata_blood,metadata_fecal[metadata_fecal$estimated_dob=="n",],metadata_urine)
metadata_pca <- metadata_pca[match(rownames(pca_all_scores),metadata_pca$sample_name),]
  
identical(metadata_pca$sample_name,rownames(pca_all_scores))

# Fig. 2B
ggplot(metadata_pca, aes(x=pca_all_scores[,1], y=pca_all_scores[,2], fill = species, shape = sample_type)) +
  geom_point(size=3.5) +
  scale_shape_manual(values=c(21, 23,24)) +
  scale_fill_manual(values = c("#9370DB","lightblue"))+
  guides(fill = guide_legend(override.aes = list(shape = c(21,21))),
         # shape = guide_legend(override.aes = list(fill="lightblue")),
         )+
  labs(x=paste0("PC1 (",round(prop_var_explained_all[1]*100),"% variance explained)"),
       y=paste0("PC2 (",round(prop_var_explained_all[2]*100),"% variance explained)"),
       #caption = paste0("PCA performs on a subset of ", nrow(pmeth_all_cov5p75_common_var_full), " variable sites."),
       shape="sample source")+
  theme_bw() +
  theme(axis.text = element_text(size = 14, color="black"),
        axis.title = element_text(size = 14, color="black"),
        legend.background = element_blank())

###### EXPLORE TECHNICAL CORRELATES OF PC1 AND PC2

#### Model of covariates on PCs

metadata_pca
metadata_pca_model = cbind(metadata_pca, data.frame("PC1"=pca_all_scores[,1],"PC2"=pca_all_scores[,2]))

metadata_pca_model <- metadata_pca_model %>% dplyr::select(-c(date_of_birth,date_of_collection))
metadata_pca_model$batch <- as.factor(metadata_pca_model$batch)
metadata_pca_model$species <- as.factor(metadata_pca_model$species)
metadata_pca_model$sample_type <- as.factor(metadata_pca_model$sample_type)
metadata_pca_model$sex <- as.factor(metadata_pca_model$sex)
metadata_pca_model$z.age <- as.vector(scale(metadata_pca_model$age))
metadata_pca_model$z.map_efficiency_perc <- as.vector(scale(metadata_pca_model$map_efficiency_perc))
metadata_pca_model$z.mean_pmeth <- as.vector(scale(metadata_pca_model$mean_pmeth))
metadata_pca_model$z.perc_chg_meth <- as.vector(scale(metadata_pca_model$perc_chg_meth))
metadata_pca_model$z.perc_chh_meth <- as.vector(scale(metadata_pca_model$perc_chh_meth))

# Drop CHH because it correlates strongly with CHG and display lower variability.
# cor.test(metadata_pca_model$perc_chg_meth, metadata_pca_model$perc_chh_meth, method = "spearman")


###### PC1 model comparison
model_all_pc1 <- lm(PC1 ~ batch + species + sex + sample_type + z.age + z.map_efficiency_perc + z.mean_pmeth + z.perc_chg_meth, data = metadata_pca_plot_model, na.action = na.fail)
car::vif(model_all_pc1) #Diagnostic of variance inflation factors

### Impossible to have batch and the continuous covariates together!

# Model using continuous covariates for batch effects

model_all_pc1_CovTech <- lm(PC1 ~ species + sex + sample_type + z.age + z.map_efficiency_perc + z.mean_pmeth + z.perc_chg_meth, data = metadata_pca_plot_model, na.action = na.fail)
car::vif(model_all_pc1_CovTech)

# Perform model selection using dredge

all_models_PC1_CovTech <- MuMIn::dredge(model_all_pc1_CovTech, rank = "AICc")
print(all_models_PC1_CovTech)
all_models_PC1_CovTech_df = as.data.frame(all_models_PC1_CovTech %>% filter(delta < 10))

# Model using categorical batch effects

model_all_pc1_Batch <- lm(PC1 ~ species + sex + sample_type + z.age + batch, data = metadata_pca_plot_model, na.action = na.fail)
car::vif(model_all_pc1_Batch)

# Perform model selection using dredge

all_models_PC1_Batch <- MuMIn::dredge(model_all_pc1_Batch, rank = "AICc")
print(all_models_PC1_Batch)
all_models_PC1_Batch_df = as.data.frame(all_models_PC1_Batch %>% filter(delta < 10))

###### PC2 model comparison

# Model using continuous covariates for batch effects

model_all_pc2_CovTech <- lm(PC2 ~ species + sex + sample_type + z.age + z.map_efficiency_perc + z.mean_pmeth + z.perc_chg_meth, data = metadata_pca_plot_model, na.action = na.fail)
car::vif(model_all_pc2_CovTech)

# Perform model selection using dredge
all_models_PC2_CovTech <- MuMIn::dredge(model_all_pc2_CovTech, rank = "AICc")
print(all_models_PC2_CovTech)
all_models_PC2_CovTech = as.data.frame(all_models_PC2_CovTech %>% filter(delta < 10))

# Model using categorical batch effects
model_all_pc2_Batch <- lm(PC2 ~ species + sex + sample_type + z.age + batch, data = metadata_pca_plot_model, na.action = na.fail)
car::vif(model_all_pc2_Batch)

# Perform model selection using dredge
all_models_PC2_Batch <- MuMIn::dredge(model_all_pc2_Batch, rank = "AICc")
print(all_models_PC2_Batch)
all_models_PC2_Batch_df = as.data.frame(all_models_PC2_Batch %>% filter(delta < 10))

### For plotting purposes, run one model for PC1 and PC2 to extract a pvalue

model_all_pc1_plot <- lm(PC1 ~ species + sex + sample_type + z.age + batch, data = metadata_pca_plot_model)
pval_pc1 <- data.frame(predictor = names(summary(model_all_pc1_plot)$coeff[,4]),pvalue=summary(model_all_pc1_plot)$coeff[,4])

# Perform pairwise comparisons with Tukey adjustment and plot

######### PC1
## Batch
posthoc_pc1_batch <- emmeans(model_all_pc1_plot, pairwise ~ batch, adjust = "tukey")

ggplot(metadata_pca,aes(x=batch, y=pca_all_scores[,1]))+
  geom_violin(trim = F)+
  geom_boxplot(width = 0.1, outlier.shape = NA, fill = "white", color = "black") +
  labs(y=paste0("PC1 (",round(prop_var_explained_all[1]*100),"% variance explained)"))+
  geom_bracket(xmin = "b", xmax = "c", y.position = max(pca_all_scores[, 1]) + 22, label = "**",
              color = "black", size = 1) +
  geom_text(x = mean(c( which(metadata_pca$batch == "b"), which(metadata_pca$batch == "c"))),
            y = max(pca_all_scores[, 1]) + 22, label = "", vjust = -1, color = "red")+
  theme_bw()+
  theme(axis.title = element_text(size=15, color="black"),axis.text = element_text(size=15, color="black"))

## Species
ggplot(metadata_pca,aes(x=species, y=pca_all_scores[,1]))+
  geom_violin(trim = F)+
  geom_boxplot(width = 0.1, outlier.shape = NA, fill = "white", color = "black") +
  labs(y=paste0("PC1 (",round(prop_var_explained_all[1]*100),"% variance explained)"),x="")+
  geom_bracket(xmin = "Cebus", xmax = "Sapajus", y.position = max(pca_all_scores[, 1]) + 15, label = "***",
               color = "black", size = 1) +
  geom_text(x = mean(c( which(metadata_pca$species == "Cebus"), which(metadata_pca$species == "Sapajus"))),
            y = max(pca_all_scores[, 1]) + 15, label = "", vjust = -1)+
  theme_bw()+
  theme(axis.title = element_text(size=15, color="black"),axis.text = element_text(size=15, color="black"))

######### PC2
model_all_pc2_plot <- lm(PC2 ~ species + sex + sample_type + z.age + batch, data = metadata_pca_plot_model)
pval_pc2 <- data.frame(predictor = names(summary(model_all_pc2_plot)$coeff[,4]),pvalue=summary(model_all_pc2_plot)$coeff[,4])

## Batch
posthoc_pc2_batch <- emmeans(model_all_pc2_plot, pairwise ~ batch, adjust = "tukey")

ggplot(metadata_pca,aes(x=batch, y=pca_all_scores[,2]))+
  geom_violin(trim = F)+
  geom_boxplot(width = 0.1, outlier.shape = NA, fill = "white", color = "black") +
  labs(y=paste0("PC2 (",round(prop_var_explained_all[2]*100),"% variance explained)"))+
  geom_bracket(xmin = "b", xmax = "c", y.position = max(pca_all_scores[, 2]) + 22, label = "*",
               color = "black", size = 1) +
  geom_text(x = mean(c( which(metadata_pca$batch == "b"), which(metadata_pca$batch == "c"))),
            y = max(pca_all_scores[, 2]) + 22, label = "", vjust = -1)+
  geom_bracket(xmin = "a", xmax = "d", y.position = max(pca_all_scores[, 2]) + 30, label = "*",
               color = "black", size = 1) +
  geom_text(x = mean(c( which(metadata_pca$batch == "a"), which(metadata_pca$batch == "d"))),
            y = max(pca_all_scores[, 2]) + 30, label = "", vjust = -1)+
  theme_bw()+
  theme(axis.title = element_text(size=15, color="black"),axis.text = element_text(size=15, color="black"))

## Sample type
posthoc_pc2_sample_type <- emmeans(model_all_pc2_plot, pairwise ~ sample_type, adjust = "tukey")

ggplot(metadata_pca,aes(x=sample_type, y=pca_all_scores[,2]))+
  geom_violin(trim = F)+
  geom_boxplot(width = 0.1, outlier.shape = NA, fill = "white", color = "black") +
  labs(y=paste0("PC2 (",round(prop_var_explained_all[2]*100),"% variance explained)"),x="")+
  geom_bracket(xmin = "blood", xmax = "fecal", y.position = max(pca_all_scores[, 2]) + 22, label = "**",
               color = "black", size = 1) +
  geom_text(x = mean(c( which(metadata_pca$sample_type == "blood"), which(metadata_pca$sample_type == "fecal"))),
            y = max(pca_all_scores[, 2]) + 22, label = "", vjust = -1)+
  geom_bracket(xmin = "fecal", xmax = "urine", y.position = max(pca_all_scores[, 2]) + 30, label = "***",
               color = "black", size = 1) +
  geom_text(x = mean(c( which(metadata_pca$sample_type == "fecal"), which(metadata_pca$sample_type == "urine"))),
            y = max(pca_all_scores[, 2]) + 30, label = "", vjust = -1)+
  theme_bw()+
  theme(axis.title = element_text(size=15, color="black"),axis.text = element_text(size=15, color="black"))

############ END OF THE CODE ###############"
