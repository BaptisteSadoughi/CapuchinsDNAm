####### DNAm in capuchins
## Exploring the average methylation at different genomic regions associated with
## specific tissues in the HumanMethylationAtlas.
## Aim: compare intestinal epithelium cells to all other cell types in blood methylation profiles.

rm(list = ls())

library_list <- c("tidyverse","stringr","comethyl","GenomicRanges","bsseq","ggpubr")
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



#### Section 1: Fecal DNA methylation profiles are evocative of epithelial cell types

# Load cell and tissue type markers from the HumanMethylationAtlas Loyfer et al.2024 Supplementary Table S4C. List of 50286 cell type-specific unmethylated markers (top 1000, hg19).
human_top1000markers = read.csv("/path/to/top1000markers_HumanMethylationAtlas_Loyfer2024_SupTables.csv")[,1:16]

## intestinal epithelium
target_markers <- as.vector(c("Colon-Ep","Colon-Ep:Gastric-Ep:Small-Int-Ep","Colon-Ep:Small-Int-Ep","Gastric-Ep:Small-Int-Ep","Small-Int-Ep"))

## grabbing promoter markers for all tissues
target_markers_promoters = human_top1000markers %>%
  filter(Genomic_class %in% c("promoter-TSS","intro_promoter-TSS","Intergenic_promoter-TSS",
                              "TTS_promoter-TSS","exon_promoter-TSS",
                              "intron_promoter-TSS","exon_intron_promoter-TSS",
                              "Intergenic_intron_promoter-TSS"))

### Cebus imitator gtf
library(GenomicFeatures)
library(rtracklayer)

# Read the GTF file.
gtf_capu <- import('/path/to/Cebus_imitator.Cebus_imitator-1.0.113.gtf')

# Extract gene_id and gene_name.
gene_data <- data.frame(gene_id = elementMetadata(gtf_capu)$gene_id,
                        Gene = elementMetadata(gtf_capu)$gene_name,
                        stringsAsFactors = FALSE)
gene_data = gene_data %>% distinct_all()

# Create granges object
txdb = makeTxDbFromGFF('/path/to/Cebus_imitator.Cebus_imitator-1.0.113.gtf')
txdb

# pull out genes/promoters/etc. (any annotation you are interested in)
genes_ = genes(txdb)
promoters_=promoters(genes_,upstream=2000,downstream=200)

promoters_df <- as.data.frame(promoters_) %>% 
  mutate(gene_id = elementMetadata(promoters_)$gene_id)

promoters_df <- merge(promoters_df,gene_data, all.x=TRUE)

promoters_df <- promoters_df %>% dplyr::rename(start_ci = start, end_ci = end)

target_markers_promoters_HuCi <- left_join(target_markers_promoters, promoters_df, by = "Gene")
 
target_markers_promoters_HuCi_oi <- target_markers_promoters_HuCi[!is.na(target_markers_promoters_HuCi$seqnames),]

## how many intestinal epithelium still available?
nrow(target_markers_promoters_HuCi_oi[target_markers_promoters_HuCi_oi$Type %in% target_markers,])


##### Filter data for coverage
data_ttms_fecal_cebus = data_ttms_fecal[,sampleNames(data_ttms_fecal) %in%
                                          metadata_fecal$sample_name[metadata_fecal$species=="Cebus"]]

data_ttms_fecal_cebus_cov5p75 <- filterCpGs(data_ttms_fecal_cebus, cov = 5, perSample = 0.75,
                                            verbose = FALSE, save = FALSE, file = NULL)

data_ttms_fecal_sapajus = data_ttms_fecal[,sampleNames(data_ttms_fecal) %in%  metadata_fecal$sample_name[metadata_fecal$species=="Sapajus"]]

data_ttms_fecal_sapajus_cov5p75 <- filterCpGs(data_ttms_fecal_sapajus, cov = 5, perSample = 0.75,
                                              verbose = FALSE, save = FALSE, file = NULL)

source("Functions_intersect_bsseq.R")

# Apply function to the bsseq objects
result_intersect_fecal_species <- extract_common_cpgs(
  bsseq_1 = data_ttms_fecal_cebus_cov5p75, 
  bsseq_2 = data_ttms_fecal_sapajus_cov5p75,
  bsseq_common_1_name = "fecal_cebus_common", 
  bsseq_common_2_name = "fecal_sapajus_common"
)


data_ttms_fecal_cebus_cov5p75_common <- result_intersect_fecal_species$fecal_cebus_common
data_ttms_fecal_sapajus_cov5p75_common <- result_intersect_fecal_species$fecal_sapajus_common

# Extract promoter coordinates
bed_promoters_oi <- target_markers_promoters_HuCi_oi %>% dplyr::select(seqnames,start_ci,end_ci) %>% 
  dplyr::rename(start=start_ci,end=end_ci)

### Remove duplicates
bed_promoters_oi = bed_promoters_oi %>% distinct_all()

# Calculate promoter total coverage and methylation
coverage_promoters_oi_cebus <- bsseq::getCoverage(data_ttms_fecal_cebus_cov5p75_common, type = "Cov",
                                            regions = bed_promoters_oi,
                                            what="perRegionTotal", withDimnames = TRUE)

coverage_promoters_oi_sapa <- bsseq::getCoverage(data_ttms_fecal_sapajus_cov5p75_common, type = "Cov",
                                                  regions = bed_promoters_oi,
                                                  what="perRegionTotal", withDimnames = TRUE)

coverage_promoters_oi <- cbind(coverage_promoters_oi_cebus,coverage_promoters_oi_sapa)

methcount_promoters_oi_cebus <- bsseq::getCoverage(data_ttms_fecal_cebus_cov5p75_common, type = "M",
                                             regions = bed_promoters_oi,
                                             what="perRegionTotal", withDimnames = TRUE)

methcount_promoters_oi_sapa <- bsseq::getCoverage(data_ttms_fecal_sapajus_cov5p75_common, type = "M",
                                                   regions = bed_promoters_oi,
                                                   what="perRegionTotal", withDimnames = TRUE)

methcount_promoters_oi <- cbind(methcount_promoters_oi_cebus,methcount_promoters_oi_sapa)

# Create matrix of percent methylation
pmeth_promoters_oi <- methcount_promoters_oi/coverage_promoters_oi

rownames(pmeth_promoters_oi) <- paste(bed_promoters_oi$seqnames, bed_promoters_oi$start, bed_promoters_oi$end, sep="_")

pmeth_promoters_oi <- pmeth_promoters_oi[!apply(is.na(pmeth_promoters_oi),1,all),]

### Then compare the average %methylation between promoters which are markers for intestinal and others.

mean_pmeth_promoters_oi <- rowMeans2(pmeth_promoters_oi, na.rm = TRUE)

markers_mean_pmethylation = data.frame(site = names(mean_pmeth_promoters_oi), 
                                       av_pmeth_promoter = as.numeric(mean_pmeth_promoters_oi))

markers_mean_pmethylation = target_markers_promoters_HuCi_oi %>%
  mutate(site = paste(.$seqnames, .$start_ci, .$end_ci, sep="_")) %>% 
  dplyr::select(Type,Target_meth,Background_meth,Genomic_class,Gene,seqnames,site) %>% 
  right_join(.,markers_mean_pmethylation) %>% 
  mutate(target_markers = as.factor(ifelse(.$Type %in% target_markers,"1","0")))

# write.csv(markers_mean_pmethylation, "/path/to/intestinal_hypomethylated_markers.csv", row.names = FALSE, quote = FALSE)

## how many capuchin promoter per cell type
markers_mean_pmethylation %>% group_by(site) %>% summarize(n_type = length(unique(Type))) %>% summarize(mean_t=mean(n_type), sd_t=sd(n_type))

# Fig 2A
ggplot(markers_mean_pmethylation, aes(x=target_markers, y=av_pmeth_promoter, fill=target_markers))+
  geom_boxplot(outlier.shape=NA, width=0.45)+
  geom_jitter(height = 0)+
  labs(x="",y="promoter mean percent methylation", fill="intestinal\nepithelium")+
  scale_fill_manual(values=c("lightblue","yellow"), labels = c("other cell types","intestinal epithelium"))+
  scale_x_discrete(labels = c("other cell types","intestinal epithelium"))+
  geom_bracket(xmin = "0", xmax = "1", y.position = 1.05, label = "*",
               color = "black", size = 1, inherit.aes=F) +
  theme_bw()+
  theme(axis.text = element_text(size=16, color = "black"),
                   axis.title = element_text(size=16, color = "black"),
                   legend.position="none")

### Calculate p-value for the difference in average methylation based on permutations

# Compute observed difference in means
obs_diff <- mean(markers_mean_pmethylation$av_pmeth_promoter[markers_mean_pmethylation$target_markers == 0], na.rm = TRUE) -
  mean(markers_mean_pmethylation$av_pmeth_promoter[markers_mean_pmethylation$target_markers == 1], na.rm = TRUE)

# Number of permutations
n_permutations <- 10000
perm_diffs <- numeric(n_permutations)

# Permutation test
for (i in 1:n_permutations) {
  permuted_labels <- sample(markers_mean_pmethylation$target_markers)  # Shuffle target_markers
  perm_diff <- mean(markers_mean_pmethylation$av_pmeth_promoter[permuted_labels == 0], na.rm = TRUE) -
    mean(markers_mean_pmethylation$av_pmeth_promoter[permuted_labels == 1], na.rm = TRUE)
  perm_diffs[i] <- perm_diff
}

# Compute one-tailed p-value: probability of getting as extreme or more extreme difference
p_value <- mean(perm_diffs >= obs_diff)

cat("Observed difference:", obs_diff, "\n")
cat("P-value:", p_value, "\n")

######## END OF THE CODE #########
