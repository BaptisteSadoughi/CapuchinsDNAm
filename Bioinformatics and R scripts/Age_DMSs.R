####### DNAm in capuchin research project
## This script investigates age-associated changes in DNA methylation levels.

rm(list = ls())

library_list <- c("tidyverse","stringr","comethyl","GenomicRanges","bsseq","ggrepel","GenomicFeatures","rtracklayer","msigdbr","fgsea","PQLseq")
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
#### Section 1: Age effects in fecal samples                             #
#                                                                        #
##########################################################################

# Restrict to individuals with known ages when filtering for coverage

data_ttms_fecal_cebus = data_ttms_fecal[,sampleNames(data_ttms_fecal) %in%
                                          metadata_fecal$sample_name[metadata_fecal$species=="Cebus" & metadata_fecal$estimated_dob=="n"]]
data_ttms_fecal_cebus_cov5p75 <- filterCpGs(data_ttms_fecal_cebus, cov = 5, perSample = 0.75,
                                            verbose = FALSE, save = FALSE, file = NULL)


data_ttms_fecal_sapajus = data_ttms_fecal[,sampleNames(data_ttms_fecal) %in%  metadata_fecal$sample_name[metadata_fecal$species=="Sapajus"]]
data_ttms_fecal_sapajus_cov5p75 <- filterCpGs(data_ttms_fecal_sapajus, cov = 5, perSample = 0.75,
                                              verbose = FALSE, save = FALSE, file = NULL)

# Intersect datasets
source("Functions_intersect_bsseq.R") #it may be necessary to adapt the function to pass the desired number of bsseq objects

# Apply function to the bsseq objects
result_intersect_fecal_species <- extract_common_cpgs(
  bsseq_1 = data_ttms_fecal_cebus_cov5p75, 
  bsseq_2 = data_ttms_fecal_sapajus_cov5p75,
  bsseq_common_1_name = "fecal_cebus_common", 
  bsseq_common_2_name = "fecal_sapajus_common"
)

data_ttms_fecal_cebus_cov5p75_common <- result_intersect_fecal_species$fecal_cebus_common
data_ttms_fecal_sapajus_cov5p75_common <- result_intersect_fecal_species$fecal_sapajus_common

data_ttms_fecal_cov5p75 <- combine(data_ttms_fecal_cebus_cov5p75_common, data_ttms_fecal_sapajus_cov5p75_common)

# Extract the coordinate of CpG sites
cpg_loc_fecal_cov5p75 <- granges(data_ttms_fecal_cov5p75)
cpg_loc_fecal_cov5p75 <- paste(seqnames(cpg_loc_fecal_cov5p75), start(cpg_loc_fecal_cov5p75), sep = "_")

cov_fecal_both_cov5p75 <- as.matrix(getCoverage(data_ttms_fecal_cov5p75, type = "Cov"))
meth_fecal_both_cov5p75 <- as.matrix(getCoverage(data_ttms_fecal_cov5p75, type = "M"))

rownames(cov_fecal_both_cov5p75) <- cpg_loc_fecal_cov5p75
rownames(meth_fecal_both_cov5p75) <- cpg_loc_fecal_cov5p75

# Find nonvariable sites
pmeth_fecal_both_cov5p75=meth_fecal_both_cov5p75/cov_fecal_both_cov5p75
rownames(pmeth_fecal_both_cov5p75) <- cpg_loc_fecal_cov5p75

site_var_ <- apply(pmeth_fecal_both_cov5p75, 1, function(x) mean(x,na.rm = TRUE))

# Subset to variable sites
cov_fecal_both_cov5p75_var <- cov_fecal_both_cov5p75[site_var_ > 0.1 & site_var_ <0.9,]
meth_fecal_both_cov5p75_var <- meth_fecal_both_cov5p75[site_var_ > 0.1 & site_var_ <0.9,]

identical(rownames(cov_fecal_both_cov5p75_var),rownames(meth_fecal_both_cov5p75_var))

metadata_fecal_pqlseq <- metadata_fecal %>% filter(estimated_dob=="n") %>% dplyr::select(sample_name, sex, batch, species, age, monkey_name)

# match metadata and methylation matrices
metadata_fecal_pqlseq <- metadata_fecal_pqlseq[match(colnames(cov_fecal_both_cov5p75_var),metadata_fecal_pqlseq$sample_name),]

age <- as.vector(metadata_fecal_pqlseq$age)

covariates_fecal_both <- model.matrix(~ sex + species + batch, data = metadata_fecal_pqlseq)
covariates_fecal_both <- covariates_fecal_both[,-which(colnames(covariates_fecal_both) %in%  c("(Intercept)","batchd"))]

# Create a mapping from monkey_name to sample_name
monkey_name_map <- metadata_fecal_pqlseq %>%
  dplyr::select(monkey_name, sample_name) %>%
  dplyr::group_by(monkey_name) %>%
  dplyr::summarise(sample_names = list(sample_name))

# Initialize the identity matrix
identity_matrix <- diag(nrow(metadata_fecal_pqlseq))
rownames(identity_matrix) <- metadata_fecal_pqlseq$sample_name
colnames(identity_matrix) <- metadata_fecal_pqlseq$sample_name

# Fill in the identity matrix based on repeated monkey_names
for (i in seq_along(monkey_name_map$sample_names)) {
  sample_names <- monkey_name_map$sample_names[[i]]
  sample_names <- unlist(sample_names)
  
  for (s1 in sample_names) {
    for (s2 in sample_names) {
      if (s1 != s2) {
        identity_matrix[s1, s2] <- 1
        identity_matrix[s2, s1] <- 1
      }
    }
  }
}

if (!identical(rownames(identity_matrix),colnames(cov_fecal_both_cov5p75_var))|
    !identical(rownames(identity_matrix),colnames(meth_fecal_both_cov5p75_var)) |
    !identical(colnames(identity_matrix),metadata_fecal_pqlseq$sample_name) |
    !identical(colnames(cov_fecal_both_cov5p75_var),metadata_fecal_pqlseq$sample_name)){
  stop("STOP there are some mismatches across matrices")
}

fit_fecal_both = PQLseq::pqlseq(
  RawCountDataSet=meth_fecal_both_cov5p75_var,
  Phenotypes=age,
  RelatednessMatrix=identity_matrix,
  Covariates = covariates_fecal_both,
  LibSize=cov_fecal_both_cov5p75_var,
  fit.model="BMM",
  numCore = 16)

fit_fecal_both$site <-rownames(fit_fecal_both)
fit_fecal_both <- fit_fecal_both %>% filter(converged==TRUE)
qresults_fecal_both <- qvalue::qvalue(fit_fecal_both$pvalue)
fit_fecal_both$qval <- qresults_fecal_both$qvalues
fit_fecal_both$padj <- stats::p.adjust(fit_fecal_both$pvalue, method = "BH")

# write.table(fit_fecal_both, "/path/to/pqlseq_fecal_both_age.txt", sep = "\t", row.names = FALSE, quote = FALSE)
fit_fecal_both <-  read.table("/path/to/pqlseq_fecal_both_age.txt", sep = "\t", h=T)

### Read Cebus imitator GTF file.
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

fit_fecal_both_intersect = fit_fecal_both %>% mutate(site_temp = site) %>% separate(site_temp, into = c("seqnames","start"), sep="_") %>% 
  mutate(start = as.numeric(start))

# Ensure the 'start' column is numeric in both dataframes
promoters_df$start <- as.numeric(promoters_df$start)
promoters_df$end <- as.numeric(promoters_df$end)

promoters_capu_location_gr <- with(promoters_df, GRanges(seqnames, IRanges(start, end), strand = strand, gene_id))
pqlseq_gr <- with(fit_fecal_both_intersect, GRanges(seqnames, IRanges(start, start), beta = beta, qval = qval, site = site))

overlaps_ <- findOverlaps(pqlseq_gr, promoters_capu_location_gr)

overlap_data <- cbind(fit_fecal_both_intersect[queryHits(overlaps_), ], promoters_df[subjectHits(overlaps_), ])

names(overlap_data)[14] <- "seqnames_annotation"
names(overlap_data)[15] <- "start_annotation"
names(overlap_data)[16] <- "end_annotation"


### Annotate with the panmammalian data Lu et al. 2023 Supplementary Table S6.1_TopMetaEWASinAllTissues

panmammal_ewas=readxl::read_xlsx("/path/to/S6.1_TopMetaEWASinAllTissues.xlsx")

# Select promoters
panmammal_ewas_prom <- panmammal_ewas[panmammal_ewas$annotation=="Promoter (<=1kb)",c("Meta.Z","Gene","Gene_mm10")]

# Extract the direction of age-associated change for each "Gene"
panmammal_ewas_prom <- panmammal_ewas_prom %>% group_by(Gene) %>% dplyr::slice(1) %>% ungroup()

# Intersect with capuchin dataset
overlap_data_panmamm <- overlap_data %>% left_join(., panmammal_ewas_prom)
overlap_data_panmamm=overlap_data_panmamm %>% dplyr::rename("MetaZ"="Meta.Z")

## Annotate the x strongest positive and negative Genes

# Create a dataframe for top 20 genes with positive and negative beta separately
top_genes <- overlap_data_panmamm %>%
  mutate(beta_sign = ifelse(beta > 0, "positive", "negative")) %>%
  arrange(beta_sign, desc(abs(beta))) %>%
  group_by(beta_sign, Gene) %>%
  filter(row_number() == 1 & !is.na(Gene)) %>%
  ungroup() %>%
  group_by(beta_sign) %>%
  filter(row_number() <= 10) %>% 
  ungroup() %>% 
  dplyr::select(-beta_sign)

# Add label variable to the entire dataframe based on top_genes
overlap_data_panmamm <- overlap_data_panmamm %>% 
  left_join(top_genes %>% dplyr::select(Gene,beta) %>% mutate(label = Gene), by = c("Gene", "beta"))

#Fig. 4A
ggplot(overlap_data_panmamm %>% mutate(alpha=ifelse(qval<0.1,"sig","nonsig")),
       aes(x=beta,y=-log10(qval), color=alpha, alpha=alpha))+
  geom_point(shape=21,fill=NA)+
  scale_color_manual(values=c("nonsig"="lightgrey","sig"="black"))+
  scale_alpha_manual(values=c("sig"=1,"nonsig"=0.1))+
  labs(x="age associated effect sizes")+
  geom_hline(yintercept = -log10(0.1), linetype="dashed")+
  # geom_text(aes(label=label), check_overlap = TRUE, size=3, hjust=1, vjust=1)+
  geom_label_repel(aes(label=label), box.padding = 0.35, point.padding = 0.5, segment.color = 'grey50', max.overlaps = Inf)+
  theme_bw()+
  theme(axis.text = element_text(size=14, color="black"),
        axis.title = element_text(size=14, color="black"),
        legend.position = "none")

################################################################################
#
# GENE SET ENRICHMENT ANALYSIS
#
################################################################################

# Load pqlseq output clean
fit_fecal_both <-  read.table("/path/to/pqlseq_fecal_both_age.txt", sep = "\t", h=T)

# pull out genes clean
genes_ = genes(txdb)

genes_df <- as.data.frame(genes_) %>% 
  mutate(gene_id = elementMetadata(genes_)$gene_id)

genes_df <- merge(genes_df,gene_data, all.x=TRUE)

fit_fecal_both_intersect = fit_fecal_both %>%
  mutate(site_temp = site) %>%
  separate(site_temp, into = c("seqnames","start"), sep="_") %>% 
  mutate(start = as.numeric(start))

# Ensure the 'start' column is numeric in both dataframes
genes_df$start <- as.numeric(genes_df$start)
genes_df$end <- as.numeric(genes_df$end)

genes_capu_location_gr <- with(genes_df, GRanges(seqnames, IRanges(start, end), strand = strand, gene_id))
pqlseq_gr <- with(fit_fecal_both_intersect, GRanges(seqnames, IRanges(start, start), beta = beta, qval = qval, site = site))

overlaps_genes <- findOverlaps(pqlseq_gr, genes_capu_location_gr)

overlap_data_genes <- cbind(fit_fecal_both_intersect[queryHits(overlaps_genes), ], genes_df[subjectHits(overlaps_genes), ])

names(overlap_data_genes)[14] <- "seqnames_annotation"
names(overlap_data_genes)[15] <- "start_annotation"
names(overlap_data_genes)[16] <- "end_annotation"

# Generate rank-ordered vector by pqlseq coefficient
age_genes<- overlap_data_genes %>% 
  dplyr::select(Gene, beta) %>%
  filter(!is.na(Gene)) %>%
  arrange(desc(beta))
  
# Assign the mean beta across all CpGs overlapping the same gene
age_genes_unique <- age_genes %>%
  group_by(Gene) %>%
  summarise(mean_beta = mean(beta), .groups = 'drop') %>% 
  arrange(desc(mean_beta)) #very important to reorder

age_genes2 <- age_genes_unique$mean_beta
names(age_genes2) = age_genes_unique$Gene


######## Generate gene ontology set

## 1) Biological Processes

go_set_bp = msigdbr(species = "Homo sapiens", subcollection = "GO:BP") #human

go_set_bp = split(x = go_set_bp$gene_symbol, f = go_set_bp$gs_name)

matched_genes <- intersect(names(age_genes2), unlist(go_set_bp))
length(matched_genes)  # Check how many genes match

#Enrichment for Hallmark set
age_gsea_bp<- fgsea(pathways = go_set_bp, 
                 stats    = age_genes2,
                 minSize  = 15,
                 maxSize  = 500,
                 nPermSimple = 10000,
                 eps = 0.0)

#Assign pathway variable as factor for plotting
age_gsea_bp <- age_gsea_bp %>%
  mutate_at(vars(pathway), as.factor) %>%
  mutate(pathway = fct_reorder(pathway, ES))

age_gsea_bp <- age_gsea_bp %>% 
  mutate(pathway = gsub("GOBP_", "", pathway),    # remove "GOBP_"
         pathway = gsub("_", " ", pathway), # replace "_" with " "
         pathway = tolower(pathway))        # convert to lower-case

print(age_gsea_bp %>% dplyr::select(-leadingEdge) %>% as.data.frame())

### Optional: Collapse pathways (!careful if the names of the pathways were edited above, it's necessary to run fgsea again).

set.seed(123)
age_gsea_bp<- fgsea(pathways = go_set_bp, 
                    stats    = age_genes2,
                    minSize  = 15,
                    maxSize  = 500,
                    nPermSimple = 10000,
                    eps = 0.0)

collapsed_age_gsea_bp <- collapsePathways(fgseaRes = age_gsea_bp[order(pval)][padj < 0.05],
                                          pathways = go_set_bp,
                                          stats = age_genes2)

main_age_gsea_bp <- age_gsea_bp[pathway %in% collapsed_age_gsea_bp$mainPathways][order(-NES),]

#Assign pathway variable as factor for plotting
main_age_gsea_bp <- main_age_gsea_bp %>%
  mutate_at(vars(pathway), as.factor) %>%
  mutate(pathway = fct_reorder(pathway, ES))

main_age_gsea_bp <- main_age_gsea_bp %>% 
  mutate(pathway = gsub("GOBP_", "", pathway),    # remove "GOBP_"
         pathway = gsub("_", " ", pathway), # replace "_" with " "
         pathway = tolower(pathway))        # convert to lower-case

# edit entries
main_age_gsea_bp <- main_age_gsea_bp %>% mutate(pathway = recode(pathway, "response to type i interferon" = "response to type I interferon",
                                                                 "neurotransmitter receptor localization to postsynaptic specialization membrane" = "neurotransmitter recep. loc. to postsynaptic spe. membrane",
                                                                 "regulation of tumor necrosis factor mediated signaling pathway"="regulation of tumor necrosis factor mediated signaling path."))

# Fig. 4B
main_go_bp_plot <- ggplot(main_age_gsea_bp %>%
                            arrange(desc(abs(NES))) %>% 
                            dplyr::slice(c(1:20)) %>% 
                            arrange(desc(NES)) %>% 
                            mutate(sign_nes = ifelse(NES>0,"Hyper","Hypo")),
                            aes(reorder(pathway, abs(NES)), abs(NES))) +
  geom_point(aes(size=size, fill=sign_nes),shape=21) +
  coord_flip() +
  labs(x="", y="normalized enrichment score", fill="Direction\nof change", size="set size") +
  theme_bw() +
  theme(axis.text.y = element_text(size=14, color="black"),
        axis.text.x = element_text(size=14, color="black"),
        axis.title.x = element_text(size=14, color="black"),
        legend.margin = margin(l=1),
        legend.text = element_text(size=14, color="black"),
        legend.title = element_text(size=14, color="black")
        )

### Similar procedure used for Cellular Component and Molecular Functions

# 2) Cellular Component

go_set_cc = msigdbr(species = "Homo sapiens", subcollection = "GO:CC") #human
go_set_cc = split(x = go_set_cc$gene_symbol, f = go_set_cc$gs_name)

## 3) Molecular Function

go_set_mf = msigdbr(species = "Homo sapiens", subcollection = "GO:MF") #human
go_set_mf = split(x = go_set_mf$gene_symbol, f = go_set_mf$gs_name)

############## END OF THE CODE ################
