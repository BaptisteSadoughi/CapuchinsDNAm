#### This script follows up on the PQLSEQ analysis for the effect of age in fecal samples to test whether results match those reported in a PanMammalianMethylation array EWAS by Lu et al. 2023

#### Extracting the most age-associated CpGs from the PanMammalianEWAS https://www.nature.com/articles/s43587-023-00462-6

rm(list = ls())

library_list <- c("tidyverse","stringr","comethyl","GenomicRanges","bsseq","GenomicFeatures","rtracklayer")
lapply(library_list, require, character.only=TRUE)

#Load the top 1000 CpGs in the pan mammalian EWAS for age Lu et al. 2023 Supplementary Table S6.1_TopMetaEWASinAllTissues

panmammal_ewas=readxl::read_xlsx("/path/to/S6.1_TopMetaEWASinAllTissues.xlsx")

# Assign direction of change with age to every site
panmammal_ewas_top1000 <- panmammal_ewas %>% mutate(sign = ifelse(Meta.Z<0,"neg","pos")) #%>% group_by(sign) %>% dplyr::slice(1:50)
panmammal_ewas_top1000 <- panmammal_ewas_top1000 %>% dplyr::select(Gene,GeneRegionID,annotation,sign)

# Capuchin GTF file.
gtf_capu <- import('/path/to/Cebus_imitator.Cebus_imitator-1.0.113.gtf')

# Extract gene_id and gene_name.
gene_data <- data.frame(gene_id = elementMetadata(gtf_capu)$gene_id,
                        Gene = elementMetadata(gtf_capu)$gene_name,
                        stringsAsFactors = FALSE)
gene_data = gene_data %>% distinct_all()

# Extract genes present in the top 1000 EWAS
selected_ewas_genes <- gene_data %>% filter(Gene %in% panmammal_ewas_top1000$Gene)

# Create granges object
txdb = makeTxDbFromGFF('/path/to/Cebus_imitator.Cebus_imitator-1.0.113.gtf')
txdb

#Build objects for promoters and exons
genes_=genes(txdb)
promoters_ = promoters(genes_,upstream=2000,downstream=200)
exons_by_gene <- exonsBy(txdb, by = "gene")

#Subset promoters and exons using selected_ewas_genes
promoters_target_ewas <- promoters_[selected_ewas_genes$gene_id]
exons_target_ewas <- exons_by_gene[names(exons_by_gene) %in% selected_ewas_genes$gene_id]

#Turn to dataframes and bind
promoters_df = as.data.frame(promoters_target_ewas) %>% 
  mutate(gene_id = elementMetadata(promoters_target_ewas)$gene_id,
         annotation = "promoter",
         exon_id = NA,
         exon_name = NA)

exon_ranges <- unlist(exons_target_ewas, use.names=FALSE)
exons_df <- as.data.frame(exon_ranges)
exons_df$gene_id <- rep(names(exons_target_ewas), elementNROWS(exons_target_ewas))
exons_df$annotation <- "exon"

selected_ewas_genes_capu_location <- rbind(promoters_df, exons_df)
rownames(selected_ewas_genes_capu_location) <- seq(1:nrow(selected_ewas_genes_capu_location))
selected_ewas_genes_capu_location = merge(selected_ewas_genes_capu_location, gene_data)

#Load result from pqlseq on fecal samples
pqlseq_fecal_both_age = read.table("/path/to/pqlseq_fecal_both_age.txt", header=TRUE, sep = "\t")

pqlseq_fecal_both_age = pqlseq_fecal_both_age %>% mutate(site_temp = site) %>% separate(site_temp, into = c("seqnames","start"), sep="_") %>% 
  mutate(start = as.numeric(start))

# Ensure the 'start' column is numeric
selected_ewas_genes_capu_location$start <- as.numeric(selected_ewas_genes_capu_location$start)
selected_ewas_genes_capu_location$end <- as.numeric(selected_ewas_genes_capu_location$end)

selected_ewas_genes_capu_location_gr <- with(selected_ewas_genes_capu_location, GRanges(seqnames, IRanges(start, end), strand = strand, gene_id))
pqlseq_gr <- with(pqlseq_fecal_both_age, GRanges(seqnames, IRanges(start, start), beta = beta, qval = qval, site = site))

overlaps_ <- findOverlaps(pqlseq_gr, selected_ewas_genes_capu_location_gr)

overlap_data <- cbind(pqlseq_fecal_both_age[queryHits(overlaps_), ], selected_ewas_genes_capu_location[subjectHits(overlaps_), ])

names(overlap_data)[14] <- "seqnames_annotation"
names(overlap_data)[15] <- "start_annotation"
names(overlap_data)[16] <- "end_annotation"

# CpGs were duplicated when overlapping several exons of a Gene.

# As we can only compare with the pan mammalian data at the "gene" level, we aggregate data.

overlap_data_gene <- overlap_data %>%
  dplyr::select(-c(seqnames_annotation, start_annotation, end_annotation, annotation, width, exon_id,exon_name)) %>% 
  distinct()

### Note: many "genes" in the panmammalian dataset are associated to both positive and negative effects (because they actually report single CpGs and not the entire gene trend).

### To avoid confusion, drop any gene that is not uniquely either increasing or decreasing with age.

Gene_ambiguous_pan <- panmammal_ewas_top1000 %>% group_by(Gene) %>% summarize(n=length(unique(sign))) %>% filter(n>1) %>% pull(Gene)

# Some of these might be still clearly overall changing positively or negatively.
proportion_pan <-panmammal_ewas_top1000 %>%
  filter(Gene %in% Gene_ambiguous_pan) %>% 
  mutate(beta_sign = case_when(
    sign == "pos" ~ "Positive",
    sign == "neg" ~ "Negative",
    TRUE ~ "Neutral"                                      # Handle beta = 0 if applicable
  )) %>%
  group_by(Gene, beta_sign) %>%                            # Group by Gene and beta_sign
  summarize(count = n(), .groups = 'drop') %>%            # Count occurrences
  group_by(Gene) %>%                                       # Regroup by Gene to calculate proportions
  mutate(proportion = count / sum(count)) %>%             # Compute proportion
  filter(beta_sign != "Neutral") %>%                      # Optionally exclude neutral (beta = 0) if not needed
  ungroup()

# Check the distribution of proportion positive per Gene in the panmammalian
proportion_pan %>%
  filter(beta_sign == "Positive") %>%
  dplyr::select(Gene, proportion_positive = proportion) %>% 
  ggplot(., aes(x = proportion_positive)) +
  geom_histogram(binwidth = 0.05, fill = "steelblue", color = "black", boundary = 0) +
  labs(
    # title = "Histogram of Proportion of Positive Betas per Gene",
    x = "proportion of positive effect sizes",
    y = "number of genes"
  ) +
  theme_minimal()+theme(axis.text = element_text(size=12, color="black"),
                        axis.title = element_text(size=12, color="black"))

## Based on this distribution, keep Genes with more than 75% of either positive or negative signs
keep_genes_clear_pan <- proportion_pan %>% filter(beta_sign == "Positive"& proportion >=0.75 |
                                                 beta_sign == "Positive"& proportion <=0.25) %>% pull(Gene) %>% unique()

genes_to_drop_pan <- setdiff(Gene_ambiguous_pan,keep_genes_clear_pan)

# Finally to avoid duplicating data, keep only the dominant sign for the panmammalian data (it is less problematic to keep all the betas from the capu even if a Gene is associated to both positive and negative, that's just what the data is. Yes it means one "pan sign" is associated to several betas, potentially opposite, but at least it avoids duplicating the betas and inflating the dataset).

# Step 1: Filter out genes in genes_to_drop_pan
filtered_data <- panmammal_ewas_top1000 %>%
  filter(!Gene %in% genes_to_drop_pan)

# Step 2: Identify and keep dominant sign entries for genes in keep_genes_clear_pan
dominant_sign_data <- filtered_data %>%
  filter(Gene %in% keep_genes_clear_pan) %>%
  group_by(Gene) %>%
  filter(sign == names(sort(table(sign), decreasing = TRUE)[1])) %>%
  ungroup()

# Step 3: Keep genes that are neither in genes_to_drop_pan nor in keep_genes_clear_pan unchanged
unique_sign_data <- filtered_data %>%
  filter(!Gene %in% keep_genes_clear_pan)

# Combine all the data
clear_gene_pan <- bind_rows(dominant_sign_data, unique_sign_data)

# Finally prepare data for plotting
overlap_data_unique = overlap_data_gene %>%
  left_join(., clear_gene_pan %>%
              dplyr::select(Gene,sign) %>%
              distinct(),
            join_by(Gene==Gene))

# Assign sign based on pan and on capu results, and drop Genes which could not be assigned a panmammalian sign
overlap_data_unique_plot <- overlap_data_unique %>%
  mutate(sign_pan = recode(sign, "neg" = "negative",
                           "pos"="positive"),
         sign_capu = ifelse(beta>0,"positive","negative")) %>% 
  filter(!is.na(sign_pan))

# Enrichment for consistent age-associated changes
overlap_capu_panmam <-fisher.test(table(overlap_data_unique_plot$sign_capu,overlap_data_unique_plot$sign_pan))
overlap_capu_panmam

## Plot the proportion of pan mammalian direction of change per binned capuchin effect sizes
summary(overlap_data_unique_plot$beta[overlap_data_unique_plot$beta<0.4]) #one clear outlier

overlap_data_unique_plot_sub <- overlap_data_unique_plot %>% filter(beta<0.4) #One clear outlier

overlap_data_unique_plot_sub$beta_bin <- cut(overlap_data_unique_plot_sub$beta, breaks = 10)

summary_df <- overlap_data_unique_plot_sub %>%
  group_by(beta_bin, sign_pan) %>%
  summarise(count = n()) %>%
  ungroup() %>%
  group_by(beta_bin) %>%
  mutate(proportion = count / sum(count))

ggplot(summary_df, aes(x = beta_bin, y = proportion, fill = sign_pan)) +
  geom_bar(stat = "identity") +
  labs(x="Binned age associated effect sizes in capuchins", y ="Proportion", fill="pan mammalian\ndirection of change\nwith age")+
  theme_minimal() +
  theme_bw()+theme(axis.text = element_text(size=14, color="black"),
                   axis.text.x = element_text(size=7, color="black"),
                   axis.title = element_text(size=14, color="black"),
                   legend.text = element_text(size=12, color="black"),
                   legend.title = element_text(size=12, color="black"))

########## END OF THE CODE ##############
