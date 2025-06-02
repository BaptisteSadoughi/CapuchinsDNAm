## Generating differentially methylation regions based on clustering of CpG sites exhibiting consistent age-associated changes
## in preparation for Transcription Factor Motif Enrichment Analysis.

# This procedure is adapted from https://github.com/mwatowich/Immune-gene-regulation-is-associated-with-age-and-environmental-adversity-in-a-nonhuman-primate/tree/main/scripts

rm(list=ls())

library_list <- c("tidyverse","stringr","GenomicRanges")
lapply(library_list, require, character.only=TRUE)

dnam_out <-  read.table("path/to/pqlseq_fecal_both_age.txt", sep = "\t", h=T)

# Prep DMR df for variable

# get the beta, se, stdbeta, pvalue, and qvalue from the site-level model, as well as the chromosome and start and stop positions of the CpG

dnam_dss_fmt <- dnam_out[,c("beta", "se_beta","pvalue","site","qval")]
rownames(dnam_dss_fmt) <- dnam_dss_fmt$site
dnam_dss_fmt$chr <- gsub("_.*","", dnam_dss_fmt$site)
dnam_dss_fmt$start <- gsub(".*_","", dnam_dss_fmt$site)
dnam_dss_fmt$pos <- dnam_dss_fmt$start
dnam_dss_fmt <- dnam_dss_fmt %>% dplyr::select(-site)
colnames(dnam_dss_fmt) <- c("beta", "se_beta", "raw_pval","pval","chr","start","pos")
dnam_dss_fmt$pos<-as.integer(dnam_dss_fmt$pos)
dnam_dss_fmt$start<-as.integer(dnam_dss_fmt$start)


### Call DMR function ---------
dmr_function <- function(input_df, fdr_cutoff_strict=0.05, fdr_cutoff_relaxed=0.1, fdr_col_name="pval", effect_size_name="beta", minCG=4, sep=1000, prop_sig=0.5) {
  input_df=input_df[order(input_df$chr,input_df$pos),]
  sig_df=input_df[input_df[,fdr_col_name]<fdr_cutoff_relaxed,]
  pos=sig_df$pos
  pos.diff <- c(as.integer(0), diff(pos))
  idx.jump <- which(pos.diff>sep|pos.diff<0)
  regions <- rbind(c(1, idx.jump), c(idx.jump-1, length(pos)))
  regions_bed=data.frame(chr = sig_df[regions[1,],"chr"],
                         start = sig_df[regions[1,],"pos"],
                         stop = sig_df[regions[2,],"pos"],
                         n_CG_relaxed = (regions[2,]-regions[1,])+1)
  regions_bed=regions_bed[regions_bed$n_CG_relaxed>=minCG,]
  for (i in 1:nrow(regions_bed)){
    chr=regions_bed[i,"chr"]
    start=regions_bed[i,"start"]
    stop=regions_bed[i,"stop"]
    tmp_df=input_df[input_df$chr==chr & input_df$pos>=start & input_df$pos<=stop,]
    n_CG=nrow(tmp_df)
    n_CG_strict=sum(as.numeric(tmp_df[,fdr_col_name])<fdr_cutoff_strict)
    median_effect=median(tmp_df[,effect_size_name])
    n_concordant=sum(sign(tmp_df[,effect_size_name])==sign(median_effect))
    n_sig_concordant=sum(sign(tmp_df[tmp_df[,fdr_col_name]<fdr_cutoff_relaxed,effect_size_name])==sign(median_effect))
    regions_bed[i,c("n_CG","n_CG_strict","median_effect","n_concordant","n_sig_concordant")]=c(n_CG,n_CG_strict,median_effect,n_concordant,n_sig_concordant)
  }
  regions_bed$p_sig=regions_bed$n_CG_relaxed/regions_bed$n_CG
  regions_bed$p_sig_concordant=regions_bed$n_sig_concordant/regions_bed$n_CG_relaxed
  regions_bed=regions_bed[regions_bed$p_sig>=prop_sig&regions_bed$n_CG_strict>0,]
  regions_bed$length=regions_bed$stop-regions_bed$start+1
  regions_bed
}


### Permutations to establish DMR CpG cutoff --------
# Use minCG=2 in  permutations: focal site and X number of other 10% sites.
perms_out <- do.call(rbind, lapply(X = 1:10, function(x) { 
  set.seed((1+1))
  dnam_dss_fmt_perm <- dnam_dss_fmt
  dnam_dss_fmt_perm$pval <- sample(x = dnam_dss_fmt_perm$pval, size = length(dnam_dss_fmt$pval), replace = F)
  dmrs2k_perm <- dmr_function(input_df=dnam_dss_fmt_perm, 
                              fdr_cutoff_strict=0.05, 
                              fdr_cutoff_relaxed=0.1, 
                              fdr_col_name="pval", 
                              effect_size_name="beta", 
                              minCG=2,
                              sep=1000, prop_sig=0.01)
  out <- data.frame(median_n_CG = median(dmrs2k_perm$n_CG), 
                    median_n_CG_relaxed = median(dmrs2k_perm$n_CG_relaxed), 
                    median_n_CG_strict = median(dmrs2k_perm$n_CG_strict), 
                    perm = x)
  return(out)
}))
perms_out$median_n_CG_relaxed
# --> 3 so we define DMR as counting at least 4 CpGs.


### Find DMRs in the real data ---------
fecal_both_dmr <- dmr_function(dnam_dss_fmt, fdr_cutoff_strict = 0.05, fdr_cutoff_relaxed = 0.1,
             fdr_col_name = "pval", effect_size_name = "beta", minCG = 4, sep = 1000,
             prop_sig = 0.01)

median(fecal_both_dmr$n_CG_relaxed)
median(fecal_both_dmr$n_CG_strict)
fecal_both_dmr$p_concordant <- (fecal_both_dmr$n_concordant/fecal_both_dmr$n_CG)

# Filter for 75% in same direction 
hist(fecal_both_dmr$p_sig_concordant,breaks = 50)
hist(fecal_both_dmr$p_concordant,breaks = 50)
fecal_both_dmr_filt <- fecal_both_dmr[fecal_both_dmr$p_concordant>=0.75 & fecal_both_dmr$p_sig_concordant>=0.75,]

# Check length of DMRs 
hist(fecal_both_dmr_filt$length,breaks=100)
summary(fecal_both_dmr_filt$length)
sd(fecal_both_dmr_filt$length)
quantile(fecal_both_dmr_filt$length, 0.99)

# How many extreme DMRs?
sum(fecal_both_dmr_filt$length>quantile(fecal_both_dmr_filt$length, 0.99))

fecal_both_dmr_filt <- fecal_both_dmr_filt %>% filter(length<quantile(fecal_both_dmr_filt$length, 0.99))

summary(fecal_both_dmr_filt$length)
sd(fecal_both_dmr_filt$length)
summary(fecal_both_dmr_filt$n_CG)
sd(fecal_both_dmr_filt$n_CG)
summary(fecal_both_dmr_filt)
sd(fecal_both_dmr_filt$n_CG_relaxed)
sd(fecal_both_dmr_filt$n_CG_strict)

# write an identifier region column
fecal_both_dmr_filt$region <- paste0(fecal_both_dmr_filt$chr, "_", fecal_both_dmr_filt$start, "_", fecal_both_dmr_filt$stop)


### Call background DMRs ----------
dmrs_bkgd_tmp <- dmr_function(dnam_dss_fmt, fdr_cutoff_strict = 1, fdr_cutoff_relaxed = 1,
                              fdr_col_name = "pval", effect_size_name = "beta", minCG = 4, sep = 1000,
                              prop_sig = 0.01)

# Keep minCG at the same level, and set pvalue to 1 (akin to permuting) 
dmrs_bkgd_tmp$p_concordant <- (dmrs_bkgd_tmp$n_concordant/dmrs_bkgd_tmp$n_CG)

# Create GRanges objects
fecal_both_dmr_filt_gr <- GRanges(
  seqnames = fecal_both_dmr_filt$chr,
  ranges = IRanges(start = fecal_both_dmr_filt$start, end = fecal_both_dmr_filt$stop)
)

dmrs_bkgd_tmp_gr <- GRanges(
  seqnames = dmrs_bkgd_tmp$chr,
  ranges = IRanges(start = dmrs_bkgd_tmp$start, end = dmrs_bkgd_tmp$stop)
)

# Find overlaps
overlaps <- findOverlaps(dmrs_bkgd_tmp_gr, fecal_both_dmr_filt_gr)

# Get the indices of the background regions that overlap with the real DMRs
overlap_indices <- queryHits(overlaps)

# Remove the overlapping regions from the background set
dmrs_bkgd_tmp_gr_filt <- dmrs_bkgd_tmp_gr[-overlap_indices]

dmrs_bkgd_tmp_filt <- as.data.frame(dmrs_bkgd_tmp_gr_filt)

# Check length of background 
hist(dmrs_bkgd_tmp_filt$width,breaks=100)
summary(dmrs_bkgd_tmp_filt$width)
sd(dmrs_bkgd_tmp_filt$width)

## remove DMRs which are much longer than the significantly associated DMRs
sum(dmrs_bkgd_tmp_filt$width>max(fecal_both_dmr_filt$length))
dmrs_bkgd_tmp_filt <- dmrs_bkgd_tmp_filt %>% filter(width<=max(fecal_both_dmr_filt$length))
summary(dmrs_bkgd_tmp_filt$width)
sd(dmrs_bkgd_tmp_filt$width)


## Before saving check for overlaps (within contig)
overlaps <- fecal_both_dmr_filt %>%
  group_by(chr) %>%
  arrange(start) %>% # Order by start within each group
  mutate(next_start = lead(start),
         next_stop = lead(stop)) %>%
  filter(stop > next_start) %>% # Check for overlaps
  dplyr::select(chr, start, stop, next_start, next_stop)

# Print overlapping regions if any
if(nrow(overlaps) > 0) {
  print("Overlapping regions:")
  print(overlaps)
} else {
  print("No overlapping regions found.")
}

overlaps <- dmrs_bkgd_tmp_filt %>%
  group_by(seqnames) %>%
  arrange(start) %>% # Order by start within each group
  mutate(next_start = lead(start),
         next_end = lead(end)) %>%
  filter(end > next_start) %>% # Check for overlaps
  dplyr::select(seqnames, start, end, next_start, next_end)

# Print overlapping regions if any
if(nrow(overlaps) > 0) {
  print("Overlapping regions:")
  print(overlaps)
} else {
  print("No overlapping regions found.")
}

write.table(dmrs_bkgd_tmp_filt, "/path/to/aDMRs_background.txt", sep = "\t", row.names = F,quote = F)
write.table(fecal_both_dmr_filt, "/path/to/aDMRs_grouping.txt", sep = "\t", row.names = F,quote = F)

####### END OF THE CODE ############
