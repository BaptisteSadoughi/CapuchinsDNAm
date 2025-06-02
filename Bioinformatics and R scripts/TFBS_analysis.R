###### TRANSCRIPTION FACTOR BINDING SITES ENRICHMENT ANALYSIS

## For an introduction to monaLisa see https://fmicompbio.github.io/monaLisa/articles/monaLisa.html#quick

rm(list=ls())

library_list <- c("tidyverse","stringr","GenomicRanges","SummarizedExperiment","TFBSTools","JASPAR2020","monaLisa")
lapply(library_list, require, character.only=TRUE)

### Download a database of vetebrate known motifs

pwms <- getMatrixSet(JASPAR2020,
                     opts = list(matrixtype = "PWM",
                                 tax_group = "vertebrates"))

### Load background and differentially methylated regions

dmrs_bkgd <- read.table("/path/to/aDMRs_background.txt", sep = "\t", header= TRUE)
dmrfile <- read.table("/path/to/aDMRs_grouping.txt", sep = "\t", header = TRUE)

# Turn to a GRanges
dmr <- GRanges(
  seqnames = dmrfile$chr,
  ranges = IRanges(start = dmrfile$start, end = dmrfile$stop),
  deltaMeth=dmrfile$median_effect
)

bkgd <- GRanges(
  seqnames = dmrs_bkgd$seqnames,
  ranges = IRanges(start = dmrs_bkgd$start, end = dmrs_bkgd$end)
)

# Create sets of upmethylated, downmethylated, and join with background sets
dmr_pos <- dmr[dmr$deltaMeth > 0]
length(dmr_pos)

dmr_neg <- dmr[dmr$deltaMeth < 0]
length(dmr_neg)

# combine the two sets or genomic regions
dmr_pos_test <- c(bkgd, dmr_pos)

dmr_neg_test <- c(bkgd, dmr_neg)

# Define binning vector for the regions
bins_pos <- rep(c("unchanged", "up"), c(length(bkgd), length(dmr_pos)))
bins_pos <- factor(bins_pos)
table(bins_pos)

bins_neg <- rep(c("unchanged", "down"), c(length(bkgd), length(dmr_neg)))
bins_neg <- factor(bins_neg)
table(bins_neg)

### Extract sequences corresponding to the regions

# Specify your FASTA file
fasta_file <- "/path/to/Cebus_imitator.Cebus_imitator-1.0.dna.toplevel.fa"

# Create a FaFile object
fa_file <- FaFile(fasta_file)

## Index the FASTA file (only needs to be done once)
# indexFa(fasta_file)

# Use getSeq to extract the sequences
dmr_pos_test_seqs <- Biostrings::getSeq(fa_file, dmr_pos_test)

dmr_neg_test_seqs <- Biostrings::getSeq(fa_file, dmr_neg_test)

# The resulting seqs only use contigs as names, which creates duplicated names and leads to an ERROR in monaLisa.

# Create unique names by combining regions' seqnames, start, and end
unique_names_pos <- paste0(
  seqnames(dmr_pos_test), "_", 
  start(dmr_pos_test), "-", 
  end(dmr_pos_test)
)

# Assign these unique names to the DNAStringSet
names(dmr_pos_test_seqs) <- unique_names_pos

unique_names_neg <- paste0(
  seqnames(dmr_neg_test), "_", 
  start(dmr_neg_test), "-", 
  end(dmr_neg_test)
)
names(dmr_neg_test_seqs) <- unique_names_neg

### Run TBFS enrichement

##### DMR gaining methylation with age

# use background
ME_pos <- calcBinnedMotifEnrR(seqs = dmr_pos_test_seqs, bins = bins_pos,
                           pwmL = pwms)


## Note we get a warning because they are "N" bases 
{
# # Find sequences with non-ACGT characters
# problematic_seqs <- dmr_pos_test_seqs[grep("[^ACGT]", as.character(dmr_pos_test_seqs))]
# if (length(problematic_seqs) > 0) {
#   print(problematic_seqs)
# } else {
#   print("All sequences contain valid ACGT characters.")
# }
}

dim(ME_pos) # motifs-by-bins
rowData(ME_pos)
colData(ME_pos)
assayNames(ME_pos)
assay(ME_pos, "log2enr")[1:3, 1:2]

# select strongly enriched motifs
sel_pos <- apply(assay(ME_pos, "negLog10Padj"), 1, 
             function(x) max(abs(x), 0, na.rm = TRUE)) > 3.0
sum(sel_pos)

# Calculate motif similarity based on pearson 
SimMatSel <- motifSimilarity(rowData(ME_pos[sel_pos,])$motif.pfm)
range(SimMatSel)

# create hclust object, similarity defined by 1 - Pearson correlation
hcl <- hclust(as.dist(1 - SimMatSel), method = "average")

# reminder the heatmap display the columns in the order of
levels(bins_pos) #unchanged, up

# Add dummy second column to force ncol == nlevels(bin)
ME_workaround <- ME_pos[, c("unchanged", "up"), drop = FALSE]

# Replace all values in the dummy "unchanged" column with NA
for (assay_name in assayNames(ME_workaround)) {
  assay(ME_workaround, assay_name)[, "unchanged"] <- NA
}

# Set bin column explicitly
colData(ME_workaround)$bin <- factor(c("unchanged", "up"), levels = c("unchanged", "up"))

plotMotifHeatmaps(
  x = ME_workaround[sel_pos,],
  which.plots = c("log2enr", "negLog10Padj"), 
  width = 1.8, maxEnr = 2, maxSig = 10,
  cluster = hcl,
  show_dendrogram = TRUE,
  show_seqlogo = TRUE,
  show_motif_GC = FALSE,
  width.seqlogo = 1.2
)
dev.off()

# ##### DMR losing methylation with age

# use background
ME_neg <- calcBinnedMotifEnrR(seqs = dmr_neg_test_seqs, bins = bins_neg,
                              pwmL = pwms)

dim(ME_neg) # motifs-by-bins
rowData(ME_neg)
colData(ME_neg)
assayNames(ME_neg)
assay(ME_neg, "log2enr")[1:3, 1:2]
assay(ME_neg, "negLog10Padj")[1:100, 1:2]

# select strongly enriched motifs
sel_neg <- apply(assay(ME_neg, "negLog10Padj"), 1, 
                 function(x) max(abs(x), 0, na.rm = TRUE)) > 2.0
sum(sel_neg)

########### END OF THE CODE ##########
