#### Two support functions used to intersect several bbseq objects and retain only those sites common to all of them. 

# FUNCTION to convert bsseq list to GRangesList
bsseq_to_gr <- function(bsseq) {
  chrom <- seqnames(bsseq)
  start <- start(bsseq)
  end <- end(bsseq)
  gr <- GRanges(seqnames = chrom, ranges = IRanges(start = start, end = end))
  return(gr)
}


# FUNCTION to intersect and filter common sites among three bsseq objects:
extract_common_cpgs <- function(bsseq_1, bsseq_2,
                                     bsseq_common_1_name = "bsseq_1_common", 
                                     bsseq_common_2_name = "bsseq_2_common") {
  
  # Create GRanges objects from bsseq_1, bsseq_2, and bsseq_3
  bsseq_gr_1 <- bsseq_to_gr(bsseq_1)
  bsseq_gr_2 <- bsseq_to_gr(bsseq_2)
  
  # Find common CpG sites among bsseq_gr_1, bsseq_gr_2, and bsseq_gr_3
  common_sites <- GenomicRanges::intersect(bsseq_gr_1, bsseq_gr_2)
  
  # Extract bsseq objects for the common sites
  bsseq_common_1 <- subsetByOverlaps(bsseq_1, common_sites)
  bsseq_common_2 <- subsetByOverlaps(bsseq_2, common_sites)
  
  # Return the three resulting bsseq objects as a named list
  result <- list(bsseq_common_1 = bsseq_common_1, bsseq_common_2 = bsseq_common_2)
  names(result) <- c(bsseq_common_1_name, bsseq_common_2_name)
  
  return(result)
}
