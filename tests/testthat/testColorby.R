library(GenomicRanges)
devtools::load_all()
# Create two simple GRanges objects
gr1 <- GRanges(seqnames = "chr1",
               ranges = IRanges(start = c(1, 10, 20),
                                end = c(5, 15, 25)),
               strand = c("+", "+", "+"))

gr2 <- GRanges(seqnames = "chr1",
               ranges = IRanges(start = c(100, 200),
                                end = c(150, 250)),
               strand = c("-", "-"))
gr3 <- GRanges(seqnames = "chr1",
               ranges = IRanges(start = c(100, 200),
                                end = c(150, 250)),
               strand = c("-", "-"))
# Now combine them into a GRangesList
grl <- GRangesList(gr1 = gr1, gr2 = gr2, gr3 = gr3)

# View the GRangesList
grl

#' transcript_annotations Data frame with at least three columns: transcript_id, gene_name, strand.
#' Used to construct transcript labels. (default: NULL)

transcript_annotations <- dplyr::tibble(
    transcript_id = names(grl),
    strand = sapply(grl, function(gr) {
      strands <- as.character(strand(gr))
      unique_strands <- unique(strands)
        
        if (length(unique_strands) == 1) {
          unique_strands
        } else {
          "*"  # if multiple strands found (inconsistent), mark "*"
        }
      }),
    color_by = c('lightgreen','#fdae61','red')
)

wiggleplotr::plotTranscripts(grl, transcript_annotations = transcript_annotations)
