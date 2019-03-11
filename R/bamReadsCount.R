
.bamReadsCount <- function(meripsp, file, tx_genes,
                           minimal_alignment_MAPQ = 30,
                           fragment_length = 100,
                           bin_width = 51) {

  binOfPsuedoGene <- unlist(meripsp)
  psuedoGene <- tx_genes
  names(psuedoGene) <- paste("psGene", 1:length(psuedoGene), sep = "")
  # prepare bam parameters
  what <- c("rname","strand", "pos","mapq","qwidth")
  param <- ScanBamParam(what=what)

  # read bam file
  ba <- scanBam(file, param=param)
  ba <- ba[[1]]
  total_reads_count <- length(ba$pos)

  # MAPQ filter
  ba$mapq[which(is.na(ba$mapq))] <- 255
  ba$rname <- ba$rname[ba$mapq > minimal_alignment_MAPQ]
  ba$strand <- ba$strand[ba$mapq > minimal_alignment_MAPQ]
  ba$pos <- ba$pos[ba$mapq > minimal_alignment_MAPQ]
  ba$qwidth <- ba$qwidth[ba$mapq > minimal_alignment_MAPQ]
  ba$mapq <- ba$mapq[ba$mapq > minimal_alignment_MAPQ]
  total_reads_count_filtered <- length(ba$pos)

  gc()

  ID_negative <- (ba$strand=="-")
  ba$pos[ID_negative] <- ba$pos[ID_negative] + ba$qwidth[ID_negative] - 1

  rm(ID_negative)
  gc()

  gr <- GRanges(seqnames = ba$rname,
                ranges = IRanges(start=ba$pos, width = 1),
                strand = ba$strand)

  rm(ba)
  gc()


  # gene reads count
  gene_reads_count <- countOverlaps(psuedoGene, gr, ignore.strand = T)

  # unimaped reads
  gr <- gr[countOverlaps(gr, psuedoGene, ignore.strand = T) > 0]
  total_reads_count_unimapped <- length(gr)

  # shift
  grt <- mapToTranscripts(gr, psuedoGene, ignore.strand = T)
  strand(grt) <- strand(gr)[mcols(grt)$xHits]
  grt <- resize(grt, width = fragment_length, fix="start", ignore.strand = F)

  # bin reads count
  gr0 <- resize(binOfPsuedoGene, width = bin_width, fix = "center")
  base_reads_count <- as.numeric(countOverlaps(binOfPsuedoGene, grt, ignore.strand = T))
  bin_reads_count <- as.numeric(countOverlaps(gr0, grt, ignore.strand = T))

  rm(gr, grt)
  gc()

  # result
  # total_reads_count <- c(total_reads_count, total_reads_count_filtered, total_reads_count_unimapped)
  # names(total_reads_count) <- c("total", "filtered", "unimapped")
  # total_reads_count <- data.frame(total_reads_count)

  result <- list(
    total_reads_count = total_reads_count,
    gene_reads_count = gene_reads_count,
    base_reads_count = base_reads_count,
    bin_reads_count = bin_reads_count
  )
  return(result)
  gc()
}
