

.makePeak2Motif <- function(peakGL, parameter, extend = 50) {

  ## get peak mapped transcript
  peak_name <- mcols(peakGL)$V4
  peak_name[is.na(peak_name)] <- "Unk"

  psuedoGene <- parameter$tx_genes
  BSgenome <- parameter$BSgenome
  motif <- parameter$motif
  ps_name <- names(psuedoGene)
  names(psuedoGene) <- paste("psGene", 1:length(psuedoGene), sep = "")

  ind_tr <- findOverlaps(peakGL, psuedoGene, type = "within")
  qname <- peak_name[queryHits(ind_tr)]
  sname <- ps_name[subjectHits(ind_tr)]
  ind_tr <- ind_tr[qname == sname]

  ps_width <- sum(width(psuedoGene))
  ind_width <- ps_width[subjectHits(ind_tr)]
  ind_map <- tapply(ind_width, queryHits(ind_tr), function(x){names(x)[which.max(x)]})

  ## get peak motif location
  genome <- BSgenome

  ind_strand <- unlist(unique(strand(peakGL))) == "-"
  peakGL[ind_strand] <- sort(peakGL[ind_strand], decreasing = T)

  peak_seq <- getSeq(genome, peakGL)
  peak_seq_con <- lapply(peak_seq, unlist)
  peak_seq_con <- DNAStringSet(peak_seq_con)

  cag_loc <- vmatchPattern(motif[1], peak_seq_con)
  motif_start <- start(cag_loc)

  for(i in 2:length(motif)) {
    cag_loc0 <- vmatchPattern(motif[i], peak_seq_con)
    motif_start0 <- start(cag_loc0)
    motif_start <- mapply(c, motif_start, motif_start0, SIMPLIFY=FALSE)
  }

  ## get motif centered sequance
  tx <- psuedoGene[as.character(ind_map)]

  ind_strand <- unlist(unique(strand(tx))) == "-"
  tx[ind_strand] <- sort(tx[ind_strand], decreasing = T)

  motif_ind <- rep(1:length(peakGL), times = unlist(lapply(motif_start, length)))
  motif_loci <- unlist(motif_start) + 2
  motif_strand <- as.character(unique(strand(peakGL)))[motif_ind]
  motif_name <- names(peakGL)[motif_ind]

  motif_tr <- GRanges(motif_name,
                      IRanges(start = motif_loci, width = 1),
                      motif_strand)

  motif_gr <- mapFromTranscripts(motif_tr, peakGL)
  motif_tx <- pmapToTranscripts(motif_gr, tx[motif_ind])

  len = length(motif_tx)
  motifM <- rbind(motif_ind, start(motif_tx))
  motifL <- as.vector(motifM)
  ind <- rep(1:len, each = 2)
  motifL <- split(motifL, ind)

  tx_seq <- extractTranscriptSeqs(genome, tx)
  seq_extend <- DNAString(paste(rep("N", extend), collapse = ""))
  tx_seq <- lapply(tx_seq, .extendSeq, seq_extend)

  motif_seq <-  lapply(motifL, .getMotifSeq, tx_seq, extend)

  motifBin <- list(motif_A_loci = motif_gr,
                   motif_seq = motif_seq,
                   motif_ind = motif_ind)

  return(motifBin)
}

.getMotifSeq <- function(x, tx_seq, extend) {
  gr_seq <- tx_seq[[x[1]]]
  return(gr_seq[x[2]:(x[2] + 2*extend)])
}

.extendSeq <- function(x, seq_extend) {
  y <- c(seq_extend, x, seq_extend)
}



