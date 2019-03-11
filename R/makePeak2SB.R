
.makePeak2SB <- function(peakGL, psuedoGene, extend = 200) {

  peak_name <- mcols(peakGL)$V4
  peak_name[is.na(peak_name)] <- "Unk"

  ps_name <- names(psuedoGene)
  names(psuedoGene) <- paste("psGene", 1:length(psuedoGene), sep = "")

  ind_tr <- findOverlaps(peakGL, psuedoGene, type = "within")
  qname <- peak_name[queryHits(ind_tr)]
  sname <- ps_name[subjectHits(ind_tr)]
  ind_tr <- ind_tr[qname == sname]

  ps_width <- sum(width(psuedoGene))
  ind_width <- ps_width[subjectHits(ind_tr)]
  ind_map <- tapply(ind_width, queryHits(ind_tr), function(x){names(x)[which.max(x)]})

  peakStart <- lapply(start(peakGL), min)
  peakEnd <- lapply(end(peakGL), max)
  peakChr <- lapply(seqnames(peakGL), unique)
  peakChr <- as.character(unlist(peakChr))
  peakStrand <- lapply(strand(peakGL), unique)
  peakStrand <- as.character(unlist(peakStrand))

  peakSG <- GRanges(seqnames = peakChr,
                    IRanges(start = unlist(peakStart), width = 1),
                    strand = peakStrand)
  names(peakSG) <- names(peakStart)

  peakEG <- GRanges(seqnames = peakChr,
                    IRanges(start = unlist(peakEnd), width = 1),
                    strand = peakStrand)
  names(peakEG) <- names(peakEnd)

  PeakTxStart <- pmapToTranscripts(peakSG, psuedoGene[as.character(ind_map)], ignore.strand = T)
  PeakTxEnd <- pmapToTranscripts(peakEG, psuedoGene[as.character(ind_map)], ignore.strand = T)

  PeakTxStart <- start(PeakTxStart)
  PeakTxEnd <- start(PeakTxEnd)
  PeakName <- names(peakGL)
  PeakTxName <- as.character(ind_map)

  len = length(peakGL)
  PeakM <- rbind(PeakTxName, PeakTxStart, PeakTxEnd, peakStrand)
  PeakL <- as.vector(PeakM)
  ind <- rep(1:len, each = 4)
  PeakL <- split(PeakL, ind)
  names(PeakL) <- PeakName

  PeakSB <- lapply(PeakL, .makeSBpeak, extend)

  PeakSB <- GRangesList(c(unlist(PeakSB)))

  return(PeakSB)
}

.makeSBpeak <- function(x, extend) {
  gr0 <- GRanges(x[1],
                 IRanges((as.numeric(x[2]) - extend):(as.numeric(x[3]) + extend), width = 1),
                 strand = x[4])
}

