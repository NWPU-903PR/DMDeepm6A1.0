

###################################################################################################################
## extract component from txdb
.extractComponent <- function(txdb,
                              maximalAmbiguity = NA,
                              minimalComponentLength = 100,
                              minimalNcRNALength = 200){

  parameter = list()
  parameter$txdb <- txdb
  parameter$maximalAmbiguity <- maximalAmbiguity # whether overlap with another transcript
  parameter$minimalComponentLength <- minimalComponentLength # minimal length required for each component
  parameter$minimalNcRNALength <- minimalNcRNALength

  txdb <- parameter$txdb
  # ambiguity filter
  exons <- exonsBy(txdb, by = "tx",use.names=TRUE)
  noTx <- length(exons)
  print(paste("total",noTx,"transcripts extracted ..."));

  if (!is.na(maximalAmbiguity)) {
    temp <- countOverlaps(exons, exons)
    ambiguityFilteredTx <- names(exons[temp < (parameter$maximalAmbiguity+2)])
    noTxLeft <- length(ambiguityFilteredTx)
    print(paste("total",noTxLeft,"transcripts left after ambiguity filter ..."))
    exons <- exons[ambiguityFilteredTx]
  }

  # extract important components
  cds <- cdsBy(txdb, by = "tx",use.names=TRUE)
  utr5 <- fiveUTRsByTranscript(txdb, use.names=TRUE)
  utr3 <- threeUTRsByTranscript(txdb, use.names=TRUE)

  # extract mRNAs
  flag_utr5 <- (sum(width(utr5)) > parameter$minimalComponentLength)
  name_utr5 <- names(utr5)[flag_utr5]
  flag_utr3 <- (sum(width(utr3)) > parameter$minimalComponentLength)
  name_utr3 <- names(utr3)[flag_utr3]
  flag_cds <- (sum(width(cds)) > parameter$minimalComponentLength)
  name_cds <- names(cds)[flag_cds]
  name_mRNA <- unique(c(name_utr5, name_utr3, name_cds))
  name_filtered_mRNA <- intersect(name_mRNA,names(exons))
  cds_filtered <- cds[intersect(name_filtered_mRNA, name_cds)]
  utr5_filtered <- utr5[intersect(name_filtered_mRNA, name_utr5)]
  utr3_filtered <- utr3[intersect(name_filtered_mRNA, name_utr3)]
  print(paste("total",length(cds_filtered),"mRNAs left after component length filter ..."))

  # extract mRNAs
  all_mRNA <- unique(c(names(utr5),names(utr3),names(cds)))
  name_ncRNA <- setdiff(names(exons),all_mRNA)
  ncRNA <- exons[name_ncRNA]
  flag_ncRNA <-
    (sum(width(ncRNA)) > parameter$minimalComponentLength) &
    (sum(width(ncRNA)) > parameter$minimalNcRNALength)
  name_ncRNA <- names(ncRNA)[flag_ncRNA]
  ncRNA_filtered <- ncRNA[name_ncRNA]
  print(paste("total",length(ncRNA_filtered),"ncRNAs left after ncRNA length filter ..."))

  # return the result
  comp <- list(cds=cds_filtered,utr3=utr3_filtered,utr5=utr5_filtered,ncRNA=ncRNA_filtered)
  return(comp)
}

## get peak position
.getpeakposition <- function(peak, comp, txdb, egSYMBOL, outfilepath){

  peak_gr <- GRanges(as.character(peak$chr),
                     IRanges(as.numeric(peak$chromEnd), width = 1),
                     as.character(peak$strand))


  tx <- exonsBy(txdb, by = "tx")
  id <- findOverlaps(peak_gr, tx)
  tx <- tx[unique(subjectHits(id))]

  y <- mapToTranscripts(peak_gr, tx)
  ind <- tapply(start(y), mcols(y)$xHits, min)

  st <- 1:nrow(peak)
  st <- as.numeric(ind[match(st, as.numeric(names(ind)))])

  TxStart <- st
  xls <- peak
  xls <- cbind(xls, TxStart)

  utr3_peak <- countOverlaps(peak_gr, comp$utr3)
  utr5_peak <- countOverlaps(peak_gr, comp$utr5)
  cds_peak <- countOverlaps(peak_gr, comp$cds)
  ncrna_peak <- countOverlaps(peak_gr, comp$ncRNA)
  position <- data.frame(UTR5 = utr5_peak, CDS = cds_peak, UTR3 = utr3_peak, LncRNA = ncrna_peak)
  position[position != 0] <- 1
  xls <- cbind(xls, position)

  x = egSYMBOL
  mapped_genes <- mappedkeys(x)
  result <- as.list(x[mapped_genes])
  entrez_id <- as.numeric(names(result))  # entrez ID
  gene_symbol <- as.character(result)     # gene symbol
  ID_convert <- data.frame(entrez_id,gene_symbol)
  id <- xls$name
  gene_symbol <- as.vector(ID_convert$gene_symbol[match(id,ID_convert$entrez_id)])

  if (sum(is.na(gene_symbol)) == length(gene_symbol)) {gene_symbol <- id}

  xls <- cbind(xls, gene_symbol)
  names(xls)[ncol(xls)] <- "GeneSymbol"
  write.table(xls, file =  paste(outfilepath, "CandidateSingleBasePeak.xls" , sep = "/"),
              sep = "\t", row.names = FALSE, quote = FALSE)
  return(xls)
}















