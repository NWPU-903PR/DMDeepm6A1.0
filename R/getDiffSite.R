

## get thresh conind
.conind <- function(bedpath, CONDITION, sig_thresh) {

  treated <- bedpath[CONDITION == "treated"]
  untreated <- bedpath[CONDITION == "untreated"]

  len_tr <- length(treated)
  grl_tr <- GRangesList()
  for (i in 1:len_tr) {
    grl_tr0 <- import.bed(paste(treated[i], "CandidateSingleBasePeak.bed", sep = "/"))
    grl_tr0 <- grl_tr0[mcols(grl_tr0)$score >= sig_thresh]
    grl_tr[[i]] <- grl_tr0
  }
  names(grl_tr) <- paste("treated", 1:len_tr)

  len_un <- length(untreated)
  grl_un <- GRangesList()
  for (i in 1:len_un) {
    grl_un0 <- import.bed(paste(untreated[i], "CandidateSingleBasePeak.bed", sep = "/"))
    grl_un0 <- grl_un0[mcols(grl_un0)$score >= sig_thresh]
    grl_un[[i]] <- grl_un0
  }
  names(grl_un) <- paste("untreated", 1:len_un)

  gr_unique <- unique(unlist(c(grl_tr, grl_un)))

  con_ind_tr <- countOverlaps(gr_unique, grl_tr)
  con_ind_un <- countOverlaps(gr_unique, grl_un)
  con_ind <- (con_ind_tr > 0) & (con_ind_un > 0)

  hyper_ind <- (con_ind_tr == len_tr) & (con_ind_un == 0)
  hypo_ind <- (con_ind_tr == 0) & (con_ind_un == len_un)
  candi_ind <- (!hyper_ind) & (!hypo_ind) & (con_ind)

  xls <- data.frame(hyper_ind = hyper_ind, hypo_ind = hypo_ind, candi_ind = candi_ind)
  mcols(gr_unique) <- xls

  return(list(grl_tr = grl_tr, grl_un = grl_un, gr_unique = gr_unique))
}


## test con site
.diffmethytest <- function(parameter) {

  diffsite_ind <- parameter$Diff$diffsite_ind
  tx_genes <- parameter$tx_genes
  bams <- parameter$Diff$bams
  bam_condition <- parameter$Diff$bam_condition
  sf <- parameter$Diff$sf
  output_filepath <- parameter$output_filepath
  diff_bin_width <- parameter$diff_bin_width

  gr <- diffsite_ind$gr_unique
  gr <- gr[mcols(gr)$candi_ind]

  psuedoGene <- tx_genes
  names(psuedoGene) <- paste("psGene", 1:length(psuedoGene), sep = "")
  grt <- mapToTranscripts(gr, psuedoGene, ignore.strand = T)
  psuedoGene <- psuedoGene[unique(mcols(grt)$transcriptsHits)]

  tx_len <- width(psuedoGene)
  tx_len <- sum(tx_len)

  tx_ind <- tapply(tx_len[mcols(grt)$transcriptsHits], mcols(grt)$xHits, max)
  tx_ind <- paste(names(tx_ind), tx_ind)
  gr_ind <- paste(mcols(grt)$xHits, tx_len[mcols(grt)$transcriptsHits])

  grt <- mapToTranscripts(gr, psuedoGene, ignore.strand = T)
  grt_u <- grt[match(tx_ind, gr_ind)]

  bamreads <- matrix(0, length(grt_u), 1)
  sf0 <- numeric()
  len <- length(bams)
  for (i in 1:len) {
    print(paste("Count reads for candidate DM peak in bam", i))
    bamreads0 <- .bamReadsCount(GRangesList(grt_u), as.character(bams[i]), tx_genes, bin_width = diff_bin_width)
    bamreads <- cbind(bamreads, bamreads0$bin_reads_count)
    sf0[i] <- as.numeric(bamreads0$total_reads_count)
  }

  bamreads <- bamreads[,-1]

  if (parameter$diff_method == "QNB") {

    if (!is.na(sf[1])) {

      if (sum(sf) > 1000) {
        standard_library_size <- exp(mean(log(sf)))
        sf <- sf/standard_library_size
      }

      names(sf) <- bam_condition

      sfl <- list(control_ip = as.numeric(sf[names(sf) == "untreated_ip"]),
                  treated_ip = as.numeric(sf[names(sf) == "treated_ip"]),
                  control_input = as.numeric(sf[names(sf) == "untreated_input"]),
                  treated_input = as.numeric(sf[names(sf) == "treated_input"]))
    } else {
      sfl <- NA
    }

    control_ip <- bamreads[,bam_condition == "untreated_ip"]
    treated_ip <- bamreads[,bam_condition == "treated_ip"]
    control_input <- bamreads[,bam_condition == "untreated_input"]
    treated_input <- bamreads[,bam_condition == "treated_input"]

    dir.create(paste(output_filepath, "QNB_out", sep = "/"))
    result <- qnbtest(control_ip, treated_ip, control_input, treated_input,
                      size.factor = sfl, output.dir = paste(output_filepath, "QNB_out", sep = "/"))
  }

  if (parameter$diff_method == "exomepeak") {

    untreated_ip <- bamreads[,bam_condition == "untreated_ip"]
    treated_ip <- bamreads[,bam_condition == "treated_ip"]
    untreated_input <- bamreads[,bam_condition == "untreated_input"]
    treated_input <- bamreads[,bam_condition == "treated_input"]

    # get reads count
    if (ncol(untreated_ip)>1) {untreated_ip=rowSums(untreated_ip)}
    if (ncol(untreated_input)>1) {untreated_input=rowSums(untreated_input)}
    if (ncol(treated_ip)>1) {treated_ip=rowSums(treated_ip)}
    if (ncol(treated_input)>1) {treated_input=rowSums(treated_input)}


    if (!is.na(sf[1])) {

      sf <- round(sf*30/200)
      names(sf) <- bam_condition

    } else {
      sf <- sf0
      names(sf) <- bam_condition
    }

    # get total
    untreated_ip_total=sum(as.numeric(sf[names(sf) == "untreated_ip"]))
    untreated_input_total=sum(as.numeric(sf[names(sf) == "untreated_input"]))
    treated_ip_total=sum(as.numeric(sf[names(sf) == "treated_ip"]))
    treated_input_total=sum(as.numeric(sf[names(sf) == "treated_input"]))

    result <- rhtest(untreated_ip,untreated_input,treated_ip,treated_input,
                     untreated_ip_total,untreated_input_total,treated_ip_total,treated_input_total)

    result <- cbind(NA, NA, NA, log2(exp(result$log.fc)), exp(result$log.p), NA, exp(result$log.fdr))
  }

  return(result)
}


## annosite
.annosite <- function(parameter, sig_site_thresh) {

  diffsite_ind <- parameter$Diff$diffsite_ind
  diffsite_test_res <- parameter$Diff$diffsite_test_res
  bedpath <- parameter$Diff$bedpath
  CONDITION <- parameter$sample_conditions


  gr <- diffsite_ind$gr_unique
  gr_ind <- mcols(gr)
  gr_hyper <- gr[gr_ind$hyper_ind]
  gr_hypo <- gr[gr_ind$hypo_ind]
  gr_diff <- gr[gr_ind$candi_ind]

  treated <- bedpath[CONDITION == "treated"]
  untreated <- bedpath[CONDITION == "untreated"]

  xls_hyper <- .getxls(paste(treated, "CandidateSingleBasePeak.xls", sep = "/"), gr_hyper, sig_site_thresh)
  xls_hypo <- .getxls(paste(untreated, "CandidateSingleBasePeak.xls", sep = "/"), gr_hypo, sig_site_thresh)
  xls_diff <- .getxls(paste(bedpath, "CandidateSingleBasePeak.xls", sep = "/"), gr_diff, sig_site_thresh)

  anno_ind <- paste(xls_diff$chr, xls_diff$chromEnd)
  diff_ind <- paste(as.character(seqnames(gr_diff)), start(gr_diff))
  anno_ind <- match(anno_ind, diff_ind)
  diff_anno <- diffsite_test_res[anno_ind,]
  colnames(diff_anno) <- c("p.methy.tr", "p.methy.ctrl", "diff.log2.RR",
                        "diff.log2.OR", "diff.pvalue",  "diff.q", "diff.padj")

  xls_diff <- cbind(xls_diff[,1:11], diff_anno, xls_diff[,12:17])

  return(list(xls_hyper = xls_hyper, xls_hypo = xls_hypo, xls_diff = xls_diff))
}

.getxls <- function(xlspath, site_gr, sig_site_thresh) {

  len <- length(xlspath)
  xls <- data.frame()
  for (i in 1:len) {
    xls0 <- read.table(xlspath[i], header = T, stringsAsFactors = F)
    xls0 <- xls0[xls0$Probability >= sig_site_thresh,]
    # Replicate <- rep(paste(rep_label, i, sep = "_"), nrow(xls0))
    # xls0 <- cbind(xls0, Replicate)
    xls <- rbind(xls, xls0)
  }

  gr_ind <- paste(as.character(seqnames(site_gr)), start(site_gr))
  xls_ind <- paste(xls$chr, xls$chromEnd)
  xls <- xls[is.element(xls_ind, gr_ind),]

  xls <- xls[order(xls$fold_enrchment, decreasing = T),]

  site_ind <- paste(xls$chr, xls$chromEnd, xls$name)
  site_ind_unique <- sort(unique(site_ind))
  site_ind <- match(site_ind_unique, site_ind)
  xls <- xls[site_ind,]

  xls1 <- xls[!is.na(xls$name),]
  xls1 <- xls1[order(xls1$chr),]
  xls2 <- xls[is.na(xls$name),]
  xls2 <- xls1[order(xls2$chr),]

  xls <- rbind(xls1, xls2)

  return(xls)
}

.getconsensussite <- function(bedpath) {

  len <- length(bedpath)
  xlspath <- paste(bedpath, "SigSingleBasePeak.xls", sep = "/")
  bedpath <- paste(bedpath, "SigSingleBasePeak.bed", sep = "/")

  xls <- data.frame()
  grl <- GRangesList()
  for (i in 1:len) {
    xls0 <- read.table(xlspath[i], header = T, stringsAsFactors = F)
    xls <- rbind(xls, xls0)
    grl[[i]] <- import.bed(bedpath[i])
  }

  xls <- xls[order(xls$fold_enrchment, decreasing = T),]

  site_ind <- paste(xls$chr, xls$chromEnd, xls$name)
  site_ind_unique <- sort(unique(site_ind))
  site_ind <- match(site_ind_unique, site_ind)
  xls <- xls[site_ind,]

  xls1 <- xls[!is.na(xls$name),]
  xls1 <- xls1[order(xls1$chr),]
  xls2 <- xls[is.na(xls$name),]
  xls2 <- xls1[order(xls2$chr),]

  xls <- rbind(xls1, xls2)

  gru <- GRanges(as.character(xls$chr),
                 IRanges(start = as.numeric(xls$chromEnd), width = 1),
                 as.character(xls$strand))
  con_ind <- countOverlaps(gru, grl) == len

  xls_con <- xls[con_ind,]

  return(list(xls = xls, xls_con = xls_con))
}

## write site anno

.writeSiteAnno <- function(xls_anno, outfilepath, xls_name) {

  xls_anno$chromStart <- as.integer(xls_anno$chromStart)
  xls_anno$chromEnd <- as.integer(xls_anno$chromEnd)
  write.table(xls_anno, file =  paste(outfilepath, "/", xls_name, ".xls", sep = ""),
              sep = "\t", row.names = FALSE, quote = FALSE)

  peakbed <- xls_anno[,1:6]
  colnames(peakbed)[1] <- "# chr"

  peakbed$chromStart <- as.integer(peakbed$chromStart)
  peakbed$chromEnd <- as.integer(peakbed$chromEnd)
  write.table(peakbed, file =  paste(outfilepath, "/", xls_name, ".bed", sep = ""),
              sep = "\t", row.names = FALSE, quote = FALSE)
}


















