#
# load("F:/DMDeepm6A/data/psuedoGene.RData")
#
# bedpath <- list.files("F:/DMDeepm6A/result/hesc/hEND8.0", pattern = "hEND8.0", full.names = T, recursive = F)
# sample_conditions <- c("untreated", "untreated", "treated", "treated")
#
# ## get bam input
# bams <- list.files("F:/DMDeepm6A/data/bam/hESC_hEND", pattern = "accepted_hits.bam$", recursive = T, full.names = T)
# ipbams <- bams[c(2,4,6,8)]
# inputbams <- bams[c(1,3,5,7)]
#
# sig_site_thresh <- 0.8
# sig_diff_thresh <- 0.1

## diff methy test function

dmtest <- function(ipbams,
                   inputbams,
                   bedpath,
                   tx_genes,
                   sample_conditions,
                   output_filepath = NA,
                   experiment_name = "DMtest_out",
                   sig_site_thresh = 0.8,
                   exomepeak_path = NA,
                   diff_method = "exomepeak",
                   diff_normalize = "TotalReads",
                   diff_bin_width = 201,
                   sig_diff_thresh = 0.01) {

  parameter <- list(tx_genes = tx_genes,
                    sample_conditions = sample_conditions,
                    output_filepath = output_filepath,
                    diff_method = diff_method,
                    diff_bin_width = diff_bin_width)
  if (is.na(output_filepath)) {
    output_filepath <- getwd()
    parameter$output_filepath <- output_filepath
  }

  output_filepath <- paste(output_filepath, experiment_name, sep = "/")
  parameter$output_filepath <- output_filepath
  if (!dir.exists(output_filepath)) {dir.create(output_filepath)}

  Diff <- list()
  Diff$bedpath <- bedpath

  ## get diff site index
  diffsite_ind <- .conind(bedpath, sample_conditions, sig_site_thresh)
  Diff$diffsite_ind <- diffsite_ind

  ## test diff methy for con sites
  bams <- c(ipbams, inputbams)
  bam_condition <- c(paste(sample_conditions, "ip", sep = "_"),
                     paste(sample_conditions, "input", sep = "_"))
  Diff$bams <- bams
  Diff$bam_condition <- bam_condition

  ## get size factor
  if (is.na(exomepeak_path[1])) {diff_normalize <- "DMSites"}
  if (diff_normalize == "DMSites") {sf <- NA}
  if (diff_normalize == "deseq") {sf <- .getdeseqsizefactor(exomepeak_path, sample_conditions)}
  if (diff_normalize == "TotalReads") {
    size_factor <- .getsf(exomepeak_path, sample_conditions)
    sf <- c(size_factor[,1], size_factor[,2])
  }

  Diff$sf <- sf
  parameter$Diff <- Diff

  diffsite_test_res <- .diffmethytest(parameter)
  parameter$Diff$diffsite_test_res <- diffsite_test_res

  ## annotate site
  site_anno <- .annosite(parameter, sig_site_thresh)
  diff_site <- site_anno$xls_diff
  diff_site <- diff_site[diff_site$diff.padj <= sig_diff_thresh,]
  diff_site <- diff_site[!is.na(diff_site[,1]),]

  .writeSiteAnno(site_anno$xls_hyper, output_filepath, "HyperMethySite")
  .writeSiteAnno(site_anno$xls_hypo, output_filepath, "HypoMethySite")
  .writeSiteAnno(site_anno$xls_diff, output_filepath, "CandidateDiffMethySite")
  .writeSiteAnno(diff_site, output_filepath, "SigDiffMethySite")

}
