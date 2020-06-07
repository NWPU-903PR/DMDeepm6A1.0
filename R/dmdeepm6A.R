

dmdeepm6A <- function(ip_bams,
                      input_bams,
                      output_filepath = NA,
                      experiment_name = "DMDeepm6A_out",
                      sample_conditions = NA,
                      exomepeak_path = NA,
                      gft_genome = NA,
                      model_filepath = NA,
                      default_genome = TRUE,
                      tx_genes = NA,
                      txdb = NA,
                      BSgenome = NA,
                      egSYMBOL = NA,
                      minimal_gene_length = 0,
                      sig_site_thresh = 0.8,
                      diff_method = "exomepeak",
                      diff_normalize = "TotalReads",
                      diff_bin_width = 201,
                      sig_diff_thresh = 0.01) {

  ## get parameter
  parameter <- list(
    ip_bams = ip_bams,
    input_bams = input_bams,
    output_filepath = output_filepath,
    experiment_name = experiment_name,
    sample_conditions = sample_conditions,
    exomepeak_path = exomepeak_path,
    gft_genome = gft_genome,
    model_filepath = model_filepath,
    default_genome = default_genome,
    tx_genes = tx_genes,
    txdb = txdb,
    BSgenome = BSgenome,
    egSYMBOL = egSYMBOL,
    minimal_gene_length = minimal_gene_length,
    sig_site_thresh = sig_site_thresh,
    diff_method = diff_method,
    diff_bin_width = diff_bin_width,
    sig_diff_thresh = sig_diff_thresh
  )
  motif <- c("GGACA", "GGACC", "GGACT", "AGACA", "AGACC", "AGACT",
             "GAACA", "GAACC", "GAACT", "AAACA", "AAACC", "AAACT",
             "TGACA", "TGACC", "TGACT", "TAACA", "TAACC", "TAACT")
  parameter$motif <- motif


  ## set parameter
  if (is.na(model_filepath)) {
    model_filepath <- system.file("extdata", package="DMDeepm6A")
  }
  if (is.na(output_filepath)) {output_filepath <- getwd()}
  output_filepath <- paste(output_filepath, experiment_name, sep = "/")
  if (!dir.exists(output_filepath)) {dir.create(output_filepath)}
  parameter$output_filepath <- output_filepath
  parameter$model_filepath <- model_filepath

  ## get genome
  if (!is.na(txdb)) {default_genome <- FALSE}
  if (default_genome) {
    txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
    BSgenome <- BSgenome.Hsapiens.UCSC.hg19
    egSYMBOL <- org.Hs.egSYMBOL
    parameter$txdb <- txdb
    parameter$BSgenome <- BSgenome
    parameter$egSYMBOL <- egSYMBOL
  }

  if (!is.na(gft_genome)) {
    txdb <- makeTxDbFromGFF(gft_genome, format = "gtf")
    parameter$txdb <- txdb
  }

  if (sum(sum(is.na(tx_genes))) == 1)
  {
    print("making transcriptome, this may take several minutes...")
    tx_genes <- .makePsuedoGeneFromTXDB(txdb, minimal_gene_length = 0)
    txnames <- names(tx_genes)
    txnames[is.na(txnames)] <- "Unk"
    names(tx_genes) <- txnames
    parameter$tx_genes <- tx_genes
  }

  save(tx_genes, file = paste(output_filepath, "psuedoGene.RData", sep = "/"))

  if (is.na(sample_conditions[1])) {sample_conditions <- rep("untreated", length(ip_bams))}
  parameter$sample_conditions <- sample_conditions


  ## get exomepeak peak
  if (is.na(exomepeak_path))
  {
    exomepeak(IP_BAM = ip_bams[sample_conditions == "untreated"],
              INPUT_BAM = input_bams[sample_conditions == "untreated"],
              TXDB = txdb,
              OUTPUT_DIR = parameter$output_filepath,
              EXPERIMENT_NAME = "exomepeak_untreated")
    exomepeak_path_un <- paste(output_filepath, "exomepeak_untreated", sep = "/")

    if(is.element("treated", sample_conditions)) {
      exomepeak(IP_BAM = ip_bams[sample_conditions == "treated"],
                INPUT_BAM = input_bams[sample_conditions == "treated"],
                TXDB = txdb,
                OUTPUT_DIR = parameter$output_filepath,
                EXPERIMENT_NAME = "exomepeak_treated")
      exomepeak_path_tr <- paste(output_filepath, "exomepeak_treated", sep = "/")

      exomepeak_path <- rep(exomepeak_path_un, length(sample_conditions))
      exomepeak_path[sample_conditions == "treated"] <- exomepeak_path_tr
    } else {
      exomepeak_path <- rep(exomepeak_path_un, length(sample_conditions))
    }
  }
  parameter$exomepeak_path <- exomepeak_path

  gc()


  ## get size factor
  size_factor <- .getsf(exomepeak_path, sample_conditions)
  parameter$size_factor <- size_factor

  exomepeak_input <- list()
  ## peak calling
  for (i in 1:length(ip_bams)) {

    print(paste("Detecting sites for IP/Input paire", i, "----------------------------"))
    ip_bam <- ip_bams[i]
    input_bam <- input_bams[i]
    sample_name <- paste(experiment_name, letters[i], sample_conditions[i], sep = "_")
    if(!is.element("treated", sample_conditions)) {
      sample_name <- paste(experiment_name, "Rep", letters[i], sep = "_")
    }

    if (i == 1) {
      exomepeak_input[[i]] <- .getexomepeakinput(paste(exomepeak_path[i], "peak.bed", sep = "/"), parameter)
    } else if (exomepeak_path[i] == exomepeak_path[i-1]) {
      exomepeak_input[[i]] <- exomepeak_input[[i-1]]
    } else {
      exomepeak_input[[i]] <- .getexomepeakinput(paste(exomepeak_path[i], "peak.bed", sep = "/"), parameter)
    }

    sig_peak <- deepm6A(IP_bam = ip_bam,
                        Input_bam = input_bam,
                        output_filepath = output_filepath,
                        experiment_name = sample_name,
                        exomepeak_path = exomepeak_path[i],
                        exomepeak_input = exomepeak_input[[i]],
                        model_filepath = model_filepath,
                        default_genome = FALSE,
                        tx_genes = tx_genes,
                        txdb = txdb,
                        BSgenome = BSgenome,
                        egSYMBOL = egSYMBOL,
                        minimal_gene_length = minimal_gene_length,
                        sig_site_thresh = sig_site_thresh,
                        size_factor = size_factor[i,])

    rm(sig_peak)
    gc()

  }

  Diff <- list()
  bedpath <- list.files(output_filepath, pattern = experiment_name, full.names = T, recursive = F)
  Diff$bedpath <- bedpath

  if(is.element("treated", sample_conditions)) {
    ## get diff site index
    diffsite_ind <- .conind(bedpath, sample_conditions, sig_site_thresh)

    ## test diff methy for con sites
    bams <- c(ip_bams, input_bams)
    bam_condition <- c(paste(sample_conditions, "ip", sep = "_"),
                       paste(sample_conditions, "input", sep = "_"))

    ## get size factor
    if (diff_normalize == "DMSites") {sf <- NA}
    if (diff_normalize == "deseq") {sf <- .getdeseqsizefactor(exomepeak_path, sample_conditions)}
    if (diff_normalize == "TotalReads") {
      size_factor <- .getsf(exomepeak_path, sample_conditions)
      sf <- c(size_factor[,1], size_factor[,2])
    }


    Diff$diffsite_ind <- diffsite_ind
    Diff$bams <- bams
    Diff$bam_condition <- bam_condition
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
  } else {
    site_anno <- .getconsensussite(bedpath)
    .writeSiteAnno(site_anno$xls, output_filepath, "SigMethySite")
    .writeSiteAnno(site_anno$xls_con, output_filepath, "ConSigMethySite")
  }


  return(site_anno)
}









