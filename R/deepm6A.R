
deepm6A <- function(IP_bam,
                    Input_bam,
                    output_filepath = NA,
                    experiment_name = "Deepm6A_out",
                    exomepeak_path = NA,
                    exomepeak_input = NA,
                    gft_genome = NA,
                    model_filepath = NA,
                    default_genome = TRUE,
                    tx_genes = NA,
                    txdb = NA,
                    BSgenome = NA,
                    egSYMBOL = NA,
                    minimal_gene_length = 0,
                    sig_site_thresh = 0.907,
                    size_factor = NA) {

  ## get parameter
  parameter <- list(
    ip_bam = IP_bam,
    input_bam = Input_bam,
    output_filepath = output_filepath,
    experiment_name = experiment_name,
    exomepeak_path = exomepeak_path,
    exomepeak_input = exomepeak_input,
    gft_genome = gft_genome,
    model_filepath = model_filepath,
    default_genome = default_genome,
    tx_genes = tx_genes,
    txdb = txdb,
    BSgenome = BSgenome,
    egSYMBOL = egSYMBOL,
    minimal_gene_length = minimal_gene_length,
    sig_site_thresh = sig_site_thresh
  )
  motif <- c("GGACA", "GGACC", "GGACT", "AGACA", "AGACC", "AGACT",
             "GAACA", "GAACC", "GAACT", "AAACA", "AAACC", "AAACT",
             "TGACA", "TGACC", "TGACT", "TAACA", "TAACC", "TAACT")
  parameter$motif <- motif


  ## set parameter
  ip_bam <- IP_bam
  input_bam <- Input_bam
  if (is.na(model_filepath)) {
    model_filepath <- system.file("extdata", package="DMDeepm6A")
  }
  if (is.na(output_filepath)) {output_filepath <- getwd()}
  output_filepath <- paste(output_filepath, experiment_name, sep = "/")
  if (!dir.exists(output_filepath)) {dir.create(output_filepath)}
  parameter$output_filepath <- output_filepath
  parameter$model_filepath <- model_filepath

  ## get genome
  if (!(is.na(txdb) & is.na(gft_genome))) {
    if (is.na(BSgenome)) {stop("BSgenome should not be NA if the genome is not defalt hg19")}
    default_genome <- FALSE
  }
  if (default_genome) {
    txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
    BSgenome <- BSgenome.Hsapiens.UCSC.hg19
    egSYMBOL <- org.Hs.egSYMBOL
    parameter$txdb <- txdb
    parameter$BSgenome <- BSgenome
    parameter$egSYMBOL <- egSYMBOL
  }

  if (!is.na(gft_genome) & is.na(txdb)) {
    txdb <- makeTxDbFromGFF(gft_genome, format = "gtf")
    parameter$txdb <- txdb
  }

  ## get exomepeak peak
  if (is.na(exomepeak_path))
  {
    exomepeak(ip_bam, input_bam,
              TXDB = txdb,
              OUTPUT_DIR = parameter$output_filepath,
              EXPERIMENT_NAME = "exomepeak_out")
    exomepeak_path <- paste(output_filepath, "exomepeak_out", sep = "/")
    parameter$exomepeak_path <- exomepeak_path
  }


  ## get genome infor
  if (is.na(tx_genes[[1]][1]))
  {
    print("making transcriptomes, this may take several minutes...")
    tx_genes <- .makePsuedoGeneFromTXDB(txdb, minimal_gene_length = 0)
    txnames <- names(tx_genes)
    txnames[is.na(txnames)] <- "Unk"
    names(tx_genes) <- txnames
    parameter$tx_genes <- tx_genes
    save(tx_genes, file = paste(parameter$output_filepath, "psuedoGene.RData", sep = "/"))
  }


  ## get single base m6A sites

  ## get input
  peak_bed <- paste(exomepeak_path, "peak.bed", sep = "/")
  if (!file.exists(peak_bed)) {peak_bed <- paste(exomepeak_path, "diff_peak.bed", sep = "/")}
  parameter$peak_bed <- peak_bed

  print("Making input from peak region...")
  input <- .getinput(parameter)


  ## get sites
  modelfile <- model_filepath

  if (is.na(size_factor)) {
    load(paste(exomepeak_path, "exomePeak.Rdata", sep = "/"))
    size_factor <- colSums(tmp_rs$READS_COUNT[,1:2])
    size_factor <- as.numeric(size_factor)
  }

  peak_xls <- strsplit(peak_bed, "bed")
  peak_xls <- paste(peak_xls[[1]], "xls", sep = "")

  print("Detecting single base m6A sites...")
  peak <- .getpredict(input, modelfile, size_factor, peak_xls, output_filepath)

  ## anno peak
  sigthresh <- sig_site_thresh

  ## anno site position
  component <- .extractComponent(txdb)

  peak <- .getpeakposition(peak, component, txdb, egSYMBOL, output_filepath)
  sigpeak <- peak[peak$score >=sigthresh,]

  write.table(sigpeak, file =  paste(output_filepath, "SigSingleBasePeak.xls" , sep = "/"),
              sep = "\t", row.names = FALSE, quote = FALSE)

  peakbed <- peak[,1:6]
  colnames(peakbed)[1] <- "# chr"
  sigpeakbed <- peakbed[peak$score >=sigthresh,]

  write.table(peakbed, file =  paste(output_filepath, "CandidateSingleBasePeak.bed" , sep = "/"),
              sep = "\t", row.names = FALSE, quote = FALSE)
  write.table(sigpeakbed, file =  paste(output_filepath, "SigSingleBasePeak.bed" , sep = "/"),
              sep = "\t", row.names = FALSE, quote = FALSE)

  return(sigpeak)

}


