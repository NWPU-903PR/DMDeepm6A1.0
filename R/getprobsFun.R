

## select model --------------------------------------------------------------

.modelselect <- function(filepath) {

  mod_list <- list.files(path = filepath,
                         pattern = "hdf5$", full.names = T, recursive = F)

  return(mod_list)

}



## sample input ----------------------------------------------------------------------------------

.sampleinput <- function(input) {

  ## get cnn input
  cnn_input <- input$cnn_input
  cnn_x_test <- lapply(cnn_input, t)
  cnn_x_test <- array_reshape(unlist(cnn_x_test),
                              c(length(cnn_x_test), ncol(cnn_x_test[[1]]), nrow(cnn_x_test[[1]]), 1))

  return(cnn_x_test)
}


## get model probs --------------------------------------------------------------------------------

.modelpredict <- function(model, test) {

  model <- load_model_hdf5(model, compile = F)

  ## compile model
  model %>% compile(
    loss = loss_categorical_crossentropy,
    optimizer = optimizer_adadelta(),
    metrics = c('accuracy')
  )

  ## predict
  y <- model %>% predict(test)
  probs <- y[,2]

  return(probs)
}


.combineprobs <- function(x) {
  n <- length(x[[1]])
  y <- matrix(unlist(x), nrow = n)
  y <- rowMeans(y)
  y <- data.frame(predict = y)
  return(y)
}

## get predict result -----------------------------------------------------------------------------

.getpredict <- function(input, filepath, sf, peak_xls, outputfile) {

  mod_list <- .modelselect(filepath)
  print("get input ...")
  cnn_input <- .sampleinput(input)

  print("get probs ...")
  probs <- lapply(as.list(mod_list), .modelpredict, cnn_input)
  pred_prob <- .combineprobs(probs)

  peakloci <- input$motif_A_loci

  meripp <- input$meripp
  anno <- mcols(meripp)
  gene <- anno$V4
  peakname <- gene[input$motif_ind]

  print("get motifs ...")
  motif <- input$motif_seq
  motif <- lapply(motif, "[", 49:53)

  motif <- DNAStringSet(motif)
  motif <- as.character(motif)

  ## get fold change
  ipreads <- exp(input$readsbymotif_ip)
  # ipreads <- colMeans(ipreads)
  ipreads <- ipreads*input$readsip$total_reads_count/(10^8)

  inputreads <- exp(input$readsbymotif_input)
  # inputreads <- colMeans(inputreads)
  inputreads <- inputreads*input$readsinput$total_reads_count/(10^8)

  fold_enrchment <- (ipreads/sf[1])/(inputreads/sf[2])
  fold_enrchment <- colMaxs(fold_enrchment)
  Probability <- pred_prob$predict

  ## get exomepeak infor
  exomepeak_peak <- read.table(peak_xls, header = T, stringsAsFactors = F)
  exomepeak_peak <- exomepeak_peak[input$exomePeak_filter,]
  ind <- input$motif_ind
  exome.lg.p <-  exomepeak_peak$lg.p[ind]
  exome.lg.fdr <- exomepeak_peak$lg.fdr[ind]
  exome.fold_enrchment <- exomepeak_peak$fold_enrchment[ind]
  fold_enrchment[fold_enrchment <= 1] <- exome.fold_enrchment[fold_enrchment <= 1]


  print("get results ...")
  xls <- data.frame(chr = as.character(seqnames(peakloci)),
                    chromStart = as.numeric(start(peakloci)) - 1,
                    chromEnd = as.numeric(end(peakloci)),
                    name = peakname,
                    score = pred_prob$predict,
                    strand = as.character(strand(peakloci)),
                    motif = motif,
                    peak.lg.p = exome.lg.p,
                    peak.lg.fdr = exome.lg.fdr,
                    Probability = Probability,
                    fold_enrchment = fold_enrchment)

  # id <- unique(peakloci)
  # ind <- findOverlaps(peakloci, id)
  #
  # probid <- pred_prob$predict + (1:length(pred_prob$predict)) * 10^-9
  # probs <- tapply(probid, subjectHits(ind), max)
  # ind <- match(probs, probid)
  #
  # xls <- xls[ind,]
  # peakloci <- peakloci[ind]

  if (!dir.exists(outputfile)) {dir.create(outputfile)}

  write.table(xls, paste(outputfile, "CandidateSingleBasePeak.xls", sep = "/"),
              col.names = T, row.names = F, quote = F, sep = "\t")
  return(xls)
}
















