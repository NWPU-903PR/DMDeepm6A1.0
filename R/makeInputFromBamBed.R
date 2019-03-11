

.getinput <- function(parameter) {

  ## get merippp
  if (length(parameter$exomepeak_input) == 1) {
    exomepeak_input <- .getexomepeakinput(parameter$peak_bed, parameter)
    parameter$exomepeak_input <- exomepeak_input
  }

  input <- .makeInput(parameter, extend = 50, min_reads = 5)
  input$exomePeak_filter <- parameter$exomepeak_input$id!=0

  save(input, file = paste(parameter$output_filepath, "input.RData", sep = "/"))

  return(input)
}

.getexomepeakinput <- function(peak_bed, parameter) {

  meripp <- BED12toGRangesList(peak_bed)
  id <- countOverlaps(meripp, parameter$tx_genes, type = "within")
  meripp <- meripp[id!=0]

  print("Cutting peak to Single Base...")
  tic <- Sys.time()
  peakSB <- .makePeak2SB(meripp, parameter$tx_genes, extend = 50)
  print(Sys.time() - tic)


  ## get motif info
  print("Searching peak for motifs...")
  tic <- Sys.time()
  peakMF <- .makePeak2Motif(meripp, parameter, extend = 50)
  print(Sys.time() - tic)

  return(list(meripp = meripp, id = id, peakSB = peakSB, peakMF = peakMF))
}
