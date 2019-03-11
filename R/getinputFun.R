

## make input -----------------------------------------------------------------------------------


.makeInput <- function(parameter, extend = 50, min_reads = 5) {

  tx_genes <- parameter$tx_genes
  meripp <- parameter$exomepeak_input$meripp

  ## get bin of peak
  peakSB <- parameter$exomepeak_input$peakSB
  peakMF <- parameter$exomepeak_input$peakMF

  ## get bam reads count
  print("Counting reads for single base...")
  tic <- Sys.time()
  readsip <- .bamReadsCount(peakSB, parameter$ip_bam, tx_genes)
  readsinput <- .bamReadsCount(peakSB, parameter$input_bam, tx_genes)
  print(Sys.time() - tic)

  readsbypeak_ip <- .makeInputMat(readsip$base_reads_count, peakSB)
  readsbypeak_input <- .makeInputMat(readsinput$base_reads_count, peakSB)

  ## get motif reads
  print("Making input...")
  motif_A_loci <- peakMF$motif_A_loci
  motif_ind <- peakMF$motif_ind
  motif_peak_loci <- pmapToTranscripts(motif_A_loci, meripp[motif_ind], ignore.strand = T)

  ind <- rbind(motif_ind, start(motif_peak_loci))
  ind <- as.vector(ind)
  id <- rep(1:length(motif_peak_loci), each = 2)

  readsbymotif_ip <- tapply(ind, id, .getmotifreads, readsbypeak_ip, extend)
  readsbymotif_ip <- matrix(unlist(readsbymotif_ip),
                            nrow = (extend*2 + 1), ncol = length(readsbymotif_ip))

  readsbymotif_input <- tapply(ind, id, .getmotifreads, readsbypeak_input, extend)
  readsbymotif_input <- matrix(unlist(readsbymotif_input),
                               nrow = (extend*2 + 1), ncol = length(readsbymotif_input))


  ## get reads sizefactor
  sf <- readsip$total_reads_count

  ## normalize reads count
  ind0 <- as.character(strand(motif_A_loci)) == "-"
  readsbymotif_ip <- .normalizeReads(readsbymotif_ip, sf, min_reads, ind0)

  ## get input reads sizefactor
  sf <- readsinput$total_reads_count

  ## normalize reads count
  readsbymotif_input <- .normalizeReads(readsbymotif_input, sf, min_reads, ind0)

  ## get sequance encode
  motif_seq <- peakMF$motif_seq
  tic <- Sys.time()
  seq_input <- lapply(motif_seq, .seq2vec, extend)
  print(Sys.time() - tic)


  ## make input
  cnn_input <- list()
  for(i in 1:length(seq_input)) {
    input_mat <- seq_input[[i]]
    input_mat <- t(t(input_mat)*readsbymotif_ip[,i])
    cnn_input[[i]] <- input_mat
  }

  re <- list(cnn_input = cnn_input,
             readsbymotif_ip = readsbymotif_ip,
             readsbymotif_input = readsbymotif_input,
             motif_seq = motif_seq,
             motif_A_loci = motif_A_loci,
             motif_peak_loci = motif_peak_loci,
             motif_ind = motif_ind,
             readsip = readsip,
             readsinput = readsinput,
             peakSB = peakSB,
             meripp = meripp)
  return(re)
}


.makeInputMat <- function(x, peakSB) {
  readsbypeak <- split(x, names(unlist(peakSB)))
  readsbypeak <- readsbypeak[names(peakSB)]
  return(readsbypeak)
}

.getmotifreads <- function(x, reads, extend) {
  reads[[x[1]]][x[2]:(x[2] + 2*extend)]
}

.normalizeReads <- function(reads, sf, min_reads, ind) {
  reads[reads < min_reads] <- min_reads
  reads[,ind] <- reads[nrow(reads):1, ind]
  reads <- round(reads/sf*10^8)
  reads <- log(reads)
  return(reads)
}

.seq2vec <- function(x, extend) {
  a <- matchPattern("A", x)
  a <- start(a)

  t <- matchPattern("T", x)
  t <- start(t)

  c <- matchPattern("C", x)
  c <- start(c)

  g <- matchPattern("G", x)
  g <- start(g)

  m <- matrix(0, nrow = 4, ncol = (2*extend + 1))
  m[1, a] <- 1
  m[2, t] <- 1
  m[3, c] <- 1
  m[4, g] <- 1
  return(m)
}
