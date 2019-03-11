
.getsf <- function(exomepeak_path, sample_conditions) {

  ## for condition 1
  con <- sample_conditions[1]
  exomepeak_path0 <- unique(exomepeak_path[sample_conditions == con])
  sf1 <- .gettotalsf(exomepeak_path0)
  ind1 <- 1:(length(sf1)/2)

  ## for condition 2
  if (is.element("treated", sample_conditions)) {
    exomepeak_path0 <- unique(exomepeak_path[sample_conditions != con])
    sf2 <- .gettotalsf(exomepeak_path0)
    ind2 <- 1:(length(sf2)/2)

    sf_ip <- c(sf1[ind1], sf2[ind2])
    sf_input <- c(sf1[-ind1], sf2[-ind2])

    size_factor <- cbind(as.numeric(sf_ip), as.numeric(sf_input))
  } else {
    size_factor <- cbind(as.numeric(sf1[ind1]), as.numeric(sf1[-ind1]))
  }


  return(size_factor)
}

## deseq size factor

.getdeseqsizefactor <- function(exomepeak_path, sample_conditions) {

  ## for condition 1
  con <- sample_conditions[1]
  exomepeak_path0 <- unique(exomepeak_path[sample_conditions == con])
  sf_m1 <- .getsfm(exomepeak_path0)
  ind1 <- 1:(ncol(sf_m1)/2)

  ## for condition 2
  exomepeak_path0 <- unique(exomepeak_path[sample_conditions != con])
  sf_m2 <- .getsfm(exomepeak_path0)
  ind2 <- 1:(ncol(sf_m2)/2)

  ## get sf matrix
  sf_m <- cbind(sf_m1[,ind1], sf_m2[,ind2], sf_m1[,-ind1], sf_m2[,-ind2])
  sf_m <- sf_m[rowSums(sf_m) > ncol(sf_m),]

  ## get deseq size factor
  size_factor <- estimateSizeFactorsForMatrix(sf_m)

  return(as.numeric(size_factor))
}

.getsfm <- function(exomepeak_path0) {

  load(paste(exomepeak_path0, "exomePeak.Rdata", sep = "/"))
  sf_m <- tmp_rs$READS_COUNT
  len <- c(ncol(sf_m), ncol(sf_m) - 1)
  sf_m <- sf_m[,-len]
  sf_m_ip <- sf_m[,tmp_rs$SAMPLE_ID$untreated_ip]
  sf_m_input <- sf_m[,tmp_rs$SAMPLE_ID$untreated_input]
  sf_m <- cbind(sf_m_ip, sf_m_input)

  return(sf_m)
}

.gettotalsf <- function(exomepeak_path0) {

  load(paste(exomepeak_path0, "exomePeak.Rdata", sep = "/"))
  sf_m <- tmp_rs$READS_COUNT
  len <- c(ncol(sf_m), ncol(sf_m) - 1)
  sf_m <- sf_m[,-len]
  size_factor <- colSums(sf_m)
  sf_ip <- size_factor[tmp_rs$SAMPLE_ID$untreated_ip]
  sf_input <- size_factor[tmp_rs$SAMPLE_ID$untreated_input]

  return(c(sf_ip, sf_input))
}
