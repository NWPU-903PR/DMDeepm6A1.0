



.makePsuedoGeneFromTXDB <- function(txdb, minimal_gene_length = 0) {

  # get exons
  tx_exons <- exonsBy(txdb, by="tx", use.name=FALSE)

  # get cluster info
  cluster_ID_by_overlap <- .findTranscriptCluster(tx_exons)
  txid <- as.numeric(names(tx_exons))
  tx2gene <- .tx2gene(txdb)
  tx_gene <- names(tx2gene)[txid]
  combined_id <- paste(tx_gene,cluster_ID_by_overlap)
  cluster_name <- unique(combined_id)
  cluster_id <- 1:length(cluster_name)
  names(cluster_id) <- cluster_name
  tx_cluster_id <- as.numeric(cluster_id[combined_id])
  tx_cluster_index <- split(1:length(tx_cluster_id),tx_cluster_id)

  # find cluster genes
  cluster2tx <- match(unique(tx_cluster_id),tx_cluster_id)
  cluster_gene <- tx_gene[cluster2tx]
  cluster_annotation <- data.frame(cluster_index=names(tx_cluster_index), gene=cluster_gene)

  # construct psuedo gene transcripts
  psuedoGene <- lapply(tx_cluster_index, .constructPsuedoGene, tx_exons=tx_exons)

  # generate check points
  psuedoGene <- GRangesList(c(unlist(psuedoGene)))
  names(psuedoGene) <- as.character(cluster_annotation[[2]])
  ps_name <- names(psuedoGene)
  names(psuedoGene)[is.na(ps_name)] <- paste("Unk", 1:sum(is.na(ps_name)))


  ## gene level
  # get gene cluster info
  cluster_ID_by_overlap <- .findTranscriptCluster(psuedoGene)
  tx_gene <- names(psuedoGene)
  combined_id <- paste(tx_gene,cluster_ID_by_overlap)
  cluster_name <- unique(combined_id)
  cluster_id <- 1:length(cluster_name)
  names(cluster_id) <- cluster_name
  tx_cluster_id <- as.numeric(cluster_id[combined_id])
  tx_cluster_index <- split(1:length(tx_cluster_id),tx_cluster_id)

  # find cluster genes
  cluster2tx <- match(unique(tx_cluster_id),tx_cluster_id)
  cluster_gene <- ps_name[cluster2tx]
  cluster_annotation <- data.frame(cluster_index=names(tx_cluster_index), gene=cluster_gene)

  # construct psuedo gene transcripts
  ind <- (lapply(tx_cluster_index, length) > 1)
  if (sum(ind) > 0) {
    psuedoGene_old <- psuedoGene[unlist(tx_cluster_index[!ind])]
    psuedoGene_new <- lapply(tx_cluster_index[ind], .constructPsuedoGene, tx_exons = psuedoGene)
    psuedoGene_new <- GRangesList(c(unlist(psuedoGene_new)))
    psuedoGene <- c(psuedoGene_old, psuedoGene_new)
    psuedoGene[ind] <- psuedoGene_new
    psuedoGene[!ind] <- psuedoGene_old
  }
  names(psuedoGene) <- as.character(cluster_annotation[[2]])

  # length filter
  gene_width_filter <- (sum(width(psuedoGene)) >= minimal_gene_length)

  # return result
  return(psuedoGene[gene_width_filter])

}

.findTranscriptCluster <- function(tx) {
  # find the hits
  hits <- findOverlaps(tx, tx)
  tx_hits <- selectHits(hits, select="first")

  old_tx_hits <- tx_hits + 1
  flag <- TRUE
  while (flag) {
    new_tx_hits <- tx_hits[tx_hits]
    old_tx_hits <- tx_hits
    tx_hits <- new_tx_hits
    flag <- (sum( abs( tx_hits - old_tx_hits ) )>0)
  }
  return(tx_hits)
}

.tx2gene <- function(txdb) {
  ID = keys(txdb, "TXID")
  temp = select(txdb, ID, c("TXID", "GENEID"), "TXID")
  txid <- temp[[1]]
  names(txid) <- temp[[2]]
  return(txid)
}

.constructPsuedoGene <- function(cluster_list, tx_exons) {
  psuedoGene <- reduce(unlist(tx_exons[cluster_list]))
}







