

BED12toGRangesList <- function(filepath, header=FALSE) {
  
  # message
  print("Converting BED12 to GRangesList")
  print("It may take a few minutes")
  
  # read bed file
  a = read.table(filepath,sep="\t",header=header,stringsAsFactors =FALSE)
  # mcols_info =a[,13:length(a[1,])]
  a = a[,1:12]
  
  # get transcripts
  no_tx = length(a[,1])
  tx_id = 1:no_tx;
  tx_name = paste("line_",1:no_tx,sep="")
  tx_chrom = a[,1]
  tx_strand = a[,6]
  tx_start = a[,2]+1
  tx_end = a[,3]
  transcripts= data.frame(tx_id,tx_name,tx_chrom,tx_strand,tx_start,tx_end)
  head(transcripts)
  

  
  # get genes
  tx_name = tx_name
  gene_id = as.character(a[,4])
  gene_id[is.na(gene_id)]="NA"
  gene=data.frame(tx_name,gene_id)
  
  # 
  splicing <- lapply(1:no_tx, .spliceSingleTrans, a=a, tx_start=tx_start)
  splicing <- .combineListOfSplicing(splicing)
  
  # make txdb
  peaks = suppressWarnings(
    makeTxDb(transcripts=transcripts, 
             splicings=splicing,
             genes=gene))
  
  # generate GRangesList
  tx <- exonsBy(peaks, "tx",use.names=TRUE)
  mcols(tx) <- a
  
  return(tx)
}

.combineListOfSplicing <- function(t){
  
  a <- paste("t[[",1:length(t),"]]", sep="")
  a <- paste(a,collapse =",")
  a <- paste("rbind(",a,")",sep="")
  c <- parse(text=a)
  b <- suppressWarnings(eval(c))
  
  return(b)
}

.spliceSingleTrans <- function(i,a,tx_start) {
  tx = a[i,]
  tx_id = i
  exon_rank=1:as.integer(tx[10])
  
  # get start
  temp = as.integer(strsplit(as.character(tx[12]), ",")[[1]]) + tx_start[i]
  exon_start=temp
  
  # get end
  temp = as.integer(strsplit(as.character(tx[11]), ",")[[1]])
  temp2 = temp + exon_start - 1
  exon_end=temp2
  
  # get CDS
  cds_start = exon_start
  cds_end = exon_end
  
  # get data frame
  splicing_tx = data.frame(tx_id,exon_rank,exon_start,exon_end,cds_start,cds_end)
  return(splicing_tx)
}