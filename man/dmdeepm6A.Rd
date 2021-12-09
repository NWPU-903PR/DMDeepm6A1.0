\name{dmdeepm6A}
\alias{dmdeepm6A}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Single base m6A and DM m6A dientification
}
\description{
This function is main function of DMDeepm6A which is used to predict m6A site and DM m6A site in single base resolution from MeRIP-Seq (m6A-Seq) data. The main features of the function includes:

1. Identification m6A site in single base resolution from exomePeak detected peak region.

2. Annotation the m6A site located gene name and transcript position, e.g. 5'UTR, CDS, 3'UTR or LncRNA.

3. Identify differential methylation m6A sites in single base resolution.
}
\usage{
dmdeepm6A(ip_bams,
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
          sig_diff_thresh = 0.01)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ip_bams}{
A vector of characters denote the file path of IP samples in bam formate
}
  \item{input_bams}{
A vector of characters denote the file path of input samples in bam formate, should be paired with the IP samples
}
  \item{output_filepath}{
The file path where to output the result
}
  \item{experiment_name}{
The name of the experiment
}
  \item{sample_conditions}{
A vector of characters denote the condition of IP and Input samples. Should be the same length as the ip_bams and the values should be "untreated" or "treated"
}
  \item{exomepeak_path}{
A vector of characters denote the file paths where you save the pre-generated exomePeak result for each IP-Input sample paire
}
  \item{gft_genome}{
The gtf genome annotation, this is should be consistent with the BSgenome.
}
  \item{model_filepath}{
The CNN models used to do site prediction, usually do not need to input, only if you trained several CNN model by yourself
}
  \item{default_genome}{
A logical parameter denotes whether use the default genome, default is TRUE. The default genome is hg19 for human
}
  \item{tx_genes}{
The transcriptome used to do the analysis,  usually do not need to input
}
  \item{txdb}{
The txdb famate genome annotation. You need to input it similar to "TxDb.Hsapiens.UCSC.hg19.knownGene" if DefaultGenome is FALSE.
}
  \item{BSgenome}{
The sequence of genome. You need to input it similar to "BSgenome.Hsapiens.UCSC.hg19" if DefaultGenome is FALSE. If there is no public BSgenome data package for your genome sequence, you can read your ".fa/.fasta" genome sequence file using readDNAStringSet function (e.g., BSgenome = readDNAStringSet(MyGenomeSeq.fa)), or forge a BSgenome data package following the instruction of "BSgenome" package.
}
  \item{egSYMBOL}{
The gene name annotation. You need to input it similar to "org.Hs.egSYMBOL" if DefaultGenome is FALSE.
}
  \item{minimal_gene_length}{
The minimal length of a trascript when make transcriptome from genome
}
  \item{sig_site_thresh}{
The probability treshold used to define a detected single base m6A site as significant
}
  \item{diff_method}{
The method used to do differential methylation test for candidate DmM sites, the defalt is the same as exomePeak, and QNB is an alternative which will improve the precision but sacrifice the sensitivity a lot.
}
  \item{diff_normalize}{
The normalization used to do QNB test for differential methylation site. The default is "DMSites" which means the readscount will be normalized with the readscount of all candidate DM sites. Alternatively, it can be set as "TotalReads" or "deseq" which means the normalization will be done based on total number of reads of each replicate or as discribed in DeSeq method.
}
  \item{diff_bin_width}{
The bin width used to extend candidate DmM sites to count the reads.
}
  \item{sig_diff_thresh}{
The thresh used to identify whether a candidate DmM site is significantly differential methylated.
}
}
\details{
tx_genes is generated from the TXDB annotation genome you are using. It is saved as "psuedoGene.RData" under the output filepath. You can use it next time when you are using the same genome to save some time of making it from TXDB.

You need to input txdb, BSgenome and egSYMBOL if you do not want to use the default hg19 genome. Generally, you need to install this 3 genome information before you use them.
}
\value{
By default, deepm6A will output results both

1. as BED/XLS files on disk (default: "Deepm6A_out") under the specified directory (default: current working directory).

2. returned the annotated DM m6A site identified in single base resolution as list object under the R environment.
}
\references{
FDMDeep-m6A: Identification and prioritization of functional differential methylation genes.
}
\author{
Songyao Zhang
}
\examples{

## example for default hg19 genome

## get input bam
ip_bam1 <- system.file("extdata", "treated_ip1.bam", package="DMDeepm6A")
ip_bam2 <- system.file("extdata", "treated_ip2.bam", package="DMDeepm6A")
ip_bam3 <- system.file("extdata", "untreated_ip1.bam", package="DMDeepm6A")
ip_bam4 <- system.file("extdata", "untreated_ip2.bam", package="DMDeepm6A")
input_bam1 <- system.file("extdata", "treated_input1.bam", package="DMDeepm6A")
input_bam2 <- system.file("extdata", "treated_input2.bam", package="DMDeepm6A")
input_bam3 <- system.file("extdata", "untreated_input1.bam", package="DMDeepm6A")
input_bam4 <- system.file("extdata", "untreated_input2.bam", package="DMDeepm6A")

ip_bams <- c(ip_bam1, ip_bam2, ip_bam3, ip_bam4)
input_bams <- c(input_bam1, input_bam2, input_bam3, input_bam4)

## get sample condition
sample_condition <- c("treated", "treated", "untreated", "untreated")

## get genome annotation, generally, you can leave this default if you use the default hg19 genome
## we use a toy gtf here to make this example run faster.
gft_genome <- system.file("extdata", "genes.gtf", package="DMDeepm6A")

## diff peak calling
re <- dmdeepm6A(ip_bams = ip_bams,
                input_bams = input_bams,
                sample_conditions = sample_condition,
                gft_genome = gft_genome,
                BSgenome = BSgenome.Hsapiens.UCSC.hg19,
                egSYMBOL = org.Hs.egSYMBOL)

## do only peak calling for several replicates

## get input bam
ip_bam1 <- system.file("extdata", "treated_ip1.bam", package="DMDeepm6A")
ip_bam2 <- system.file("extdata", "treated_ip2.bam", package="DMDeepm6A")
input_bam1 <- system.file("extdata", "treated_input1.bam", package="DMDeepm6A")
input_bam2 <- system.file("extdata", "treated_input2.bam", package="DMDeepm6A")

ip_bams <- c(ip_bam1, ip_bam2)
input_bams <- c(input_bam1, input_bam2)

## peak calling. sample_conditions should leave default (NA) or set all rep as "untreated".
re <- dmdeepm6A(ip_bams = ip_bams,
                input_bams = input_bams,
                gft_genome = gft_genome,
                BSgenome = BSgenome.Hsapiens.UCSC.hg19,
                egSYMBOL = org.Hs.egSYMBOL)

# An genome input formate example for rat rn5 genome
# not run

## get genome annotation
# library(TxDb.Rnorvegicus.UCSC.rn5.refGene)
# txdb <- TxDb.Rnorvegicus.UCSC.rn5.refGene

# library(BSgenome.Rnorvegicus.UCSC.rn5)
# BSgenome <- BSgenome.Rnorvegicus.UCSC.rn5

# library(org.Rn.eg.db)
# egSYMBOL <- org.Rn.egSYMBOL

# sigpeak <- deepm6A(ip_bam,
#                    input_bam,
#                    sample_conditions = sample_conditions,
#                    DefaultGenome = FALSE,
#                    txdb = txdb,
#                    BSgenome = BSgenome,
#                    egSYMBOL = egSYMBOL)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Statistical Inference }
