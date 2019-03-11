\name{deepm6A}
\alias{deepm6A}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Single base m6A dientification
}
\description{
This function is used to predict m6A site in single base resolution from MeRIP-Seq (m6A-Seq) data. The main features of the function includes:

1. Identification m6A site in single base resolution from exomePeak detected peak region.

2. Annotation the m6A site located gene name and transcript position, e.g. 5'UTR, CDS, 3'UTR or LncRNA.
}
\usage{
deepm6A(IP_bam,
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
        size_factor = NA)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{IP_bam}{
The file path of an IP sample in bam formate from MeRIP-Seq data.
}
  \item{Input_bam}{
The file path of an Input sample in bam formate from MeRIP-Seq data.
}
  \item{output_filepath}{
The file path where to output the result.
}
  \item{experiment_name}{
The name of the experiment.
}
  \item{exomepeak_path}{
The file path where you save the pre-generated exomePeak result.
}
  \item{exomepeak_input}{
The candidate single base m6A information extracted from exomepeak identified peaks, which is used for multiple replicates.
}
  \item{gft_genome}{
The gtf genome annotation, this is should be consistent with the BSgenome.
}
  \item{model_filepath}{
The CNN models used to do site prediction, usually do not need to input, only if you trained several CNN model by yourself.
}
  \item{default_genome}{
A logical parameter denotes whether use the default genome, default is TRUE. The default genome is hg19 for human.
}
  \item{tx_genes}{
The transcriptome used to do the analysis,  usually do not need to input.
}
  \item{txdb}{
The txdb famate genome annotation. You need to input it similar to "TxDb.Hsapiens.UCSC.hg19.knownGene" if DefaultGenome is FALSE.
}
  \item{BSgenome}{
The sequence of genome. You need to input it similar to "BSgenome.Hsapiens.UCSC.hg19" if DefaultGenome is FALSE.
}
  \item{egSYMBOL}{
The gene name annotation. You need to input it similar to "org.Hs.egSYMBOL" if DefaultGenome is FALSE.
}
  \item{minimal_gene_length}{
The minimal length of a trascript when make transcriptome from genome.
}
  \item{sig_site_thresh}{
The probability treshold used to define a detected single base m6A site as significant.
}
  \item{size_factor}{
The size_factor used to normalize the readscount of IP and Input samples when calculate the foldchange.
}
}
\details{
deepm6A is used to detect m6A site using one replicate, if you want do m6A site detection for multi replicates, please use dmdeepm6A. See example of dmdeepm6A function.

tx_genes is generated from the TXDB annotation genome you are using. It is saved as "psuedoGene.RData" under the output filepath. You can use it next time when you are using the same genome to save some time of making it from TXDB.

You need to input txdb, BSgenome and egSYMBOL if you do not want to use the default hg19 genome. Generally, you need to install this 3 genome information before you use them.
}
\value{
By default, deepm6A will output results both

1. as BED/XLS files on disk (default: "Deepm6A_out") under the specified directory (default: current working directory).

2. returned the annotated m6A site identified in single base resolution as data.frame object under the R environment.
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
ip_bam <- system.file("extdata", "untreated_ip1.bam", package="DMDeepm6A")
input_bam <- system.file("extdata", "untreated_input1.bam", package="DMDeepm6A")

## get genome annotation, generally, you can leave this default if you use the default hg19 genome
## we use a toy gtf here to make this example run faster.
gft_genome <- system.file("extdata", "genes.gtf", package="DMDeepm6A")

## peak calling
re <- deepm6A(IP_bam = ip_bam,
              Input_bam = input_bam,
              gft_genome = gft_genome)


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
#                    DefaultGenome = FALSE,
#                    txdb = txdb,
#                    BSgenome = BSgenome,
#                    egSYMBOL = egSYMBOL)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Statistical Inference }
