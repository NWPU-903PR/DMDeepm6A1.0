\name{dmtest}
\alias{dmtest}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Differential methylation test
}
\description{
This function is used to do differential methylation test.
}
\usage{
dmtest(ipbams,
       inputbams,
       bedpath,
       tx_genes,
       sample_conditions,
       output_filepath = NA,
       experiment_name = "DMtest_out",
       sig_site_thresh = 0.8,
       exomepeak_path = NA,
       diff_method = "exomepeak",
       diff_normalize = "TotalReads",
       diff_bin_width = 201,
       sig_diff_thresh = 0.01)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ipbams}{
A vector of characters denote the file path of IP samples in bam formate
}
  \item{inputbams}{
A vector of characters denote the file path of input samples in bam formate, should be paired with the IP samples
}
  \item{bedpath}{
A vector of characters denote the file path of deepm6A result of each sample should be paired with the IP samples
}
  \item{tx_genes}{
The transcriptome used to do the analysis,  usually do not need to input
}
  \item{sample_conditions}{
A vector of characters denote the condition of IP and Input samples. Should be the same length as the ip_bams and the values should be "untreated" or "treated"
}
  \item{output_filepath}{
The file path where to output the result
}
  \item{experiment_name}{
The name of the experiment
}
  \item{exomepeak_path}{
A vector of characters denote the file paths where you save the pre-generated exomePeak result for each IP-Input sample paire
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
Bin width used to count reads for each candidate site to do the test
}
  \item{sig_diff_thresh}{
Threshold used to determine whether a diff methylation site is significant
}
}
\details{
This function can be used based on the result of "dmdeepm6A" function
}
\value{
The input for this function should be the result of "dmdeepm6A" function
}
\references{
FDMDeep-m6A: Identification and prioritization of functional differential methylation genes.
}
\author{
Songyao Zhang
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
# load("F:/DMDeepm6A/data/psuedoGene.RData")
#
# bedpath <- list.files("F:/DMDeepm6A/result/hesc/hEND8.0", pattern = "hEND8.0", full.names = T, recursive = F)
# sample_conditions <- c("untreated", "untreated", "treated", "treated")
#
# ## get bam input
# bams <- list.files("F:/DMDeepm6A/data/bam/hESC_hEND", pattern = "accepted_hits.bam$", recursive = T, full.names = T)
# ipbams <- bams[c(2,4,6,8)]
# inputbams <- bams[c(1,3,5,7)]
#
# sig_site_thresh <- 0.8
# sig_diff_thresh <- 0.1

## diff methy test function
## dmtest(ipbams, inputbams, bedpath, tx_genes, sample_conditions, output_filepath = NA, sig_site_thresh = 0.907,
##         diff_method = "QNB", QNB_bin_width = 401, sig_diff_thresh = 0.1)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Statistical Inference }% use one of  RShowDoc("KEYWORDS")
