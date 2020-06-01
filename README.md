# DMDeepm6A1.0
A R package used to identify single base resolution m6A and differential m6A methylation site from MeRIP-seq data version 1.0.

Version: 1.0.3

Date: 2019-11-01

Author: Songyao Zhang zsynwpu@gmail.com, Jia Meng Jia.Meng@xjtlu.edu.cn

Maintainer: Songyao Zhang zsynwpu@gmail.com

The package is developed for the single base resolution m6A sites identification and differential analysis for MeRIP-seq data of two experimental conditions to unveil the dynamics in post-transcriptional regulation of the RNA methylome. This is to our knowledge the first tool to identify single base m6A and differential m6A methylation (DmM) sites using a deep learning and statistic test. The statistic test employed the same rhtest as exomePeak or alternatively used QNB which is more strict so that sacrifice some sensitivity. Please feel free to contact Songyao Zhang zsynwpu@gmail.com if you have any questions. Many thanks to Jia Meng for the exomePeak R package (Meng, Jia, et al. "Exome-based analysis for RNA epigenome sequencing data." Bioinformatics 29.12 (2013): 1565-1567.), which provides the base for DMDeepm6A.

License: GPL-2

Depends: Guitar, exomePeak, keras, QNB, DESeq, TxDb.Hsapiens.UCSC.hg19.knownGene, BSgenome.Hsapiens.UCSC.hg19, org.Hs.eg.db

# Installation

DMDeepm6A depends on Guitar, exomePeak, keras, QNB, DESeq, TxDb.Hsapiens.UCSC.hg19.knownGene, BSgenome.Hsapiens.UCSC.hg19, org.Hs.eg.dbr R packages and please make sure install them before installing DMDeepm6A.

1. Keras installation
Make sure Anaconda is installed for windows systerm

For Python 3.x (https://www.anaconda.com/download/#windows) before installing Keras.

Then follow the instruction on web (https://keras.rstudio.com/) to install keras or run following code in r:

install.packages("devtools")  
devtools::install_github("rstudio/keras")  ##install keras from github  
library(keras)  
install_keras()  

2.Other packages

install.packages("QNB")

if (!requireNamespace("BiocManager", quietly = TRUE))  
     install.packages("BiocManager")  
BiocManager::install("exomePeak", version = "3.8")  
BiocManager::install("Guitar", version = "3.8")  
BiocManager::install("DESeq", version = "3.8")  
BiocManager::install("BSgenome.Hsapiens.UCSC.hg19", version = "3.8")  
BiocManager::install("org.Hs.eg.db", version = "3.8")  

At last, DMDeepm6A can be installed as:

if (!requireNamespace("devtools", quietly = TRUE))    
     install.packages("devtools")    
devtools::install_github("NWPU-903PR/DMDeepm6A1.0")

# Toy Example diff sites calling

\# example for default hg19 genome

\# get input bam

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

\# get sample condition

sample_condition <- c("treated", "treated", "untreated", "untreated")

\# get genome annotation, generally, you can leave this default if you use the default hg19 genome

\# we use a toy gtf here to make this example run faster.

gft_genome <- system.file("extdata", "genes.gtf", package="DMDeepm6A")

\# diff peak calling

re <- dmdeepm6A(ip_bams = ip_bams,  
                input_bams = input_bams,  
                sample_conditions = sample_condition,    
                gft_genome = gft_genome)  

# Toy Example m6A site calling

\# do only peak calling for several replicates

\# get input bam

ip_bam1 <- system.file("extdata", "treated_ip1.bam", package="DMDeepm6A")  
ip_bam2 <- system.file("extdata", "treated_ip2.bam", package="DMDeepm6A")  
input_bam1 <- system.file("extdata", "treated_input1.bam", package="DMDeepm6A")  
input_bam2 <- system.file("extdata", "treated_input2.bam", package="DMDeepm6A")  

ip_bams <- c(ip_bam1, ip_bam2)  
input_bams <- c(input_bam1, input_bam2)

\# peak calling. sample_conditions should leave default (NA) or set all rep as "untreated".

re <- dmdeepm6A(ip_bams = ip_bams,  
                input_bams = input_bams,  
                gft_genome = gft_genome)  

# Example code for other species
\# An genome input formate example for rat rn5 genome

\# not run

\# get genome annotation

library(TxDb.Rnorvegicus.UCSC.rn5.refGene)  
txdb <- TxDb.Rnorvegicus.UCSC.rn5.refGene  

library(BSgenome.Rnorvegicus.UCSC.rn5)  
BSgenome <- BSgenome.Rnorvegicus.UCSC.rn5  

library(org.Rn.eg.db)  
egSYMBOL <- org.Rn.egSYMBOL  

sigpeak <- deepm6A(ip_bam,  
                   input_bam,  
                   sample_conditions = sample_conditions,  
                   DefaultGenome = FALSE,  
                   txdb = txdb,  
                   BSgenome = BSgenome,  
                   egSYMBOL = egSYMBOL)  




