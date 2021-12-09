# DMDeepm6A1.0
A R package used to identify single base resolution m6A and differential m6A methylation site from MeRIP-seq data version 1.0.

Version: 1.0.6

Date: 2021-12-09

Author: Songyao Zhang zsynwpu@gmail.com, Jia Meng Jia.Meng@xjtlu.edu.cn

Maintainer: Songyao Zhang zsynwpu@gmail.com

The package is developed for the single base resolution m6A sites identification and differential analysis for MeRIP-seq data of two experimental conditions to unveil the dynamics in post-transcriptional regulation of the RNA methylome. This is to our knowledge the first tool to identify single base m6A and differential m6A methylation (DmM) sites using a deep learning and statistic test. The statistic test employed the same rhtest as exomePeak or alternatively used QNB （QNB is currently unavailable, we will try to solve this soon） which is more strict so that sacrifice some sensitivity. Please feel free to contact Songyao Zhang zsynwpu@gmail.com if you have any questions. Many thanks to Jia Meng for the exomePeak R package (Meng, Jia, et al. "Exome-based analysis for RNA epigenome sequencing data." Bioinformatics 29.12 (2013): 1565-1567.), which provides the base for DMDeepm6A.

License: GPL-2

Depends: exomePeak, keras, DESeq, TxDb.Hsapiens.UCSC.hg19.knownGene, BSgenome.Hsapiens.UCSC.hg19, org.Hs.eg.db

Citation: Song-Yao Z , Shao-Wu Z , Xiao-Nan F , et al. FunDMDeep-m6A: identification and prioritization of functional differential m6A methylation genes [J]. Bioinformatics, 2019(14):14.

## 1. Installation

DMDeepm6A depends on exomePeak, keras, DESeq, TxDb.Hsapiens.UCSC.hg19.knownGene, BSgenome.Hsapiens.UCSC.hg19, org.Hs.eg.dbr R packages and please make sure install them before installing DMDeepm6A.

1.	Keras installation    
Make sure Anaconda is installed for windows system for Python 3.x (https://www.anaconda.com/download/#windows) before installing Keras ("Anaconda3-5.3.0-Windows-x86_64" is suggested). Then follow the instruction on web (https://keras.rstudio.com/) to install keras (please do not install the current version, tensorflow version 1.10 and keras version 2.2.0 is suggested) or install in R as following:

```{r, eval=FALSE}
if (!requireNamespace("devtools", quietly = TRUE))    
    install.packages("devtools")

devtools::install_version("reticulate", version = "1.10",repos = "https://cloud.r-project.org/")
devtools::install_version("tensorflow", version = "1.10",repos = "https://cloud.r-project.org/")
devtools::install_version("keras", version = "2.2.0",repos = "https://cloud.r-project.org/")

```

You may set your convenient CRAN mirrors by setting `repos`.    
Then create the required virtual environment and install corresponding version of keras and tensorflow (e.g., tensorflow1.10, keras2.2.0) using conda, for Windows users, you can set this in Anaconda Prompt:

```
conda create -n r-tensorflow tensorflow==1.10 keras==2.2.0 numpy==1.14.5
```

` r-tensorflow ` is the name of your virtual environment will be used by keras in R. For windows users, please make sure your Anaconda path is added to environment variable. Then open R to check whether the keras is available:

```
library(keras)
is_keras_available()

[1] TRUE
```
For Linux/ubuntu users, ` reticulate ` will firstly search ` ~/.virtualenvs ` for python, then to make sure your R use the python installed in your conda env, you may need to firstly setup the ` RETICULATE_PYTHON ` before library and use ` keras ` and ` DMDeepm6A `:

```
Sys.setenv(RETICULATE_PYTHON = "Path_To_conda/envs/r-tensorflow/bin/python")
library(keras)
is_keras_available()

[1] TRUE
```

2.	Other required Bioconductor packages
```{r, eval=FALSE}
if (!requireNamespace("BiocManager", quietly = TRUE))    
    install.packages("BiocManager")

BiocManager::install(c("exomePeak", "DESeq", "TxDb.Hsapiens.UCSC.hg19.knownGene",    
                       "BSgenome.Hsapiens.UCSC.hg19", "org.Hs.eg.db"))
```

The "exomePeak" and "DESeq" packages are not available in Bioconductor3.10 for R>=4.0. If you are using R4.0 or above, Please install "exomePeak" from https://github.com/ZW-xjtlu/exomePeak. And install the "DESeq" package by:

```{r, eval=FALSE}
##for Linux installation
install.packages("https://www.bioconductor.org/packages/3.11/bioc/src/contrib/DESeq_1.39.0.tar.gz", repos = NULL, type="source")
##for Windows installation
install.packages("https://www.bioconductor.org/packages/3.11/bioc/bin/windows/contrib/4.0/DESeq_1.39.0.zip", repos = NULL, type="source")
```

3.	DMDeepm6A installation    
```
if (!requireNamespace("devtools", quietly = TRUE))    
    install.packages("devtools")

devtools::install_github("NWPU-903PR/DMDeepm6A1.0")
```

## 2. Toy Example diff sites calling

\# example for default hg19 genome

\# get input bam

```
Sys.setenv(RETICULATE_PYTHON = "Path_To_conda/envs/r-tensorflow/bin/python") ## For Linux/ubuntu
library(DMDeepm6A)
library(keras)

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
```

\# get sample condition

```
sample_condition <- c("treated", "treated", "untreated", "untreated")
```

\# get genome annotation, generally, you can leave this default if you use the default hg19 genome

\# we use a toy gtf here to make this example run faster.

```
gft_genome <- system.file("extdata", "genes.gtf", package="DMDeepm6A")
```

\# diff peak calling

```
re <- dmdeepm6A(ip_bams = ip_bams,  
                input_bams = input_bams,  
                sample_conditions = sample_condition,    
                gft_genome = gft_genome)  
```
Please note that, if you are not using default hg19 genome, you need to input at least both the genome annotation (txdb or gft_genome) and sequence (BSgenome). Please see section 4 for more details.

## 3. Toy Example m6A site calling

\# do only peak calling for several replicates

\# get input bam

```
ip_bam1 <- system.file("extdata", "treated_ip1.bam", package="DMDeepm6A")  
ip_bam2 <- system.file("extdata", "treated_ip2.bam", package="DMDeepm6A")  
input_bam1 <- system.file("extdata", "treated_input1.bam", package="DMDeepm6A")  
input_bam2 <- system.file("extdata", "treated_input2.bam", package="DMDeepm6A")  

ip_bams <- c(ip_bam1, ip_bam2)  
input_bams <- c(input_bam1, input_bam2)
```

\# peak calling. sample_conditions should leave default (NA) or set all rep as "untreated".

```
re <- dmdeepm6A(ip_bams = ip_bams,  
                input_bams = input_bams,  
                gft_genome = gft_genome)  
```
Please note that, if you are not using default hg19 genome, you need to input at least both the genome annotation (txdb or gft_genome) and sequence (BSgenome). Please see section 4 for more details.

## 4. Example code for other species
\# An genome input formate example for rat rn5 genome

\# not run

\# get genome annotation

```
library(TxDb.Rnorvegicus.UCSC.rn5.refGene)  
txdb <- TxDb.Rnorvegicus.UCSC.rn5.refGene  

library(BSgenome.Rnorvegicus.UCSC.rn5)  
BSgenome <- BSgenome.Rnorvegicus.UCSC.rn5  
```
` BSgenome ` is necessary if the genome is not hg19. If there is no public BSgenome data package for your own genome sequence, you can read your ".fa/.fasta" genome sequence file using readDNAStringSet function (e.g., BSgenome = readDNAStringSet(MyGenomeSeq.fa)), or forge a BSgenome data package following the instruction of "BSgenome" package.
```
library(org.Rn.eg.db)  
egSYMBOL <- org.Rn.egSYMBOL  
```
You can set ` egSYMBOL ` as NA, if there is no public egSYMBOL data package for your own species.
```
sigpeak <- deepm6A(ip_bam,  
                   input_bam,  
                   sample_conditions = sample_conditions,  
                   DefaultGenome = FALSE,  
                   txdb = txdb,  
                   BSgenome = BSgenome,  
                   egSYMBOL = egSYMBOL)  
```
