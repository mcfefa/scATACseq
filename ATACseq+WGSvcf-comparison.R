## Goal: Compare VCF from WGS to peaks identified in ATACseq
##
## Following along with: https://nbis-workshop-epigenomics.readthedocs.io/en/latest/content/tutorials/atac-chip-downstream/lab-PeakAnnot.html#introduction
##
setwd("~/Dropbox (UFL)/GitHub/scATACseq")

## installing packages needed -- only do when first setting up environment
install.packages('BiocManager', lib="~/Dropbox (UFL)/GitHub/scATACseq/lib")
BiocManager::install('ChIPseeker', lib="~/Dropbox (UFL)/GitHub/scATACseq/lib")
BiocManager::install('TxDb.Hsapiens.UCSC.hg37.knownGene', lib="~/Dropbox (UFL)/GitHub/scATACseq/lib")
BiocManager::install('GenomicAlignments', lib="~/Dropbox (UFL)/GitHub/scATACseq/lib")
BiocManager::install('GenomicFeatures', lib="~/Dropbox (UFL)/GitHub/scATACseq/lib")
BiocManager::install('clusterProfiler', lib="~/Dropbox (UFL)/GitHub/scATACseq/lib")
BiocManager::install('biomaRt', lib="~/Dropbox (UFL)/GitHub/scATACseq/lib")
BiocManager::install('org.Hs.eg.db', lib="~/Dropbox (UFL)/GitHub/scATACseq/lib")
BiocManager::install('ReactomePA', lib="~/Dropbox (UFL)/GitHub/scATACseq/lib")

## extra downloading to do 
BiocManager::install('TxDb.Hsapiens.UCSC.hg38.knownGene', lib="~/Dropbox (UFL)/GitHub/scATACseq/lib")

## load packages
library('BiocManager', lib="~/Dropbox (UFL)/GitHub/scATACseq/lib")

library('GenomicAlignments', lib="~/Dropbox (UFL)/GitHub/scATACseq/lib")
library('GenomicFeatures', lib="~/Dropbox (UFL)/GitHub/scATACseq/lib")
library('biomaRt', lib="~/Dropbox (UFL)/GitHub/scATACseq/lib")
library('org.Hs.eg.db', lib="~/Dropbox (UFL)/GitHub/scATACseq/lib")
library('TxDb.Hsapiens.UCSC.hg38.knownGene', lib="~/Dropbox (UFL)/GitHub/scATACseq/lib")
library('purrr', lib="~/Dropbox (UFL)/GitHub/scATACseq/lib")
library('withr', lib="~/Dropbox (UFL)/GitHub/scATACseq/lib")
library('farver', lib="~/Dropbox (UFL)/GitHub/scATACseq/lib")
library('clusterProfiler', lib="~/Dropbox (UFL)/GitHub/scATACseq/lib")
## to address error when loading, had to run: install.packages("withr", lib="~/Dropbox (UFL)/GitHub/scATACseq/lib")
## to address error when loading, had to run: BiocManager::install("MPO.db", lib="~/Dropbox (UFL)/GitHub/scATACseq/lib")
## to address error when loading, had to run: install.packages("farver", lib="~/Dropbox (UFL)/GitHub/scATACseq/lib")

library('reactome.db', lib="~/Dropbox (UFL)/GitHub/scATACseq/lib")
library('ReactomePA', lib="~/Dropbox (UFL)/GitHub/scATACseq/lib")
## to address error when loading, had to run: BiocManager::install("HDO.db", lib="~/Dropbox (UFL)/GitHub/scATACseq/lib")
## to address error when loading, had to run: BiocManager::install("HPO.db", lib="~/Dropbox (UFL)/GitHub/scATACseq/lib")
## to address error when loading, had to run: install.packages("tidyverse", lib="~/Dropbox (UFL)/GitHub/scATACseq/lib")
## to address error when loading, had to run: install.packages("purrr", lib="~/Dropbox (UFL)/GitHub/scATACseq/lib")
## to address error when loading, had to run: BiocManager::install("reactome.db", lib="~/Dropbox (UFL)/GitHub/scATACseq/lib")


library('ChIPseeker', lib="~/Dropbox (UFL)/GitHub/scATACseq/lib")
## to address error when loading, had to run: install.packages("memoise")
## to address error when loading, had to run: BiocManager::install("TxDb.Hsapiens.UCSC.hg19.knownGene", lib="~/Dropbox (UFL)/GitHub/scATACseq/lib")







