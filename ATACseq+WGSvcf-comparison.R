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
library('TxDb.Hsapiens.UCSC.hg19.knownGene', lib="~/Dropbox (UFL)/GitHub/scATACseq/lib")
library('ChIPseeker', lib="~/Dropbox (UFL)/GitHub/scATACseq/lib")
## to address error when loading, had to run: install.packages("memoise")
## to address error when loading, had to run: BiocManager::install("TxDb.Hsapiens.UCSC.hg19.knownGene", lib="~/Dropbox (UFL)/GitHub/scATACseq/lib")

## Pipeline
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

##### Read in BED file #####################################
pth2peaks_bed="./data/GSE96769-single-cell/GSE96769_PeakFile_20160207.bed.gz"
## below didn't work, so followed along with command from: https://bioconductor.org/packages/devel/bioc/vignettes/ChIPseeker/inst/doc/ChIPseeker.html#chip-profiling
#peaks.bed=read.table(pth2peaks_bed, sep="\t", header=FALSE, blank.lines.skip=TRUE, fill=TRUE)
#rownames(peaks.bed)=peaks.bed[,4]
#peaks.gr <- GRanges(seqnames=peaks.bed[,1], ranges=IRanges(peaks.bed[,2], peaks.bed[,3]), strand="*", mcols=data.frame(peakID=peaks.bed[,4]))

peakFromFile <- readPeakFile(pth2peaks_bed)
peaks.gr <- peakFromFile

##### To inspect peak coverage along the chromosomes:
covplot(peaks.gr, chrs=c("chr14", "chr15")). ####<----------

#to save the image to file
pdf("./results/PeakCoverage.pdf")
covplot(peaks.gr, chrs=c("chr14", "chr15"))
dev.off()

##### Peak Annotation
bed.annot = annotatePeak(peaks.gr, tssRegion=c(-3000, 3000),TxDb=txdb, annoDb="org.Hs.eg.db")

annot_peaks=as.data.frame(bed.annot)

write.table(annot_peaks, "./results/HSC_merged_annotated.txt",
            append = FALSE,
            quote = FALSE,
            sep = "\t",
            row.names = FALSE,
            col.names = TRUE,
            fileEncoding = "")

pdf("./results/AnnotVis.pdf")
upsetplot(bed.annot, vennpie=TRUE)
dev.off()


