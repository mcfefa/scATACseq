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
BiocManager::install('TxDb.Hsapiens.UCSC.hg38.knownGene', lib="~/Dropbox (UFL)/GitHub/scATACseq/lib")
BiocManager::install('GenomicRanges', lib="~/Dropbox (UFL)/GitHub/scATACseq/lib")


## load packages
library('BiocManager', lib="~/Dropbox (UFL)/GitHub/scATACseq/lib")

library('farver', lib="~/Dropbox (UFL)/GitHub/scATACseq/lib")
library('ChIPseeker', lib="~/Dropbox (UFL)/GitHub/scATACseq/lib")
## to address error when loading, had to run: install.packages("memoise")
## to address error when loading, had to run: BiocManager::install("TxDb.Hsapiens.UCSC.hg19.knownGene", lib="~/Dropbox (UFL)/GitHub/scATACseq/lib")

library('GenomicAlignments', lib="~/Dropbox (UFL)/GitHub/scATACseq/lib")
library('GenomicFeatures', lib="~/Dropbox (UFL)/GitHub/scATACseq/lib")
library('biomaRt', lib="~/Dropbox (UFL)/GitHub/scATACseq/lib")
library('org.Hs.eg.db', lib="~/Dropbox (UFL)/GitHub/scATACseq/lib")
library('TxDb.Hsapiens.UCSC.hg38.knownGene', lib="~/Dropbox (UFL)/GitHub/scATACseq/lib")
library('purrr', lib="~/Dropbox (UFL)/GitHub/scATACseq/lib")
library('withr', lib="~/Dropbox (UFL)/GitHub/scATACseq/lib")
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

## Pipeline
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene

##### Read in BED file #####################################
pth2peaks_bed="./data/GSE96769-single-cell/GSE96769_PeakFile_20160207.bed.gz"
## below didn't work, so followed along with command from: https://bioconductor.org/packages/devel/bioc/vignettes/ChIPseeker/inst/doc/ChIPseeker.html#chip-profiling
#peaks.bed=read.table(pth2peaks_bed, sep="\t", header=FALSE, blank.lines.skip=TRUE, fill=TRUE)
#rownames(peaks.bed)=peaks.bed[,4]
#peaks.gr <- GRanges(seqnames=peaks.bed[,1], ranges=IRanges(peaks.bed[,2], peaks.bed[,3]), strand="*", mcols=data.frame(peakID=peaks.bed[,4]))

peaks.gr <- readPeakFile(pth2peaks_bed, header=FALSE)

##### To inspect peak coverage along the chromosomes:
library('labeling', lib="~/Dropbox (UFL)/GitHub/scATACseq/lib")
covplot(peaks.gr, chrs=c("chr14", "chr15")) 
## Errored needing a package "labeling", so ran: install.packages("labeling", lib="~/Dropbox (UFL)/GitHub/scATACseq/lib")

#to save the image to file
pdf("./results/PeakCoverage_chr1-8_2023-11-05.pdf")
covplot(peaks.gr, chrs=c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6","chr7","chr8"))
dev.off()

pdf("./results/PeakCoverage_chr9-16_2023-11-05.pdf")
covplot(peaks.gr, chrs=c("chr9", "chr10", "chr11", "chr12", "chr13", "chr14","chr15","chr16"))
dev.off()

pdf("./results/PeakCoverage_chr17-X_2023-11-05.pdf")
covplot(peaks.gr, chrs=c("chr17", "chr18", "chr19", "chr20", "chr21", "chr22","chrX"))
dev.off()

saveRDS(peaks.gr, "./results/peaks-gr_largeGRanges-object_2023-11-05.rds")

##### Peak Annotation 
bed.annot_hg19 = annotatePeak(peakFromFile, tssRegion=c(-3000, 3000),TxDb=txdb, annoDb="org.Hs.eg.db")

annot_peaks_hg19=as.data.frame(bed.annot_hg19)

write.table(annot_peaks_hg19, "./results/HSC_merged_annotated_hg19_2023-11-06.txt",
            append = FALSE,
            quote = FALSE,
            sep = "\t",
            row.names = FALSE,
            col.names = TRUE,
            fileEncoding = "")

saveRDS(bed.annot_hg19, "./results/bed-annot_hg19_csAnno-object_2023-11-06.rds")
saveRDS(annot_peaks_hg19, "./results/annot-peaks_hg19_DataFrame_2023-11-06.rds")

## reload previous steps
peaks.gr <- readRDS("./results/peaks-gr_largeGRanges-object_2023-11-05.rds")
bed.annot_hg19 <- readRDS("./results/bed-annot_hg19_csAnno-object_2023-11-06.rds")
annot_peaks_hg19 <- readRDS("./results/annot-peaks_hg19_DataFrame_2023-11-06.rds")

## additional tutorial: https://www.bioconductor.org/packages/3.2/bioc/vignettes/ChIPseeker/inst/doc/ChIPseeker.html

##install.packages('GenomicRanges', lib="~/Dropbox (UFL)/GitHub/scATACseq/lib")
##install.packages('ggupset', lib="~/Dropbox (UFL)/GitHub/scATACseq/lib")
library('ggupset', lib="~/Dropbox (UFL)/GitHub/scATACseq/lib")
library('labeling', lib="~/Dropbox (UFL)/GitHub/scATACseq/lib")
pdf("./results/AnnotVis_2023-11-06.pdf")
upsetplot(bed.annot_hg19, vennpie=TRUE)
dev.off()

pdf("./results/Vis_AnnoPie_2023-11-06.pdf")
plotAnnoPie(bed.annot_hg19)
dev.off()

pdf("./results/Vis_vennpie_2023-11-06.pdf")
vennpie(bed.annot_hg19)
dev.off()

pdf("./results/Vis_plotAnnoBar_2023-11-06.pdf")
plotAnnoBar(bed.annot_hg19)
dev.off()

######### struggling to get upsetplot working, not sure why, works with example
#########   putting on hold for now and can revisit later
#########
# #install.packages("UpSetR", lib="~/Dropbox (UFL)/GitHub/scATACseq/lib")
# library("UpSetR", lib="~/Dropbox (UFL)/GitHub/scATACseq/lib")
# ####<----------
# 
# ## bed.annot_hg19c = as.GRanges(bed.annot_hg19)[,c(6,7,12:23)]
# ## bed.annot_hg19d = annotatePeak(bed.annot_hg19c, TxDb=TxDb_hg19)
# ## attempt trimming to same order of metadata columns
# 
# pdf("./results/Vis_upsetplot_2023-11-07.pdf")
# upsetplot(bed.annot_hg19)
# dev.off()
# 
# pdf("./results/Vis_upsetplot+venn_2023-11-07.pdf")
# upsetplot(bed.annot_hg19d, vennpie=TRUE)
# dev.off()

## upsetplot errored needing a package "ggupset", so ran: install.packages("ggupset", lib="~/Dropbox (UFL)/GitHub/scATACseq/lib")

## Distribution of loci with respect to TSS:
pdf("./results/TSSdist_2023-11-07.pdf")
plotDistToTSS(bed.annot_hg19, title="Distribution of ATAC-seq peaks loci\nrelative to TSS")
dev.off()

##### Functional Analysis ####
##########
########## functional analysis is also not working, can revisit later
##########
# ## Reactome pathway enrichment of genes defined as the nearest feature to the peaks:
# #finding enriched Reactome pathways using chromosome 1 and 2 genes as a background
# pathway.reac <- enrichPathway(as.data.frame(annot_peaks)$geneId)
# 
# #previewing enriched Reactome pathways
# head(pathway.reac)
# 
# ## Let’s search for enriched GO terms:
# pathway.GO <- enrichGO(as.data.frame(annot_peaks)$geneId, org.Hs.eg.db, ont = "MF")
# 
# 



