### Following along with Current Protocols Approach

############################################################################
## Basic Protocol 1
####### ChIPseeker AND EPIGENOMIC DATASET PREPARATION ######################

download.file("https://raw.githubusercontent.com/YuLab-SMU/ChIPseeker_current_protocols/master/CP_demo_data/CTCF_H1.test.bed", destfile = "./test/CTCF_H1.test.bed")
ChIPseq_CTCF_demo =readPeakFile("CTCF_H1.test.bed",header = FALSE)

## renaming columns, skipping for now
#ChIPseq_CTCF_demo$CTCF_peaks = ChIPseq_CTCF_demo$V4 
#ChIPseq_CTCF_demo$level = ChIPseq_CTCF_demo$V5 
#ChIPseq_CTCF_demo$V4 = ChIPseq_CTCF_demo$V5 = NULL 
#ChIPseq_CTCF_demo

download.file("https://raw.githubusercontent.com/YuLab-SMU/ChIPseeker_current_protocols/master/CP_demo_data/H3K4me1_H1.test.bed", destfile = "./test/H3K4me1_H1.test.bed") 
download.file("https://raw.githubusercontent.com/YuLab-SMU/ChIPseeker_current_protocols/master/CP_demo_data/H3K4me3_H1.test.bed", destfile = "./test/H3K4me3_H1.test.bed")
download.file("https://raw.githubusercontent.com/YuLab-SMU/ChIPseeker_current_protocols/master/CP_demo_data/DHSs_H1.test.bed", destfile = "./test/DHSs_H1.test.bed") 
download.file("https://raw.githubusercontent.com/YuLab-SMU/ChIPseeker_current_protocols/master/CP_demo_data/DNAmeth_H1.test.bed", destfile = "./test/DNAmeth_H1.test.bed")
download.file("https://raw.githubusercontent.com/YuLab-SMU/ChIPseeker_current_protocols/master/CP_demo_data/smRNA_H1.test.bed", destfile = "./test/smRNA_H1.test.bed")

ChIPseq_H3K4me1_demo1 = readPeakFile("H3K4me1_H1.test.bed", header=TRUE)
ChIPseq_H3K4me3_demo2 = readPeakFile("H3K4me3_H1.test.bed", header=TRUE)
DNaseseq_demo = readPeakFile("DHSs_H1.test.bed", header=TRUE)
Methylseq_demo = readPeakFile("DNAmeth_H1.test.bed", header=TRUE)

smRNA_demo = readPeakFile("smRNA_H1.test.bed", header=TRUE)

head(getGEOgenomeVersion(), 5)
getGEOInfo(genome = "anoCar2")

downloadGEObedFiles(genome = "anoCar2", destDir = "./test")
gsm = list("GSM1064688","GSM1064689") 
downloadGSMbedFiles(gsm, destDir = "./test")

#######################################################################
## Basic Protocol 2
####### ANNOTATION OF EPIGENOMIC DATASETS #############################
TxDb_hg19 = TxDb.Hsapiens.UCSC.hg19.knownGene

ChIPseq_CTCF_demo_anno_default = annotatePeak(ChIPseq_CTCF_demo, TxDb=TxDb_hg19)

ChIPseq_CTCF_demo_anno_default

## bed.annot_hg19 -- visually looks comparable

head(as.GRanges(ChIPseq_CTCF_demo_anno_default), 5)

head(as.GRanges(bed.annot_hg19), 5)

write.table(as.data.frame(ChIPseq_CTCF_demo_anno_default),file= "./test/ChIPseq_CTCF_demo_anno_default_taable.csv")

ChIPseq_CTCF_demo_anno_change_priority =annotatePeak(ChIPseq_CTCF_demo, TxDb=TxDb_hg19, genomicAnnotationPriority = c("Exon", "Intron", "5UTR", "3UTR", "Promoter", "Downstream", "Intergenic"))

options(ChIPseeker.ignore_1st_exon = TRUE) 
options(ChIPseeker.ignore_1st_intron = TRUE) 
options(ChIPseeker.ignore_downstream = TRUE) 
options(ChIPseeker.ignore_promoter_subcategory = TRUE) 
ChIPseq_CTCF_demo_anno_with_options = annotatePeak(ChIPseq_CTCF_demo, TxDb=TxDb_hg19)

ChIPseq_CTCF_demo_anno_user_defined = annotatePeak(ChIPseq_CTCF_demo, tssRegion = c(-2000,0), TxDb=TxDb_hg19)

ChIPseq_CTCF_demo_anno_gene_name = annotatePeak(ChIPseq_CTCF_demo, tssRegion=c(-2000,0), TxDb=TxDb_hg19, annoDb="org.Hs.eg.db")

head(as.GRanges(ChIPseq_CTCF_demo_anno_gene_name), 5)

ChIPseq_CTCF_demo_anno_flank_5_kb = annotatePeak(ChIPseq_CTCF_demo, tssRegion=c(-2000,0), TxDb=TxDb_hg19, addFlankGeneInfo=TRUE, flankDistance=5000)

head(as.GRanges(ChIPseq_CTCF_demo_anno_flank_5_kb), 5)

Epi_data_list = GRangesList(CTCF=ChIPseq_CTCF_demo, DHSs=DNaseseq_demo, H3K4me1=ChIPseq_H3K4me1_demo1, m5C=Methylseq_demo, smRNA=smRNA_demo)

peakAnnoList_user_defined = lapply(Epi_data_list, annotatePeak, tssRegion=c(-2000,0), TxDb=TxDb_hg19)

user_defined_GRange =GRanges(seqnames =Rle(c("chr1", "chr10", "chr1", "chr20"),  c(1, 3, 1, 5)),
                             ranges =IRanges(start = 55267513:55267522, end = 55714466:55714475),
                             strand =Rle(strand(c("-","+","*","+","-")), c(1, 2, 1, 4, 2)))

ChIPseq_CTCF_demo_anno_GR =annotatePeak(ChIPseq_CTCF_demo, TxDb = user_defined_GRange)

as.GRanges(ChIPseq_CTCF_demo_anno_GR)[as.GRanges(ChIPseq_CTCF_demo_anno_GR) $distanceToTSS == 0]

CTCF_demo_anno_with_m5C_demo = annotatePeak(Epi_data_list$CTCF, TxDb=Epi_data_list$m5C)


#######################################################################
## Basic Protocol 3
####### COMPARISON OF EPIGENOMIC DATASETS #############################
BiocManager::install("ggVennDiagram") 
library(ggVennDiagram)

vennplot(
  list(DHSs = as.data.frame(peakAnnoList_user_defined$DHSs)$geneId,
       CTCF =as.data.frame(peakAnnoList_user_defined$CTCF)$geneId),
  by = "ggVennDiagram")

peakAnnoList_user_defined_gene = lapply(peakAnnoList_user_defined, function(i) as.data.frame(i)$geneId)

vennplot(peakAnnoList_user_defined_gene, by = "ggVennDiagram",
         label_percent_digit = 2, edge_size = 1.5)

files <- getSampleFiles() 
enrichPeakOverlap(queryPeak = files[[5]], targetPeak =unlist(files[1:4]), TxDb = TxDb_hg19, pAdjustMethod = "BH", nShuffle = 10, chainFile = NULL, verbose = FALSE)

#######################################################################
## Basic Protocol 4
####### VISUALIZATION OF ANNOTATED RESULTS #############################
plotAnnoPie(ChIPseq_CTCF_demo_anno_default)

plotAnnoBar(ChIPseq_CTCF_demo_anno_user_defined)

vennpie(ChIPseq_CTCF_demo_anno_default)

upsetplot(ChIPseq_CTCF_demo_anno_default)

BiocManager::install("ggimage") 
library(ggimage)

upsetplot(ChIPseq_CTCF_demo_anno_default, vennpie=TRUE)




