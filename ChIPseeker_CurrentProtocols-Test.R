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


