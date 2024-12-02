library(ChIPseeker)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(clusterProfiler)
library(GenomicRanges)
library(IRanges)
library(dplyr)

#### D
peak_file <- "/path/to/SRR396786_peaks.narrowPeak"
peaks <- readPeakFile(peak_file)
# Annotate peaks with TxDb
txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
peak_anno <- annotatePeak(peaks, TxDb = txdb)
# Plot distribution of peaks across genomic features
plotAnnoBar(peak_anno, title= 'MyoD1')
upsetplot(peak_anno, vennpie=TRUE)