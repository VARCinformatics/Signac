#Start by loading in the necessary libraries
library(Signac)
library(Seurat)
library(GenomicRanges)
library(ggplot2)
library(patchwork)

# install.packages("harmony")
# install.packages("hdf5r")
# BiocManager::install("Rsamtools")

library(harmony)
library(hdf5r)
library(Rsamtools)
#Computing QC Metrics

unintegrated_ALS_PN <- readRDS("/home/guidubjl/Signac/results_integrated_ALS_PN/Signac_w_Harmony_Complete_Run/unintegrated_ALS_QC_PN.rds")



# TSS > 2
#unintegrated_ALS_PN_TSS2 <- subset(
#  x = unintegrated_ALS_PN,
#  subset = nCount_peaks > 9000 &
#    nCount_peaks < 100000 &
#    pct_reads_in_peaks > 40 &
#    blacklist_ratio < 0.01 &
#    nucleosome_signal < 3 &
#    TSS.enrichment > 2
#)

# TSS > 4  
#unintegrated_ALS_PN_TSS4 <- subset(
#  x = unintegrated_ALS_PN,
#  subset = nCount_peaks > 9000 &
#    nCount_peaks < 100000 &
#    pct_reads_in_peaks > 40 &
#    blacklist_ratio < 0.01 &
#    nucleosome_signal < 3 &
#    TSS.enrichment > 4
#)

# TSS > 6
#unintegrated_ALS_PN_TSS6 <- subset(
#  x = unintegrated_ALS_PN,
#  subset = nCount_peaks > 9000 &
#    nCount_peaks < 100000 &
#    pct_reads_in_peaks > 40 &
#    blacklist_ratio < 0.01 &
#    nucleosome_signal < 3 &
#    TSS.enrichment > 6
#)


# TSS > 5
unintegrated_ALS_PN_TSS5_normal <- subset(
  x = unintegrated_ALS_PN,
  subset = nCount_peaks > 9000 &
    nCount_peaks < 100000 &
    pct_reads_in_peaks > 40 &
    blacklist_ratio < 0.01 &
    nucleosome_signal < 3 &
    TSS.enrichment > 5
)

# TSS > 5 and pct_reads_in_peaks > 60
#IMPORTANT NOTE: I AM ALSO INCREASING THE PCT_READS_IN_PEAKS CUTOFF
unintegrated_ALS_PN_TSS5_MoreStrict <- subset(
  x = unintegrated_ALS_PN,
  subset = nCount_peaks > 9000 &
    nCount_peaks < 100000 &
    pct_reads_in_peaks > 60 &
    blacklist_ratio < 0.01 &
    nucleosome_signal < 3 &
    TSS.enrichment > 5
)


# Print cell counts for comparison
#cat("Cells with TSS > 2:", ncol(unintegrated_ALS_PN_TSS2), "\n")
#cat("Cells with TSS > 4:", ncol(unintegrated_ALS_PN_TSS4), "\n")
#cat("Cells with TSS > 6:", ncol(unintegrated_ALS_PN_TSS6), "\n")

#saveRDS(unintegrated_ALS_PN_TSS2, "unintegrated_ALS_PN_TSS2.rds")
#saveRDS(unintegrated_ALS_PN_TSS4, "unintegrated_ALS_PN_TSS4.rds")
#saveRDS(unintegrated_ALS_PN_TSS6, "unintegrated_ALS_PN_TSS6.rds")



cat("Cells with TSS > 5:", ncol(unintegrated_ALS_PN_TSS5_normal), "\n")
cat("Cells with TSS > 5 and pct_reads_in_peaks > 60:", ncol(unintegrated_ALS_PN_TSS5_MoreStrict), "\n")

saveRDS(unintegrated_ALS_PN_TSS5_normal, "unintegrated_ALS_PN_TSS5_normal.rds")
saveRDS(unintegrated_ALS_PN_TSS5_MoreStrict, "unintegrated_ALS_PN_TSS5_MoreStrict.rds")

