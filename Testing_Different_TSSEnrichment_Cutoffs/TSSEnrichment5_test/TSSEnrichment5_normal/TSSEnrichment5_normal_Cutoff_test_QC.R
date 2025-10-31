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

unintegrated_ALS_PN <- readRDS("/home/guidubjl/Signac/results_integrated_ALS_PN/Signac_IndividualSample_QC/Testing_Different_TSSEnrichment_Cutoffs/unintegrated_ALS_PN_TSS5_normal.rds")


#Visualizing the relationship between these different variables in the metadata

DensityScatter <- DensityScatter(unintegrated_ALS_PN, x = 'nCount_peaks', y = 'TSS.enrichment', log_x = TRUE, quantiles = TRUE)
ggsave("DensityScatter.png", plot = DensityScatter, width = 12, height = 6, dpi = 300)
#DensityScatter




#Plotting the distribution of each QC metric separately with violin plots

unintegrated_ALS_PN_QC_ViolinPlots <- VlnPlot(
  object = unintegrated_ALS_PN,
  features = c('nCount_peaks', 'TSS.enrichment', 'blacklist_ratio', 'nucleosome_signal', 'pct_reads_in_peaks'),
  pt.size = 0.1,
  ncol = 5
)

ggsave("unintegrated_ALS_PN_QC_ViolinPlots.png", plot = unintegrated_ALS_PN_QC_ViolinPlots, width = 12, height = 6, dpi = 300)

#unintegrated_ALS_PN_QC_ViolinPlots
