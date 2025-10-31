
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



#Visualizing the relationship between these different variables in the metadata

DensityScatter <- DensityScatter(unintegrated_ALS_PN, x = 'nCount_peaks', y = 'TSS.enrichment', log_x = TRUE, quantiles = TRUE)
ggsave("DensityScatter.png", plot = DensityScatter, width = 12, height = 6, dpi = 300)
DensityScatter



#Filtering out cells that are outliers for the mononucleosomal/nucleosome-free ratio

#unintegrated_ALS_PN$nucleosome_group <- ifelse(unintegrated_ALS_PN$nucleosome_signal > 3, 'NS > 3', 'NS < 3')
#FragmentHistogram(object = unintegrated_ALS_PN, group.by = 'nucleosome_group')

#Note the y-axes on these plots
#This plot looks very different from the tutorial



#Plotting the distribution of each QC metric separately with violin plots

#unintegrated_ALS_PN_QC_ViolinPlots <- VlnPlot(
 # object = unintegrated_ALS_PN,
 #features = c('nCount_peaks', 'TSS.enrichment', 'blacklist_ratio', 'nucleosome_signal', 'pct_reads_in_peaks'),
 # pt.size = 0.1,
 # ncol = 5
#)

#ggsave("unintegrated_ALS_PN_QC_ViolinPlots.png", plot = unintegrated_ALS_PN_QC_ViolinPlots, width = 12, height = 6, dpi = 300)

#unintegrated_ALS_PN_QC_ViolinPlots




# Alternative version with better formatting for many samples
#TSS_enrichment_by_sample <- VlnPlot(
#  object = unintegrated_ALS_PN,
#  features = "TSS.enrichment",
#  group.by = "orig.ident",
#  pt.size = 0.05,  # Smaller points since you have many samples
#  cols = rainbow(length(unique(unintegrated_ALS_PN$orig.ident)))  # Different colors for each sample
#) + 
#  theme(
#    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 8),
#    legend.position = "none"  # Remove legend since x-axis shows sample names
#  ) +
#  labs(x = "Sample", y = "TSS Enrichment", title = "TSS Enrichment by Sample")


TSS_enrichment_by_sample <- VlnPlot(
  object = unintegrated_ALS_PN,
  features = 'TSS.enrichment',
  split.by = "dataset",
  pt.size = 0,
  same.y.lims = TRUE
)


TSS_enrichment_by_sample

ggsave("TSS_enrichment_by_sample_NoDot.png", plot = TSS_enrichment_by_sample, width = 12, height = 6, dpi = 300)
