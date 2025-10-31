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

#unintegrated_ALS_PN <- readRDS("/home/guidubjl/Signac/results_integrated_ALS_PN/Signac_w_Harmony_Complete_Run/unintegrated_ALS_PN.rds")


#Explanation for each of these steps in the vignette
#Each one of these adds their values into the Seurat object

# compute nucleosome signal score per cell
#unintegrated_ALS_PN <- NucleosomeSignal(object = unintegrated_ALS_PN)

# compute TSS enrichment score per cell
#unintegrated_ALS_PN <- TSSEnrichment(object = unintegrated_ALS_PN)

# add fraction of reads in peaks
#unintegrated_ALS_PN$pct_reads_in_peaks <- unintegrated_ALS_PN$peak_region_fragments / unintegrated_ALS_PN$passed_filters * 100

# add blacklist ratio
#unintegrated_ALS_PN$blacklist_ratio <- FractionCountsInRegion(
#  object = unintegrated_ALS_PN, 
#  assay = 'peaks',
#  regions = blacklist_hg38_unified
#)



#Saving merged sample with QC
#saveRDS(unintegrated_ALS_PN, file = "unintegrated_ALS_QC_PN.rds")



#Visualizing the relationship between these different variables in the metadata

#DensityScatter <- DensityScatter(unintegrated_ALS_PN, x = 'nCount_peaks', y = 'TSS.enrichment', log_x = TRUE, quantiles = TRUE)
#ggsave("DensityScatter.png", plot = DensityScatter, width = 12, height = 6, dpi = 300)
#DensityScatter



#Filtering out cells that are outliers for the mononucleosomal/nucleosome-free ratio

#unintegrated_ALS_PN$nucleosome_group <- ifelse(unintegrated_ALS_PN$nucleosome_signal > 3, 'NS > 3', 'NS < 3')
#FragmentHistogram(object = unintegrated_ALS_PN, group.by = 'nucleosome_group')

#Note the y-axes on these plots
#This plot looks very different from the tutorial



#Plotting the distribution of each QC metric separately with violin plots

#unintegrated_ALS_PN_QC_ViolinPlots <- VlnPlot(
#  object = unintegrated_ALS_PN,
#  features = c('nCount_peaks', 'TSS.enrichment', 'blacklist_ratio', 'nucleosome_signal', 'pct_reads_in_peaks'),
#  pt.size = 0.1,
#  ncol = 5
#)

#ggsave("unintegrated_ALS_PN_QC_ViolinPlots.png", plot = unintegrated_ALS_PN_QC_ViolinPlots, width = 12, height = 6, dpi = 300)

#unintegrated_ALS_PN_QC_ViolinPlots

#The resolution on these is terrible in this viewer, but if you save the image as a .png they display much better.

#This plot also looks very different from the tutorial data.



#Remove cells based on outliers for these QC metrics.
#The exact thresholds for these can and should all be tailored based on your dataset

#unintegrated_ALS_PN <- subset(
#  x = unintegrated_ALS_PN,
#  subset = nCount_peaks > 9000 &
#    nCount_peaks < 100000 &
#    pct_reads_in_peaks > 40 &
#    blacklist_ratio < 0.01 &
#    nucleosome_signal < 3 &
#    TSS.enrichment > 3
#)
#unintegrated_ALS_PN

#Saving filted object
#saveRDS(unintegrated_ALS_PN, file = "unintegrated_ALS_PN_QC_Filtered.rds")



#Normalization and Linear Dimensional Reduction
#Details for both of these can be found in the vignette

#NOTE: One thing that's important to consider is the cutoff threshold for features.
#In this example, they are NOT excluding any features!
#However, in your own samples, you may want to set this threshold higher (ex. 'q75' to use only the top 25% of peaks), as it will run faster.

#unintegrated_ALS_PN <- RunTFIDF(unintegrated_ALS_PN) #This is the normalization step: Term frequency-inverse document frequency normalization
#unintegrated_ALS_PN <- FindTopFeatures(unintegrated_ALS_PN, min.cutoff = 'q0') #This is the feature selection step: min.cutoff sets the percentile threshold
#unintegrated_ALS_PN <- RunSVD(unintegrated_ALS_PN) #This is the dimension reduction step: singular value decomposition

#The combination of TF-IDF and SVD is known as "latent semantic indexing" (LSI). After running these steps, you should see that appear in the Seurat obj
#unintegrated_ALS_PN

#The output from all of these steps is *conceptually* similar to what you get after PCA in scRNAseq analysis

#Saving object
#saveRDS(unintegrated_ALS_PN, file = "unintegrated_ALS_PN_LinDimRed.rds")



#The first LSI component frequently captures sequencing depth (technical variation) rather that biological variation.
#If this is the case, it should be removed from downstream analysis.
#We can assess the correlation between each LSI component and sequencing depth with the following function:

#DepthCor <- DepthCor(unintegrated_ALS_PN)
#ggsave("DepthCor.png", plot = DepthCor, width = 12, height = 6, dpi = 300)
#DepthCor



#Looking at the above plot for the tutorial data, it is clear that the first LSI component is basically just capturing sequencing depth,
#thus it should not be used in downstream analyses.



#Non-linear dimension reduction and clustering

#These steps perform nonlinear and graph-based clustering on the now-relatively-low-dimensional data.
#These steps all come from the Seurat pacakge and are the same as what are used in scRNAseq analysis

#unintegrated_ALS_PN <- RunUMAP(object = unintegrated_ALS_PN, reduction = 'lsi', dims = 2:30)
#unintegrated_ALS_PN <- FindNeighbors(object = unintegrated_ALS_PN, reduction = 'lsi', dims = 2:30)
#unintegrated_ALS_PN <- FindClusters(object = unintegrated_ALS_PN, verbose = FALSE, algorithm = 3)
#unintegrated_ALS_PN_clusters <- DimPlot(object = unintegrated_ALS_PN, label = TRUE) + NoLegend()
#unintegrated_ALS_PN_overlap <- DimPlot(unintegrated_ALS_PN, group.by = 'dataset', pt.size = 0.1)

#ggsave("unintegrated_ALS_PN_clusters.png", plot = unintegrated_ALS_PN_clusters, width = 12, height = 6, dpi = 300)
#ggsave("unintegrated_ALS_PN_overlap.png", plot = unintegrated_ALS_PN_overlap, width = 12, height = 6, dpi = 300)

#unintegrated_ALS_PN_clusters
#unintegrated_ALS_PN_overlap

#Saving rds
#saveRDS(unintegrated_ALS_PN, file = "unintegrated_ALS_PN_NonLimDimRed.rds")


#Nelow step only required if starting from an intermediary .rds
#unintegrated_ALS_PN <- readRDS("/home/guidubjl/Signac/results_integrated_ALS_PN/unintegrated_ALS_PN_NonLimDimRed.rds")


#Batch correct using Harmony
#hm_integrated_ALS_PN <- RunHarmony(object = unintegrated_ALS_PN, group.by.vars = 'dataset', reduction.use = 'lsi', assay.use = 'peaks', project.dim = FALSE)
#hm_integrated_ALS_PN <- RunUMAP(hm_integrated_ALS_PN, dims = 2:30, reduction = 'harmony')
#hm_integrated_ALS_PN_UMAP_Plot <- DimPlot(hm_integrated_ALS_PN, group.by = 'dataset', pt.size = 0.1)

#rm(unintegrated_ALS_PN)

#ggsave("hm_integrated_ALS_PN_UMAP_Plot.png", plot = hm_integrated_ALS_PN_UMAP_Plot, width = 12, height = 6, dpi = 300)

#hm_integrated_ALS_PN_UMAP_Plot



#sessionInfo <- sessionInfo()
#saveRDS(sessionInfo, file = "sessionInfo.rds")
#saveRDS(hm_integrated_ALS_PN, file = "hm_integrated_ALS_PN.rds")

#Below step only required if starting from an intermediary .rds
#hm_integrated_ALS_PN <- readRDS("/home/guidubjl/Signac/results_integrated_ALS_PN/Signac_w_Harmony_Complete_Run/hm_integrated_ALS_PN.rds")


#gene.activities <- GeneActivity(hm_integrated_ALS_PN)


# add the gene activity matrix to the Seurat object as a new assay and normalize it
#hm_integrated_ALS_PN[['RNA']] <- CreateAssayObject(counts = gene.activities)
#hm_integrated_ALS_PN <- NormalizeData(
#  object = hm_integrated_ALS_PN,
#  assay = 'RNA',
#  normalization.method = 'LogNormalize',
#  scale.factor = median(hm_integrated_ALS_PN$nCount_RNA)
#)


#DefaultAssay(hm_integrated_ALS_PN) <- 'RNA'

#glymph_genes_clusters <- FeaturePlot(
#  object = hm_integrated_ALS_PN,
#  features = c('GFAP', 'AQP4', 'VEGFA', 'VEGFB', 'SHTN1'),
#  pt.size = 0.1,
#  max.cutoff = 'q95',
#  ncol = 3
#)

#ggsave("glymph_genes_clusters.png", plot = glymph_genes_clusters, width = 12, height = 6, dpi = 300)
#glymph_genes_clusters

#saveRDS(hm_integrated_ALS_PN, file = "hm_integrated_ALS_PN_GeneActivity.rds")

hm_integrated_ALS_PN <- readRDS("/home/guidubjl/Signac/results_integrated_ALS_PN/Signac_w_Harmony_Complete_Run/hm_integrated_ALS_PN_GeneActivity.rds")

## Integrating with scRNA-seq data


# Load the pre-processed scRNA-seq data
library(dplyr)
library(purrr)
library(SummarizedExperiment)


#ALS_rna <- readRDS("/data/l3_varc/Data/snATAC+RNAseq/2021_MIT_D002_snRNA-seq_PFC-MCX/se_agg_BA4_fil/se_agg_BA4_SALS_fil.rds")
#PN_rna <- readRDS("/data/l3_varc/Data/snATAC+RNAseq/2021_MIT_D002_snRNA-seq_PFC-MCX/se_agg_BA4_fil/se_agg_BA4_PN_fil.rds")
#merged_rna_data <- cbind(ALS_rna, PN_rna)

#merged_rna_data_counts <- SummarizedExperiment::assay(merged_rna_data, "counts")
#rownames(merged_rna_data_counts) <- rowRanges(merged_rna_data)@elementMetadata$Gene
#merged_rna_data_counts <- merged_rna_data_counts[rownames(merged_rna_data) != "N/A", ]
#merged_rna_data_counts <- merged_rna_data_counts[match(unique(rownames(merged_rna_data_counts)), rownames(merged_rna_data_counts)),]
#merged_rna_data_metadata <- as.data.frame(SummarizedExperiment::colData(merged_rna_data))

#merged_rna_data_seurat <- CreateSeuratObject(counts = merged_rna_data_counts, meta.data = merged_rna_data_metadata, assay = "RNA")

#merged_rna_data_seurat_updated <- UpdateSeuratObject(merged_rna_data_seurat)

#Adding in a portion here where I run my RNAseq data through Seurat
library(Seurat)
packageVersion('Seurat')

#Installing recommended performance-enhancing packages
#setRepositories(ind = 1:3, addURLs = c('https://satijalab.r-universe.dev', 'https://bnprks.r-universe.dev/'))
#install.packages(c("BPCells", "presto", "glmGamPoi"))
library(BPCells)
library(presto)
library(glmGamPoi)
packageVersion('BPCells')
packageVersion('presto')
packageVersion('glmGamPoi')

#It creates a bunch of issues if you load these modules later, so instead try loading them first and then
#different function calls can hopefully get overwritten later.
#If that still doesn't work, then just don't ever load them.
###I actually think that they might not be necessary to ever load, so don't unless you need to!
packageVersion('Signac')
# #devtools::install_github("satijalab/seurat-data")
#packageVersion('SeuratData')
#library(SeuratData)
# #remotes::install_github("satijalab/azimuth", ref = 'master')
#packageVersion('Azimuth')
#library(Azimuth)

library(ggplot2)

#The below two steps are only if you are starting with the second half of the code vs starting at the beginning
#From this point forward, it should only take about 15-20min to run using 64GB RAM.
#hm_integrated_ALS_PN <- readRDS("/home/guidubjl/Signac/results_integrated_ALS_PN/hm_integrated_ALS_PN_GeneActivity.rds")
#merged_rna_data_complete <- readRDS("/home/guidubjl/Signac/results_integrated_ALS_PN/merged_rna_data_complete.rds")


#There's an error that comes up in the next PercentageFeatureSet command. It is supposedly fixed with this:
#You need to have ACCRE's BioC think loaded in order for the following library calls to work
#system("ml r-bundle-bioconductor/3.20")
#library("TFBSTools")
#library("irlba")

#merged_rna_data_complete <- merged_rna_data_seurat_updated

#saveRDS(merged_rna_data_complete, file = "merged_rna_data_complete.rds")


merged_rna_data_complete <- readRDS("/home/guidubjl/Signac/results_integrated_ALS_PN/Signac_w_Harmony_Complete_Run/merged_rna_data_complete.rds")

#rm(merged_rna_data)
#rm(merged_rna_data_seurat)
#rm(merged_rna_data_seurat_updated)

# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
merged_rna_data_complete[["percent.mt"]] <- PercentageFeatureSet(merged_rna_data_complete, pattern = "^MT-")

#Subset the data based on QC cutoffs
merged_rna_data_complete <- subset(merged_rna_data_complete, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

#Normalize the data
merged_rna_data_complete <- NormalizeData(merged_rna_data_complete, normalization.method = "LogNormalize", scale.factor = 10000)


#Identification of highly variable features (feature selection)

#Find the top 2000 most variable features
merged_rna_data_complete <- FindVariableFeatures(merged_rna_data_complete, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(merged_rna_data_complete), 10)


#Scale the data
all.genes <- rownames(merged_rna_data_complete)
merged_rna_data_complete <- ScaleData(merged_rna_data_complete, features = all.genes)

#Perform linear dimensional reduction via PCA
#N.B. By default, this pulls in previously determined variable features.
#It is possible to instead run a custom list if you want - just make sure you've ran it through ScaleData() first.

merged_rna_data_complete <- RunPCA(merged_rna_data_complete, features = VariableFeatures(object = merged_rna_data_complete))




transfer.anchors <- FindTransferAnchors(
  reference = merged_rna_data_complete,
  query = hm_integrated_ALS_PN,
  reduction = 'cca'
)

predicted.labels <- TransferData(
  anchorset = transfer.anchors,
  refdata = merged_rna_data_complete$CellType,
  weight.reduction = hm_integrated_ALS_PN[['lsi']],
  dims = 2:30
)

hm_integrated_ALS_PN <- AddMetaData(object = hm_integrated_ALS_PN, metadata = predicted.labels)


plot1 <- DimPlot(
  object = merged_rna_data_complete,
  group.by = 'CellType',
  label = TRUE,
  repel = TRUE) + NoLegend() + ggtitle('scRNA-seq')

plot2 <- DimPlot(
  object = hm_integrated_ALS_PN,
  group.by = 'predicted.id',
  label = TRUE,
  repel = TRUE) + NoLegend() + ggtitle('scATAC-seq')

plot1 + plot2

combined_cluster_plots <- plot1	+ plot2

ggsave("snRNAseq_clusters.png", plot = plot1, width = 12, height = 6, dpi = 300)
ggsave("snATACseq_clusters.png", plot = plot2, width = 12, height = 6, dpi = 300)
ggsave("combined_cluster_plots.png", plot = combined_cluster_plots, width = 12, height = 6, dpi = 300)

VlnPlot(hm_integrated_ALS_PN, 'prediction.score.max', group.by = 'predicted.id')

pred_score_max_plot <- VlnPlot(hm_integrated_ALS_PN, 'prediction.score.max', group.by = 'predicted.id')
ggsave("pred_score_max_plot.png", plot = pred_score_max_plot, width = 12, height = 6, dpi = 300)
