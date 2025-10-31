#Signac + Harmony for Accre

#This is a new cumulative script for all Signac + Harmony ALS+PN Samples
#Stitching together all the working ACCRE scripts

#Data I will be using:
#ATAC data in ACCRE: /data/l3_varc/Data/snATAC+RNAseq/2021_MIT_D001_snATAC-seq_PFC-MCX_ALS-FTD-PN/191112_191114_ALS_FTLD_PN_MCX_snATAC/cellranger_out/
#RNAseq data in ACCRE: /data/l3_varc/Data/snATAC+RNAseq/2021_MIT_D002_snRNA-seq_PFC-MCX/se_agg_BA4_fil
#On local device: /Users/guidubjl/Desktop/Belzil_Lab/snRNAseq_snATACseq_Project/Signac/Data/



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



#Load in the data

#/data/l3_varc/Data/snATAC+RNAseq/2021_MIT_D001_snATAC-seq_PFC-MCX_ALS-FTD-PN/191112_191114_ALS_FTLD_PN_MCX_snATAC/cellranger_out/

PN301_counts <- Read10X_h5(filename = "/data/l3_varc/Data/snATAC+RNAseq/2021_MIT_D001_snATAC-seq_PFC-MCX_ALS-FTD-PN/191112_191114_ALS_FTLD_PN_MCX_snATAC/cellranger_out/191114_PN_301_snATAC-C3/filtered_peak_bc_matrix.h5")

PN303_counts <- Read10X_h5(filename = "/data/l3_varc/Data/snATAC+RNAseq/2021_MIT_D001_snATAC-seq_PFC-MCX_ALS-FTD-PN/191112_191114_ALS_FTLD_PN_MCX_snATAC/cellranger_out/191114_PN_303_snATAC-C4/filtered_peak_bc_matrix.h5")

PN304_counts <- Read10X_h5(filename = "/data/l3_varc/Data/snATAC+RNAseq/2021_MIT_D001_snATAC-seq_PFC-MCX_ALS-FTD-PN/191112_191114_ALS_FTLD_PN_MCX_snATAC/cellranger_out/191114_PN_304_snATAC-C5/filtered_peak_bc_matrix.h5")

PN306_counts <- Read10X_h5(filename = "/data/l3_varc/Data/snATAC+RNAseq/2021_MIT_D001_snATAC-seq_PFC-MCX_ALS-FTD-PN/191112_191114_ALS_FTLD_PN_MCX_snATAC/cellranger_out/191114_PN_306_snATAC-C6/filtered_peak_bc_matrix.h5")

PN307_counts <- Read10X_h5(filename = "/data/l3_varc/Data/snATAC+RNAseq/2021_MIT_D001_snATAC-seq_PFC-MCX_ALS-FTD-PN/191112_191114_ALS_FTLD_PN_MCX_snATAC/cellranger_out/191114_PN_307_snATAC-C7/filtered_peak_bc_matrix.h5")

PN308_counts <- Read10X_h5(filename = "/data/l3_varc/Data/snATAC+RNAseq/2021_MIT_D001_snATAC-seq_PFC-MCX_ALS-FTD-PN/191112_191114_ALS_FTLD_PN_MCX_snATAC/cellranger_out/191114_PN_308_snATAC-C8/filtered_peak_bc_matrix.h5")

PN309_counts <- Read10X_h5(filename = "/data/l3_varc/Data/snATAC+RNAseq/2021_MIT_D001_snATAC-seq_PFC-MCX_ALS-FTD-PN/191112_191114_ALS_FTLD_PN_MCX_snATAC/cellranger_out/191114_PN_309_snATAC-C9/filtered_peak_bc_matrix.h5")

PN317_counts <- Read10X_h5(filename = "/data/l3_varc/Data/snATAC+RNAseq/2021_MIT_D001_snATAC-seq_PFC-MCX_ALS-FTD-PN/191112_191114_ALS_FTLD_PN_MCX_snATAC/cellranger_out/191114_PN_317_snATAC-C10/filtered_peak_bc_matrix.h5")

PN318_counts <- Read10X_h5(filename = "/data/l3_varc/Data/snATAC+RNAseq/2021_MIT_D001_snATAC-seq_PFC-MCX_ALS-FTD-PN/191112_191114_ALS_FTLD_PN_MCX_snATAC/cellranger_out/191114_PN_318_snATAC-C11/filtered_peak_bc_matrix.h5")

PN319_counts <- Read10X_h5(filename = "/data/l3_varc/Data/snATAC+RNAseq/2021_MIT_D001_snATAC-seq_PFC-MCX_ALS-FTD-PN/191112_191114_ALS_FTLD_PN_MCX_snATAC/cellranger_out/191114_PN_319_snATAC-C12/filtered_peak_bc_matrix.h5")

PN322_counts <- Read10X_h5(filename = "/data/l3_varc/Data/snATAC+RNAseq/2021_MIT_D001_snATAC-seq_PFC-MCX_ALS-FTD-PN/191112_191114_ALS_FTLD_PN_MCX_snATAC/cellranger_out/191114_PN_322_snATAC-D1/filtered_peak_bc_matrix.h5")

PN323_counts <- Read10X_h5(filename = "/data/l3_varc/Data/snATAC+RNAseq/2021_MIT_D001_snATAC-seq_PFC-MCX_ALS-FTD-PN/191112_191114_ALS_FTLD_PN_MCX_snATAC/cellranger_out/191114_PN_323_snATAC-D2/filtered_peak_bc_matrix.h5")

PN324_counts <- Read10X_h5(filename = "/data/l3_varc/Data/snATAC+RNAseq/2021_MIT_D001_snATAC-seq_PFC-MCX_ALS-FTD-PN/191112_191114_ALS_FTLD_PN_MCX_snATAC/cellranger_out/191114_PN_324_snATAC-D3/filtered_peak_bc_matrix.h5")

PN325_counts <- Read10X_h5(filename = "/data/l3_varc/Data/snATAC+RNAseq/2021_MIT_D001_snATAC-seq_PFC-MCX_ALS-FTD-PN/191112_191114_ALS_FTLD_PN_MCX_snATAC/cellranger_out/191114_PN_325_snATAC-D4/filtered_peak_bc_matrix.h5")

PN328_counts <- Read10X_h5(filename = "/data/l3_varc/Data/snATAC+RNAseq/2021_MIT_D001_snATAC-seq_PFC-MCX_ALS-FTD-PN/191112_191114_ALS_FTLD_PN_MCX_snATAC/cellranger_out/191114_PN_328_snATAC-D5/filtered_peak_bc_matrix.h5")



ALS110_counts <- Read10X_h5(filename = "/data/l3_varc/Data/snATAC+RNAseq/2021_MIT_D001_snATAC-seq_PFC-MCX_ALS-FTD-PN/191112_191114_ALS_FTLD_PN_MCX_snATAC/cellranger_out/191112_ALS_110_snATAC-H7/filtered_peak_bc_matrix.h5")

ALS111_counts <- Read10X_h5(filename = "/data/l3_varc/Data/snATAC+RNAseq/2021_MIT_D001_snATAC-seq_PFC-MCX_ALS-FTD-PN/191112_191114_ALS_FTLD_PN_MCX_snATAC/cellranger_out/191112_ALS_111_snATAC-H8/filtered_peak_bc_matrix.h5")

ALS112_counts <- Read10X_h5(filename = "/data/l3_varc/Data/snATAC+RNAseq/2021_MIT_D001_snATAC-seq_PFC-MCX_ALS-FTD-PN/191112_191114_ALS_FTLD_PN_MCX_snATAC/cellranger_out/191112_ALS_112_snATAC-H9/filtered_peak_bc_matrix.h5")

ALS114_counts <- Read10X_h5(filename = "/data/l3_varc/Data/snATAC+RNAseq/2021_MIT_D001_snATAC-seq_PFC-MCX_ALS-FTD-PN/191112_191114_ALS_FTLD_PN_MCX_snATAC/cellranger_out/191112_ALS_114_snATAC-H10/filtered_peak_bc_matrix.h5")

ALS115_counts <- Read10X_h5(filename = "/data/l3_varc/Data/snATAC+RNAseq/2021_MIT_D001_snATAC-seq_PFC-MCX_ALS-FTD-PN/191112_191114_ALS_FTLD_PN_MCX_snATAC/cellranger_out/191112_ALS_115_snATAC-H11/filtered_peak_bc_matrix.h5")

ALS116_counts <- Read10X_h5(filename = "/data/l3_varc/Data/snATAC+RNAseq/2021_MIT_D001_snATAC-seq_PFC-MCX_ALS-FTD-PN/191112_191114_ALS_FTLD_PN_MCX_snATAC/cellranger_out/191112_ALS_116_snATAC-H12/filtered_peak_bc_matrix.h5")

ALS117_counts <- Read10X_h5(filename = "/data/l3_varc/Data/snATAC+RNAseq/2021_MIT_D001_snATAC-seq_PFC-MCX_ALS-FTD-PN/191112_191114_ALS_FTLD_PN_MCX_snATAC/cellranger_out/191112_ALS_117_snATAC-A1/filtered_peak_bc_matrix.h5")

ALS118_counts <- Read10X_h5(filename = "/data/l3_varc/Data/snATAC+RNAseq/2021_MIT_D001_snATAC-seq_PFC-MCX_ALS-FTD-PN/191112_191114_ALS_FTLD_PN_MCX_snATAC/cellranger_out/191112_ALS_118_snATAC-A2/filtered_peak_bc_matrix.h5")

ALS120_counts <- Read10X_h5(filename = "/data/l3_varc/Data/snATAC+RNAseq/2021_MIT_D001_snATAC-seq_PFC-MCX_ALS-FTD-PN/191112_191114_ALS_FTLD_PN_MCX_snATAC/cellranger_out/191112_ALS_120_snATAC-A3/filtered_peak_bc_matrix.h5")

ALS122_counts <- Read10X_h5(filename = "/data/l3_varc/Data/snATAC+RNAseq/2021_MIT_D001_snATAC-seq_PFC-MCX_ALS-FTD-PN/191112_191114_ALS_FTLD_PN_MCX_snATAC/cellranger_out/191112_ALS_122_snATAC-A4/filtered_peak_bc_matrix.h5")

ALS124_counts <- Read10X_h5(filename = "/data/l3_varc/Data/snATAC+RNAseq/2021_MIT_D001_snATAC-seq_PFC-MCX_ALS-FTD-PN/191112_191114_ALS_FTLD_PN_MCX_snATAC/cellranger_out/191112_ALS_124_snATAC-A5/filtered_peak_bc_matrix.h5")

ALS126_counts <- Read10X_h5(filename = "/data/l3_varc/Data/snATAC+RNAseq/2021_MIT_D001_snATAC-seq_PFC-MCX_ALS-FTD-PN/191112_191114_ALS_FTLD_PN_MCX_snATAC/cellranger_out/191112_ALS_126_snATAC-A6/filtered_peak_bc_matrix.h5")

ALS127_counts <- Read10X_h5(filename = "/data/l3_varc/Data/snATAC+RNAseq/2021_MIT_D001_snATAC-seq_PFC-MCX_ALS-FTD-PN/191112_191114_ALS_FTLD_PN_MCX_snATAC/cellranger_out/191112_ALS_127_snATAC-A7/filtered_peak_bc_matrix.h5")

ALS128_counts <- Read10X_h5(filename = "/data/l3_varc/Data/snATAC+RNAseq/2021_MIT_D001_snATAC-seq_PFC-MCX_ALS-FTD-PN/191112_191114_ALS_FTLD_PN_MCX_snATAC/cellranger_out/191112_ALS_128_snATAC-A8/filtered_peak_bc_matrix.h5")

ALS129_counts <- Read10X_h5(filename = "/data/l3_varc/Data/snATAC+RNAseq/2021_MIT_D001_snATAC-seq_PFC-MCX_ALS-FTD-PN/191112_191114_ALS_FTLD_PN_MCX_snATAC/cellranger_out/191112_ALS_129_snATAC-A9/filtered_peak_bc_matrix.h5")

ALS130_counts <- Read10X_h5(filename = "/data/l3_varc/Data/snATAC+RNAseq/2021_MIT_D001_snATAC-seq_PFC-MCX_ALS-FTD-PN/191112_191114_ALS_FTLD_PN_MCX_snATAC/cellranger_out/191112_ALS_130_snATAC-A10/filtered_peak_bc_matrix.h5")




PN301_metadata <- read.csv(
  file = "/data/l3_varc/Data/snATAC+RNAseq/2021_MIT_D001_snATAC-seq_PFC-MCX_ALS-FTD-PN/191112_191114_ALS_FTLD_PN_MCX_snATAC/cellranger_out/191114_PN_301_snATAC-C3/singlecell.csv",
  header = TRUE,
  row.names = 1
)

PN303_metadata <- read.csv(
  file = "/data/l3_varc/Data/snATAC+RNAseq/2021_MIT_D001_snATAC-seq_PFC-MCX_ALS-FTD-PN/191112_191114_ALS_FTLD_PN_MCX_snATAC/cellranger_out/191114_PN_303_snATAC-C4/singlecell.csv",
  header = TRUE,
  row.names = 1
)

PN304_metadata <- read.csv(
  file = "/data/l3_varc/Data/snATAC+RNAseq/2021_MIT_D001_snATAC-seq_PFC-MCX_ALS-FTD-PN/191112_191114_ALS_FTLD_PN_MCX_snATAC/cellranger_out/191114_PN_304_snATAC-C5/singlecell.csv",
  header = TRUE,
  row.names = 1
)

PN306_metadata <- read.csv(
  file = "/data/l3_varc/Data/snATAC+RNAseq/2021_MIT_D001_snATAC-seq_PFC-MCX_ALS-FTD-PN/191112_191114_ALS_FTLD_PN_MCX_snATAC/cellranger_out/191114_PN_306_snATAC-C6/singlecell.csv",
  header = TRUE,
  row.names = 1
)

PN307_metadata <- read.csv(
  file = "/data/l3_varc/Data/snATAC+RNAseq/2021_MIT_D001_snATAC-seq_PFC-MCX_ALS-FTD-PN/191112_191114_ALS_FTLD_PN_MCX_snATAC/cellranger_out/191114_PN_307_snATAC-C7/singlecell.csv",
  header = TRUE,
  row.names = 1
)

PN308_metadata <- read.csv(
  file = "/data/l3_varc/Data/snATAC+RNAseq/2021_MIT_D001_snATAC-seq_PFC-MCX_ALS-FTD-PN/191112_191114_ALS_FTLD_PN_MCX_snATAC/cellranger_out/191114_PN_308_snATAC-C8/singlecell.csv",
  header = TRUE,
  row.names = 1
)

PN309_metadata <- read.csv(
  file = "/data/l3_varc/Data/snATAC+RNAseq/2021_MIT_D001_snATAC-seq_PFC-MCX_ALS-FTD-PN/191112_191114_ALS_FTLD_PN_MCX_snATAC/cellranger_out/191114_PN_309_snATAC-C9/singlecell.csv",
  header = TRUE,
  row.names = 1
)

PN317_metadata <- read.csv(
  file = "/data/l3_varc/Data/snATAC+RNAseq/2021_MIT_D001_snATAC-seq_PFC-MCX_ALS-FTD-PN/191112_191114_ALS_FTLD_PN_MCX_snATAC/cellranger_out/191114_PN_317_snATAC-C10/singlecell.csv",
  header = TRUE,
  row.names = 1
)

PN318_metadata <- read.csv(
  file = "/data/l3_varc/Data/snATAC+RNAseq/2021_MIT_D001_snATAC-seq_PFC-MCX_ALS-FTD-PN/191112_191114_ALS_FTLD_PN_MCX_snATAC/cellranger_out/191114_PN_318_snATAC-C11/singlecell.csv",
  header = TRUE,
  row.names = 1
)

PN319_metadata <- read.csv(
  file = "/data/l3_varc/Data/snATAC+RNAseq/2021_MIT_D001_snATAC-seq_PFC-MCX_ALS-FTD-PN/191112_191114_ALS_FTLD_PN_MCX_snATAC/cellranger_out/191114_PN_319_snATAC-C12/singlecell.csv",
  header = TRUE,
  row.names = 1
)

PN322_metadata <- read.csv(
  file = "/data/l3_varc/Data/snATAC+RNAseq/2021_MIT_D001_snATAC-seq_PFC-MCX_ALS-FTD-PN/191112_191114_ALS_FTLD_PN_MCX_snATAC/cellranger_out/191114_PN_322_snATAC-D1/singlecell.csv",
  header = TRUE,
  row.names = 1
)

PN323_metadata <- read.csv(
  file = "/data/l3_varc/Data/snATAC+RNAseq/2021_MIT_D001_snATAC-seq_PFC-MCX_ALS-FTD-PN/191112_191114_ALS_FTLD_PN_MCX_snATAC/cellranger_out/191114_PN_323_snATAC-D2/singlecell.csv",
  header = TRUE,
  row.names = 1
)

PN324_metadata <- read.csv(
  file = "/data/l3_varc/Data/snATAC+RNAseq/2021_MIT_D001_snATAC-seq_PFC-MCX_ALS-FTD-PN/191112_191114_ALS_FTLD_PN_MCX_snATAC/cellranger_out/191114_PN_324_snATAC-D3/singlecell.csv",
  header = TRUE,
  row.names = 1
)

PN325_metadata <- read.csv(
  file = "/data/l3_varc/Data/snATAC+RNAseq/2021_MIT_D001_snATAC-seq_PFC-MCX_ALS-FTD-PN/191112_191114_ALS_FTLD_PN_MCX_snATAC/cellranger_out/191114_PN_325_snATAC-D4/singlecell.csv",
  header = TRUE,
  row.names = 1
)

PN328_metadata <- read.csv(
  file = "/data/l3_varc/Data/snATAC+RNAseq/2021_MIT_D001_snATAC-seq_PFC-MCX_ALS-FTD-PN/191112_191114_ALS_FTLD_PN_MCX_snATAC/cellranger_out/191114_PN_328_snATAC-D5/singlecell.csv",
  header = TRUE,
  row.names = 1
)


ALS110_metadata <- read.csv("/data/l3_varc/Data/snATAC+RNAseq/2021_MIT_D001_snATAC-seq_PFC-MCX_ALS-FTD-PN/191112_191114_ALS_FTLD_PN_MCX_snATAC/cellranger_out/191112_ALS_110_snATAC-H7/singlecell.csv", header = TRUE, row.names = 1)
ALS111_metadata <- read.csv("/data/l3_varc/Data/snATAC+RNAseq/2021_MIT_D001_snATAC-seq_PFC-MCX_ALS-FTD-PN/191112_191114_ALS_FTLD_PN_MCX_snATAC/cellranger_out/191112_ALS_111_snATAC-H8/singlecell.csv", header = TRUE, row.names = 1)
ALS112_metadata <- read.csv("/data/l3_varc/Data/snATAC+RNAseq/2021_MIT_D001_snATAC-seq_PFC-MCX_ALS-FTD-PN/191112_191114_ALS_FTLD_PN_MCX_snATAC/cellranger_out/191112_ALS_112_snATAC-H9/singlecell.csv", header = TRUE, row.names = 1)
ALS114_metadata <- read.csv("/data/l3_varc/Data/snATAC+RNAseq/2021_MIT_D001_snATAC-seq_PFC-MCX_ALS-FTD-PN/191112_191114_ALS_FTLD_PN_MCX_snATAC/cellranger_out/191112_ALS_114_snATAC-H10/singlecell.csv", header = TRUE, row.names = 1)
ALS115_metadata <- read.csv("/data/l3_varc/Data/snATAC+RNAseq/2021_MIT_D001_snATAC-seq_PFC-MCX_ALS-FTD-PN/191112_191114_ALS_FTLD_PN_MCX_snATAC/cellranger_out/191112_ALS_115_snATAC-H11/singlecell.csv", header = TRUE, row.names = 1)
ALS116_metadata <- read.csv("/data/l3_varc/Data/snATAC+RNAseq/2021_MIT_D001_snATAC-seq_PFC-MCX_ALS-FTD-PN/191112_191114_ALS_FTLD_PN_MCX_snATAC/cellranger_out/191112_ALS_116_snATAC-H12/singlecell.csv", header = TRUE, row.names = 1)
ALS117_metadata <- read.csv("/data/l3_varc/Data/snATAC+RNAseq/2021_MIT_D001_snATAC-seq_PFC-MCX_ALS-FTD-PN/191112_191114_ALS_FTLD_PN_MCX_snATAC/cellranger_out/191112_ALS_117_snATAC-A1/singlecell.csv", header = TRUE, row.names = 1)
ALS118_metadata <- read.csv("/data/l3_varc/Data/snATAC+RNAseq/2021_MIT_D001_snATAC-seq_PFC-MCX_ALS-FTD-PN/191112_191114_ALS_FTLD_PN_MCX_snATAC/cellranger_out/191112_ALS_118_snATAC-A2/singlecell.csv", header = TRUE, row.names = 1)
ALS120_metadata <- read.csv("/data/l3_varc/Data/snATAC+RNAseq/2021_MIT_D001_snATAC-seq_PFC-MCX_ALS-FTD-PN/191112_191114_ALS_FTLD_PN_MCX_snATAC/cellranger_out/191112_ALS_120_snATAC-A3/singlecell.csv", header = TRUE, row.names = 1)
ALS122_metadata <- read.csv("/data/l3_varc/Data/snATAC+RNAseq/2021_MIT_D001_snATAC-seq_PFC-MCX_ALS-FTD-PN/191112_191114_ALS_FTLD_PN_MCX_snATAC/cellranger_out/191112_ALS_122_snATAC-A4/singlecell.csv", header = TRUE, row.names = 1)
ALS124_metadata <- read.csv("/data/l3_varc/Data/snATAC+RNAseq/2021_MIT_D001_snATAC-seq_PFC-MCX_ALS-FTD-PN/191112_191114_ALS_FTLD_PN_MCX_snATAC/cellranger_out/191112_ALS_124_snATAC-A5/singlecell.csv", header = TRUE, row.names = 1)
ALS126_metadata <- read.csv("/data/l3_varc/Data/snATAC+RNAseq/2021_MIT_D001_snATAC-seq_PFC-MCX_ALS-FTD-PN/191112_191114_ALS_FTLD_PN_MCX_snATAC/cellranger_out/191112_ALS_126_snATAC-A6/singlecell.csv", header = TRUE, row.names = 1)
ALS127_metadata <- read.csv("/data/l3_varc/Data/snATAC+RNAseq/2021_MIT_D001_snATAC-seq_PFC-MCX_ALS-FTD-PN/191112_191114_ALS_FTLD_PN_MCX_snATAC/cellranger_out/191112_ALS_127_snATAC-A7/singlecell.csv", header = TRUE, row.names = 1)
ALS128_metadata <- read.csv("/data/l3_varc/Data/snATAC+RNAseq/2021_MIT_D001_snATAC-seq_PFC-MCX_ALS-FTD-PN/191112_191114_ALS_FTLD_PN_MCX_snATAC/cellranger_out/191112_ALS_128_snATAC-A8/singlecell.csv", header = TRUE, row.names = 1)
ALS129_metadata <- read.csv("/data/l3_varc/Data/snATAC+RNAseq/2021_MIT_D001_snATAC-seq_PFC-MCX_ALS-FTD-PN/191112_191114_ALS_FTLD_PN_MCX_snATAC/cellranger_out/191112_ALS_129_snATAC-A9/singlecell.csv", header = TRUE, row.names = 1)
ALS130_metadata <- read.csv("/data/l3_varc/Data/snATAC+RNAseq/2021_MIT_D001_snATAC-seq_PFC-MCX_ALS-FTD-PN/191112_191114_ALS_FTLD_PN_MCX_snATAC/cellranger_out/191112_ALS_130_snATAC-A10/singlecell.csv", header = TRUE, row.names = 1)



library(Signac)
library(Seurat)


PN301_chrom_assay <- CreateChromatinAssay(
  counts = PN301_counts,
  sep = c(":", "-"),
  fragments = "/data/l3_varc/Data/snATAC+RNAseq/2021_MIT_D001_snATAC-seq_PFC-MCX_ALS-FTD-PN/191112_191114_ALS_FTLD_PN_MCX_snATAC/cellranger_out/191114_PN_301_snATAC-C3/fragments.tsv.gz",
  min.cells = 10,
  min.features = 200
)

PN303_chrom_assay <- CreateChromatinAssay(
  counts = PN303_counts,
  sep = c(":", "-"),
  fragments = "/data/l3_varc/Data/snATAC+RNAseq/2021_MIT_D001_snATAC-seq_PFC-MCX_ALS-FTD-PN/191112_191114_ALS_FTLD_PN_MCX_snATAC/cellranger_out/191114_PN_303_snATAC-C4/fragments.tsv.gz",
  min.cells = 10,
  min.features = 200
)

PN304_chrom_assay <- CreateChromatinAssay(
  counts = PN304_counts,
  sep = c(":", "-"),
  fragments = "/data/l3_varc/Data/snATAC+RNAseq/2021_MIT_D001_snATAC-seq_PFC-MCX_ALS-FTD-PN/191112_191114_ALS_FTLD_PN_MCX_snATAC/cellranger_out/191114_PN_304_snATAC-C5/fragments.tsv.gz",
  min.cells = 10,
  min.features = 200
)

PN306_chrom_assay <- CreateChromatinAssay(
  counts = PN306_counts,
  sep = c(":", "-"),
  fragments = "/data/l3_varc/Data/snATAC+RNAseq/2021_MIT_D001_snATAC-seq_PFC-MCX_ALS-FTD-PN/191112_191114_ALS_FTLD_PN_MCX_snATAC/cellranger_out/191114_PN_306_snATAC-C6/fragments.tsv.gz",
  min.cells = 10,
  min.features = 200
)

PN307_chrom_assay <- CreateChromatinAssay(
  counts = PN307_counts,
  sep = c(":", "-"),
  fragments = "/data/l3_varc/Data/snATAC+RNAseq/2021_MIT_D001_snATAC-seq_PFC-MCX_ALS-FTD-PN/191112_191114_ALS_FTLD_PN_MCX_snATAC/cellranger_out/191114_PN_307_snATAC-C7/fragments.tsv.gz",
  min.cells = 10,
  min.features = 200
)

PN308_chrom_assay <- CreateChromatinAssay(
  counts = PN308_counts,
  sep = c(":", "-"),
  fragments = "/data/l3_varc/Data/snATAC+RNAseq/2021_MIT_D001_snATAC-seq_PFC-MCX_ALS-FTD-PN/191112_191114_ALS_FTLD_PN_MCX_snATAC/cellranger_out/191114_PN_308_snATAC-C8/fragments.tsv.gz",
  min.cells = 10,
  min.features = 200
)

PN309_chrom_assay <- CreateChromatinAssay(
  counts = PN309_counts,
  sep = c(":", "-"),
  fragments = "/data/l3_varc/Data/snATAC+RNAseq/2021_MIT_D001_snATAC-seq_PFC-MCX_ALS-FTD-PN/191112_191114_ALS_FTLD_PN_MCX_snATAC/cellranger_out/191114_PN_309_snATAC-C9/fragments.tsv.gz",
  min.cells = 10,
  min.features = 200
)

PN317_chrom_assay <- CreateChromatinAssay(
  counts = PN317_counts,
  sep = c(":", "-"),
  fragments = "/data/l3_varc/Data/snATAC+RNAseq/2021_MIT_D001_snATAC-seq_PFC-MCX_ALS-FTD-PN/191112_191114_ALS_FTLD_PN_MCX_snATAC/cellranger_out/191114_PN_317_snATAC-C10/fragments.tsv.gz",
  min.cells = 10,
  min.features = 200
)

PN318_chrom_assay <- CreateChromatinAssay(
  counts = PN318_counts,
  sep = c(":", "-"),
  fragments = "/data/l3_varc/Data/snATAC+RNAseq/2021_MIT_D001_snATAC-seq_PFC-MCX_ALS-FTD-PN/191112_191114_ALS_FTLD_PN_MCX_snATAC/cellranger_out/191114_PN_318_snATAC-C11/fragments.tsv.gz",
  min.cells = 10,
  min.features = 200
)

PN319_chrom_assay <- CreateChromatinAssay(
  counts = PN319_counts,
  sep = c(":", "-"),
  fragments = "/data/l3_varc/Data/snATAC+RNAseq/2021_MIT_D001_snATAC-seq_PFC-MCX_ALS-FTD-PN/191112_191114_ALS_FTLD_PN_MCX_snATAC/cellranger_out/191114_PN_319_snATAC-C12/fragments.tsv.gz",
  min.cells = 10,
  min.features = 200
)

PN322_chrom_assay <- CreateChromatinAssay(
  counts = PN322_counts,
  sep = c(":", "-"),
  fragments = "/data/l3_varc/Data/snATAC+RNAseq/2021_MIT_D001_snATAC-seq_PFC-MCX_ALS-FTD-PN/191112_191114_ALS_FTLD_PN_MCX_snATAC/cellranger_out/191114_PN_322_snATAC-D1/fragments.tsv.gz",
  min.cells = 10,
  min.features = 200
)

PN323_chrom_assay <- CreateChromatinAssay(
  counts = PN323_counts,
  sep = c(":", "-"),
  fragments = "/data/l3_varc/Data/snATAC+RNAseq/2021_MIT_D001_snATAC-seq_PFC-MCX_ALS-FTD-PN/191112_191114_ALS_FTLD_PN_MCX_snATAC/cellranger_out/191114_PN_323_snATAC-D2/fragments.tsv.gz",
  min.cells = 10,
  min.features = 200
)

PN324_chrom_assay <- CreateChromatinAssay(
  counts = PN324_counts,
  sep = c(":", "-"),
  fragments = "/data/l3_varc/Data/snATAC+RNAseq/2021_MIT_D001_snATAC-seq_PFC-MCX_ALS-FTD-PN/191112_191114_ALS_FTLD_PN_MCX_snATAC/cellranger_out/191114_PN_324_snATAC-D3/fragments.tsv.gz",
  min.cells = 10,
  min.features = 200
)

PN325_chrom_assay <- CreateChromatinAssay(
  counts = PN325_counts,
  sep = c(":", "-"),
  fragments = "/data/l3_varc/Data/snATAC+RNAseq/2021_MIT_D001_snATAC-seq_PFC-MCX_ALS-FTD-PN/191112_191114_ALS_FTLD_PN_MCX_snATAC/cellranger_out/191114_PN_325_snATAC-D4/fragments.tsv.gz",
  min.cells = 10,
  min.features = 200
)

PN328_chrom_assay <- CreateChromatinAssay(
  counts = PN328_counts,
  sep = c(":", "-"),
  fragments = "/data/l3_varc/Data/snATAC+RNAseq/2021_MIT_D001_snATAC-seq_PFC-MCX_ALS-FTD-PN/191112_191114_ALS_FTLD_PN_MCX_snATAC/cellranger_out/191114_PN_328_snATAC-D5/fragments.tsv.gz",
  min.cells = 10,
  min.features = 200
)



ALS110_chrom_assay <- CreateChromatinAssay(
  counts = ALS110_counts,
  sep = c(":", "-"),
  fragments = "/data/l3_varc/Data/snATAC+RNAseq/2021_MIT_D001_snATAC-seq_PFC-MCX_ALS-FTD-PN/191112_191114_ALS_FTLD_PN_MCX_snATAC/cellranger_out/191112_ALS_110_snATAC-H7/fragments.tsv.gz",
  min.cells = 10,
  min.features = 200
)

ALS111_chrom_assay <- CreateChromatinAssay(
  counts = ALS111_counts,
  sep = c(":", "-"),
  fragments = "/data/l3_varc/Data/snATAC+RNAseq/2021_MIT_D001_snATAC-seq_PFC-MCX_ALS-FTD-PN/191112_191114_ALS_FTLD_PN_MCX_snATAC/cellranger_out/191112_ALS_111_snATAC-H8/fragments.tsv.gz",
  min.cells = 10,
  min.features = 200
)

ALS112_chrom_assay <- CreateChromatinAssay(
  counts = ALS112_counts,
  sep = c(":", "-"),
  fragments = "/data/l3_varc/Data/snATAC+RNAseq/2021_MIT_D001_snATAC-seq_PFC-MCX_ALS-FTD-PN/191112_191114_ALS_FTLD_PN_MCX_snATAC/cellranger_out/191112_ALS_112_snATAC-H9/fragments.tsv.gz",
  min.cells = 10,
  min.features = 200
)

ALS114_chrom_assay <- CreateChromatinAssay(
  counts = ALS114_counts,
  sep = c(":", "-"),
  fragments = "/data/l3_varc/Data/snATAC+RNAseq/2021_MIT_D001_snATAC-seq_PFC-MCX_ALS-FTD-PN/191112_191114_ALS_FTLD_PN_MCX_snATAC/cellranger_out/191112_ALS_114_snATAC-H10/fragments.tsv.gz",
  min.cells = 10,
  min.features = 200
)

ALS115_chrom_assay <- CreateChromatinAssay(
  counts = ALS115_counts,
  sep = c(":", "-"),
  fragments = "/data/l3_varc/Data/snATAC+RNAseq/2021_MIT_D001_snATAC-seq_PFC-MCX_ALS-FTD-PN/191112_191114_ALS_FTLD_PN_MCX_snATAC/cellranger_out/191112_ALS_115_snATAC-H11/fragments.tsv.gz",
  min.cells = 10,
  min.features = 200
)

ALS116_chrom_assay <- CreateChromatinAssay(
  counts = ALS116_counts,
  sep = c(":", "-"),
  fragments = "/data/l3_varc/Data/snATAC+RNAseq/2021_MIT_D001_snATAC-seq_PFC-MCX_ALS-FTD-PN/191112_191114_ALS_FTLD_PN_MCX_snATAC/cellranger_out/191112_ALS_116_snATAC-H12/fragments.tsv.gz",
  min.cells = 10,
  min.features = 200
)

ALS117_chrom_assay <- CreateChromatinAssay(
  counts = ALS117_counts,
  sep = c(":", "-"),
  fragments = "/data/l3_varc/Data/snATAC+RNAseq/2021_MIT_D001_snATAC-seq_PFC-MCX_ALS-FTD-PN/191112_191114_ALS_FTLD_PN_MCX_snATAC/cellranger_out/191112_ALS_117_snATAC-A1/fragments.tsv.gz",
  min.cells = 10,
  min.features = 200
)

ALS118_chrom_assay <- CreateChromatinAssay(
  counts = ALS118_counts,
  sep = c(":", "-"),
  fragments = "/data/l3_varc/Data/snATAC+RNAseq/2021_MIT_D001_snATAC-seq_PFC-MCX_ALS-FTD-PN/191112_191114_ALS_FTLD_PN_MCX_snATAC/cellranger_out/191112_ALS_118_snATAC-A2/fragments.tsv.gz",
  min.cells = 10,
  min.features = 200
)

ALS120_chrom_assay <- CreateChromatinAssay(
  counts = ALS120_counts,
  sep = c(":", "-"),
  fragments = "/data/l3_varc/Data/snATAC+RNAseq/2021_MIT_D001_snATAC-seq_PFC-MCX_ALS-FTD-PN/191112_191114_ALS_FTLD_PN_MCX_snATAC/cellranger_out/191112_ALS_120_snATAC-A3/fragments.tsv.gz",
  min.cells = 10,
  min.features = 200
)

ALS122_chrom_assay <- CreateChromatinAssay(
  counts = ALS122_counts,
  sep = c(":", "-"),
  fragments = "/data/l3_varc/Data/snATAC+RNAseq/2021_MIT_D001_snATAC-seq_PFC-MCX_ALS-FTD-PN/191112_191114_ALS_FTLD_PN_MCX_snATAC/cellranger_out/191112_ALS_122_snATAC-A4/fragments.tsv.gz",
  min.cells = 10,
  min.features = 200
)

ALS124_chrom_assay <- CreateChromatinAssay(
  counts = ALS124_counts,
  sep = c(":", "-"),
  fragments = "/data/l3_varc/Data/snATAC+RNAseq/2021_MIT_D001_snATAC-seq_PFC-MCX_ALS-FTD-PN/191112_191114_ALS_FTLD_PN_MCX_snATAC/cellranger_out/191112_ALS_124_snATAC-A5/fragments.tsv.gz",
  min.cells = 10,
  min.features = 200
)

ALS126_chrom_assay <- CreateChromatinAssay(
  counts = ALS126_counts,
  sep = c(":", "-"),
  fragments = "/data/l3_varc/Data/snATAC+RNAseq/2021_MIT_D001_snATAC-seq_PFC-MCX_ALS-FTD-PN/191112_191114_ALS_FTLD_PN_MCX_snATAC/cellranger_out/191112_ALS_126_snATAC-A6/fragments.tsv.gz",
  min.cells = 10,
  min.features = 200
)

ALS127_chrom_assay <- CreateChromatinAssay(
  counts = ALS127_counts,
  sep = c(":", "-"),
  fragments = "/data/l3_varc/Data/snATAC+RNAseq/2021_MIT_D001_snATAC-seq_PFC-MCX_ALS-FTD-PN/191112_191114_ALS_FTLD_PN_MCX_snATAC/cellranger_out/191112_ALS_127_snATAC-A7/fragments.tsv.gz",
  min.cells = 10,
  min.features = 200
)

ALS128_chrom_assay <- CreateChromatinAssay(
  counts = ALS128_counts,
  sep = c(":", "-"),
  fragments = "/data/l3_varc/Data/snATAC+RNAseq/2021_MIT_D001_snATAC-seq_PFC-MCX_ALS-FTD-PN/191112_191114_ALS_FTLD_PN_MCX_snATAC/cellranger_out/191112_ALS_128_snATAC-A8/fragments.tsv.gz",
  min.cells = 10,
  min.features = 200
)

ALS129_chrom_assay <- CreateChromatinAssay(
  counts = ALS129_counts,
  sep = c(":", "-"),
  fragments = "/data/l3_varc/Data/snATAC+RNAseq/2021_MIT_D001_snATAC-seq_PFC-MCX_ALS-FTD-PN/191112_191114_ALS_FTLD_PN_MCX_snATAC/cellranger_out/191112_ALS_129_snATAC-A9/fragments.tsv.gz",
  min.cells = 10,
  min.features = 200
)

ALS130_chrom_assay <- CreateChromatinAssay(
  counts = ALS130_counts,
  sep = c(":", "-"),
  fragments = "/data/l3_varc/Data/snATAC+RNAseq/2021_MIT_D001_snATAC-seq_PFC-MCX_ALS-FTD-PN/191112_191114_ALS_FTLD_PN_MCX_snATAC/cellranger_out/191112_ALS_130_snATAC-A10/fragments.tsv.gz",
  min.cells = 10,
  min.features = 200
)




PN301 <- CreateSeuratObject(
  counts = PN301_chrom_assay,
  assay = "peaks",
  meta.data = PN301_metadata
)

PN303 <- CreateSeuratObject(
  counts = PN303_chrom_assay,
  assay = "peaks",
  meta.data = PN303_metadata
)

PN304 <- CreateSeuratObject(
  counts = PN304_chrom_assay,
  assay = "peaks",
  meta.data = PN304_metadata
)

PN306 <- CreateSeuratObject(
  counts = PN306_chrom_assay,
  assay = "peaks",
  meta.data = PN306_metadata
)

PN307 <- CreateSeuratObject(
  counts = PN307_chrom_assay,
  assay = "peaks",
  meta.data = PN307_metadata
)

PN308 <- CreateSeuratObject(
  counts = PN308_chrom_assay,
  assay = "peaks",
  meta.data = PN308_metadata
)

PN309 <- CreateSeuratObject(
  counts = PN309_chrom_assay,
  assay = "peaks",
  meta.data = PN309_metadata
)

PN317 <- CreateSeuratObject(
  counts = PN317_chrom_assay,
  assay = "peaks",
  meta.data = PN317_metadata
)

PN318 <- CreateSeuratObject(
  counts = PN318_chrom_assay,
  assay = "peaks",
  meta.data = PN318_metadata
)

PN319 <- CreateSeuratObject(
  counts = PN319_chrom_assay,
  assay = "peaks",
  meta.data = PN319_metadata
)

PN322 <- CreateSeuratObject(
  counts = PN322_chrom_assay,
  assay = "peaks",
  meta.data = PN322_metadata
)

PN323 <- CreateSeuratObject(
  counts = PN323_chrom_assay,
  assay = "peaks",
  meta.data = PN323_metadata
)

PN324 <- CreateSeuratObject(
  counts = PN324_chrom_assay,
  assay = "peaks",
  meta.data = PN324_metadata
)

PN325 <- CreateSeuratObject(
  counts = PN325_chrom_assay,
  assay = "peaks",
  meta.data = PN325_metadata
)

PN328 <- CreateSeuratObject(
  counts = PN328_chrom_assay,
  assay = "peaks",
  meta.data = PN328_metadata
)



ALS110 <- CreateSeuratObject(counts = ALS110_chrom_assay, assay = "peaks", meta.data = ALS110_metadata)
ALS111 <- CreateSeuratObject(counts = ALS111_chrom_assay, assay = "peaks", meta.data = ALS111_metadata)
ALS112 <- CreateSeuratObject(counts = ALS112_chrom_assay, assay = "peaks", meta.data = ALS112_metadata)
ALS114 <- CreateSeuratObject(counts = ALS114_chrom_assay, assay = "peaks", meta.data = ALS114_metadata)
ALS115 <- CreateSeuratObject(counts = ALS115_chrom_assay, assay = "peaks", meta.data = ALS115_metadata)
ALS116 <- CreateSeuratObject(counts = ALS116_chrom_assay, assay = "peaks", meta.data = ALS116_metadata)
ALS117 <- CreateSeuratObject(counts = ALS117_chrom_assay, assay = "peaks", meta.data = ALS117_metadata)
ALS118 <- CreateSeuratObject(counts = ALS118_chrom_assay, assay = "peaks", meta.data = ALS118_metadata)
ALS120 <- CreateSeuratObject(counts = ALS120_chrom_assay, assay = "peaks", meta.data = ALS120_metadata)
ALS122 <- CreateSeuratObject(counts = ALS122_chrom_assay, assay = "peaks", meta.data = ALS122_metadata)
ALS124 <- CreateSeuratObject(counts = ALS124_chrom_assay, assay = "peaks", meta.data = ALS124_metadata)
ALS126 <- CreateSeuratObject(counts = ALS126_chrom_assay, assay = "peaks", meta.data = ALS126_metadata)
ALS127 <- CreateSeuratObject(counts = ALS127_chrom_assay, assay = "peaks", meta.data = ALS127_metadata)
ALS128 <- CreateSeuratObject(counts = ALS128_chrom_assay, assay = "peaks", meta.data = ALS128_metadata)
ALS129 <- CreateSeuratObject(counts = ALS129_chrom_assay, assay = "peaks", meta.data = ALS129_metadata)
ALS130 <- CreateSeuratObject(counts = ALS130_chrom_assay, assay = "peaks", meta.data = ALS130_metadata)




# Check the Seurat objects
PN301
PN303
PN304
PN306
PN307
PN308
PN309
PN317
PN318
PN319
PN322
PN323
PN324
PN325
PN328


ALS110
ALS111
ALS112
ALS114
ALS115
ALS116
ALS117
ALS118
ALS120
ALS122
ALS124
ALS126
ALS127
ALS128
ALS129
ALS130

# Cleaning up memory

# counts
rm(PN301_counts)
rm(PN303_counts)
rm(PN304_counts)
rm(PN306_counts)
rm(PN307_counts)
rm(PN308_counts)
rm(PN309_counts)
rm(PN317_counts)
rm(PN318_counts)
rm(PN319_counts)
rm(PN322_counts)
rm(PN323_counts)
rm(PN324_counts)
rm(PN325_counts)
rm(PN328_counts)

# metadata
rm(PN301_metadata)
rm(PN303_metadata)
rm(PN304_metadata)
rm(PN306_metadata)
rm(PN307_metadata)
rm(PN308_metadata)
rm(PN309_metadata)
rm(PN317_metadata)
rm(PN318_metadata)
rm(PN319_metadata)
rm(PN322_metadata)
rm(PN323_metadata)
rm(PN324_metadata)
rm(PN325_metadata)
rm(PN328_metadata)

# chromatin assays
rm(PN301_chrom_assay)
rm(PN303_chrom_assay)
rm(PN304_chrom_assay)
rm(PN306_chrom_assay)
rm(PN307_chrom_assay)
rm(PN308_chrom_assay)
rm(PN309_chrom_assay)
rm(PN317_chrom_assay)
rm(PN318_chrom_assay)
rm(PN319_chrom_assay)
rm(PN322_chrom_assay)
rm(PN323_chrom_assay)
rm(PN324_chrom_assay)
rm(PN325_chrom_assay)
rm(PN328_chrom_assay)

# counts
rm(ALS110_counts)
rm(ALS111_counts)
rm(ALS112_counts)
rm(ALS114_counts)
rm(ALS115_counts)
rm(ALS116_counts)
rm(ALS117_counts)
rm(ALS118_counts)
rm(ALS120_counts)
rm(ALS122_counts)
rm(ALS124_counts)
rm(ALS126_counts)
rm(ALS127_counts)
rm(ALS128_counts)
rm(ALS129_counts)
rm(ALS130_counts)

# metadata
rm(ALS110_metadata)
rm(ALS111_metadata)
rm(ALS112_metadata)
rm(ALS114_metadata)
rm(ALS115_metadata)
rm(ALS116_metadata)
rm(ALS117_metadata)
rm(ALS118_metadata)
rm(ALS120_metadata)
rm(ALS122_metadata)
rm(ALS124_metadata)
rm(ALS126_metadata)
rm(ALS127_metadata)
rm(ALS128_metadata)
rm(ALS129_metadata)
rm(ALS130_metadata)

# chromatin assays
rm(ALS110_chrom_assay)
rm(ALS111_chrom_assay)
rm(ALS112_chrom_assay)
rm(ALS114_chrom_assay)
rm(ALS115_chrom_assay)
rm(ALS116_chrom_assay)
rm(ALS117_chrom_assay)
rm(ALS118_chrom_assay)
rm(ALS120_chrom_assay)
rm(ALS122_chrom_assay)
rm(ALS124_chrom_assay)
rm(ALS126_chrom_assay)
rm(ALS127_chrom_assay)
rm(ALS128_chrom_assay)
rm(ALS129_chrom_assay)
rm(ALS130_chrom_assay)



# Checking specifically the ATAC-seq data in the Seurat obj, stored using the custom assay "ChromatinAssay"

PN301[['peaks']]
PN303[['peaks']]
PN304[['peaks']]
PN306[['peaks']]
PN307[['peaks']]
PN308[['peaks']]
PN309[['peaks']]
PN317[['peaks']]
PN318[['peaks']]
PN319[['peaks']]
PN322[['peaks']]
PN323[['peaks']]
PN324[['peaks']]
PN325[['peaks']]
PN328[['peaks']]

ALS110[['peaks']]
ALS111[['peaks']]
ALS112[['peaks']]
ALS114[['peaks']]
ALS115[['peaks']]
ALS116[['peaks']]
ALS117[['peaks']]
ALS118[['peaks']]
ALS120[['peaks']]
ALS122[['peaks']]
ALS124[['peaks']]
ALS126[['peaks']]
ALS127[['peaks']]
ALS128[['peaks']]
ALS129[['peaks']]
ALS130[['peaks']]


# One example of what can be done with ATAC-seq data stored as a ChromatinAssay:
# `granges` to see the genomic ranges associated with each feature in the object

granges(PN301)
granges(PN303)
granges(PN304)
granges(PN306)
granges(PN307)
granges(PN308)
granges(PN309)
granges(PN317)
granges(PN318)
granges(PN319)
granges(PN322)
granges(PN323)
granges(PN324)
granges(PN325)
granges(PN328)


granges(ALS110)
granges(ALS111)
granges(ALS112)
granges(ALS114)
granges(ALS115)
granges(ALS116)
granges(ALS117)
granges(ALS118)
granges(ALS120)
granges(ALS122)
granges(ALS124)
granges(ALS126)
granges(ALS127)
granges(ALS128)
granges(ALS129)
granges(ALS130)



# With the results from granges, we can now remove features that correspond 
# to chromosome scaffolds (e.g. KI270713.1) or other sequences instead of 
# the (22+2) standard chromosomes.


peaks.keep <- seqnames(granges(PN301)) %in% standardChromosomes(granges(PN301))
PN301 <- PN301[as.vector(peaks.keep), ]
PN301

peaks.keep <- seqnames(granges(PN303)) %in% standardChromosomes(granges(PN303))
PN303 <- PN303[as.vector(peaks.keep), ]
PN303

peaks.keep <- seqnames(granges(PN304)) %in% standardChromosomes(granges(PN304))
PN304 <- PN304[as.vector(peaks.keep), ]
PN304

peaks.keep <- seqnames(granges(PN306)) %in% standardChromosomes(granges(PN306))
PN306 <- PN306[as.vector(peaks.keep), ]
PN306

peaks.keep <- seqnames(granges(PN307)) %in% standardChromosomes(granges(PN307))
PN307 <- PN307[as.vector(peaks.keep), ]
PN307

peaks.keep <- seqnames(granges(PN308)) %in% standardChromosomes(granges(PN308))
PN308 <- PN308[as.vector(peaks.keep), ]
PN308

peaks.keep <- seqnames(granges(PN309)) %in% standardChromosomes(granges(PN309))
PN309 <- PN309[as.vector(peaks.keep), ]
PN309

peaks.keep <- seqnames(granges(PN317)) %in% standardChromosomes(granges(PN317))
PN317 <- PN317[as.vector(peaks.keep), ]
PN317

peaks.keep <- seqnames(granges(PN318)) %in% standardChromosomes(granges(PN318))
PN318 <- PN318[as.vector(peaks.keep), ]
PN318

peaks.keep <- seqnames(granges(PN319)) %in% standardChromosomes(granges(PN319))
PN319 <- PN319[as.vector(peaks.keep), ]
PN319

peaks.keep <- seqnames(granges(PN322)) %in% standardChromosomes(granges(PN322))
PN322 <- PN322[as.vector(peaks.keep), ]
PN322

peaks.keep <- seqnames(granges(PN323)) %in% standardChromosomes(granges(PN323))
PN323 <- PN323[as.vector(peaks.keep), ]
PN323

peaks.keep <- seqnames(granges(PN324)) %in% standardChromosomes(granges(PN324))
PN324 <- PN324[as.vector(peaks.keep), ]
PN324

peaks.keep <- seqnames(granges(PN325)) %in% standardChromosomes(granges(PN325))
PN325 <- PN325[as.vector(peaks.keep), ]
PN325

peaks.keep <- seqnames(granges(PN328)) %in% standardChromosomes(granges(PN328))
PN328 <- PN328[as.vector(peaks.keep), ]
PN328




peaks.keep <- seqnames(granges(ALS110)) %in% standardChromosomes(granges(ALS110))
ALS110 <- ALS110[as.vector(peaks.keep), ]
ALS110

peaks.keep <- seqnames(granges(ALS111)) %in% standardChromosomes(granges(ALS111))
ALS111 <- ALS111[as.vector(peaks.keep), ]
ALS111

peaks.keep <- seqnames(granges(ALS112)) %in% standardChromosomes(granges(ALS112))
ALS112 <- ALS112[as.vector(peaks.keep), ]
ALS112

peaks.keep <- seqnames(granges(ALS114)) %in% standardChromosomes(granges(ALS114))
ALS114 <- ALS114[as.vector(peaks.keep), ]
ALS114

peaks.keep <- seqnames(granges(ALS115)) %in% standardChromosomes(granges(ALS115))
ALS115 <- ALS115[as.vector(peaks.keep), ]
ALS115

peaks.keep <- seqnames(granges(ALS116)) %in% standardChromosomes(granges(ALS116))
ALS116 <- ALS116[as.vector(peaks.keep), ]
ALS116

peaks.keep <- seqnames(granges(ALS117)) %in% standardChromosomes(granges(ALS117))
ALS117 <- ALS117[as.vector(peaks.keep), ]
ALS117

peaks.keep <- seqnames(granges(ALS118)) %in% standardChromosomes(granges(ALS118))
ALS118 <- ALS118[as.vector(peaks.keep), ]
ALS118

peaks.keep <- seqnames(granges(ALS120)) %in% standardChromosomes(granges(ALS120))
ALS120 <- ALS120[as.vector(peaks.keep), ]
ALS120

peaks.keep <- seqnames(granges(ALS122)) %in% standardChromosomes(granges(ALS122))
ALS122 <- ALS122[as.vector(peaks.keep), ]
ALS122

peaks.keep <- seqnames(granges(ALS124)) %in% standardChromosomes(granges(ALS124))
ALS124 <- ALS124[as.vector(peaks.keep), ]
ALS124

peaks.keep <- seqnames(granges(ALS126)) %in% standardChromosomes(granges(ALS126))
ALS126 <- ALS126[as.vector(peaks.keep), ]
ALS126

peaks.keep <- seqnames(granges(ALS127)) %in% standardChromosomes(granges(ALS127))
ALS127 <- ALS127[as.vector(peaks.keep), ]
ALS127

peaks.keep <- seqnames(granges(ALS128)) %in% standardChromosomes(granges(ALS128))
ALS128 <- ALS128[as.vector(peaks.keep), ]
ALS128

peaks.keep <- seqnames(granges(ALS129)) %in% standardChromosomes(granges(ALS129))
ALS129 <- ALS129[as.vector(peaks.keep), ]
ALS129

peaks.keep <- seqnames(granges(ALS130)) %in% standardChromosomes(granges(ALS130))
ALS130 <- ALS130[as.vector(peaks.keep), ]
ALS130



#You can add gene annotation data to the Seurat object, so that this information can be pulled
#directly from the object in downstream steps (as opposed to needing to pull from other sources)

#This will be done across the next few code chunks.
#This first one checks AnnotationHub for the details surrounding the specific reference accession.

#*WARNING!!!* There are lots of different patches released for each genome assembly. It is very important that
#you use/attach the same patch that was used in the initial mapping of the data when it was orignally processed.

#For example, for this tutorial code from 10x Genomics, the original data was mapped using GRCh38-2020-A (as found in the dataset summary: https://cf.10xgenomics.com/samples/cell-atac/2.1.0/10k_pbmc_ATACv2_nextgem_Chromium_Controller/10k_pbmc_ATACv2_nextgem_Chromium_Controller_web_summary.html)
#That particular reference package corresponds to the Ensembl v98 patch release, so that is what we have to add to this Seurat object.


# if (!require("BiocManager", quietly = TRUE))
#     install.packages("BiocManager")
# 
# BiocManager::install("AnnotationHub")

library(AnnotationHub)
#In the following line, you may need to give it permission to write a new file into your R directory
ah <- AnnotationHub()

# Search for the Ensembl 92 EnsDb for Homo sapiens on AnnotationHub
query(ah, "EnsDb.Hsapiens.v92")



#Save this AnnotationHub reference as an object
ensdb_v92 <- ah[["AH60977"]]



#BiocManager::install("biovizBase")

# extract gene annotations from EnsDb
annotations <- GetGRangesFromEnsDb(ensdb = ensdb_v92)

# change to UCSC style since the data was mapped to hg38
seqlevels(annotations) <- paste0('chr', seqlevels(annotations))
genome(annotations) <- "hg38"



# Add the gene information to the object


Annotation(PN301) <- annotations
head(PN301)
tail(PN301)

Annotation(PN303) <- annotations
head(PN303)
tail(PN303)

Annotation(PN304) <- annotations
head(PN304)
tail(PN304)

Annotation(PN306) <- annotations
head(PN306)
tail(PN306)

Annotation(PN307) <- annotations
head(PN307)
tail(PN307)

Annotation(PN308) <- annotations
head(PN308)
tail(PN308)

Annotation(PN309) <- annotations
head(PN309)
tail(PN309)

Annotation(PN317) <- annotations
head(PN317)
tail(PN317)

Annotation(PN318) <- annotations
head(PN318)
tail(PN318)

Annotation(PN319) <- annotations
head(PN319)
tail(PN319)

Annotation(PN322) <- annotations
head(PN322)
tail(PN322)

Annotation(PN323) <- annotations
head(PN323)
tail(PN323)

Annotation(PN324) <- annotations
head(PN324)
tail(PN324)

Annotation(PN325) <- annotations
head(PN325)
tail(PN325)

Annotation(PN328) <- annotations
head(PN328)
tail(PN328)




Annotation(ALS110) <- annotations
head(ALS110)
tail(ALS110)

Annotation(ALS111) <- annotations
head(ALS111)
tail(ALS111)

Annotation(ALS112) <- annotations
head(ALS112)
tail(ALS112)

Annotation(ALS114) <- annotations
head(ALS114)
tail(ALS114)

Annotation(ALS115) <- annotations
head(ALS115)
tail(ALS115)

Annotation(ALS116) <- annotations
head(ALS116)
tail(ALS116)

Annotation(ALS117) <- annotations
head(ALS117)
tail(ALS117)

Annotation(ALS118) <- annotations
head(ALS118)
tail(ALS118)

Annotation(ALS120) <- annotations
head(ALS120)
tail(ALS120)

Annotation(ALS122) <- annotations
head(ALS122)
tail(ALS122)

Annotation(ALS124) <- annotations
head(ALS124)
tail(ALS124)

Annotation(ALS126) <- annotations
head(ALS126)
tail(ALS126)

Annotation(ALS127) <- annotations
head(ALS127)
tail(ALS127)

Annotation(ALS128) <- annotations
head(ALS128)
tail(ALS128)

Annotation(ALS129) <- annotations
head(ALS129)
tail(ALS129)

Annotation(ALS130) <- annotations
head(ALS130)
tail(ALS130)



# Adding sample data to each Seurat object

PN301$dataset <- 'PN301'
PN303$dataset <- 'PN303'
PN304$dataset <- 'PN304'
PN306$dataset <- 'PN306'
PN307$dataset <- 'PN307'
PN308$dataset <- 'PN308'
PN309$dataset <- 'PN309'
PN317$dataset <- 'PN317'
PN318$dataset <- 'PN318'
PN319$dataset <- 'PN319'
PN322$dataset <- 'PN322'
PN323$dataset <- 'PN323'
PN324$dataset <- 'PN324'
PN325$dataset <- 'PN325'
PN328$dataset <- 'PN328'


ALS110$dataset <- 'ALS110'
ALS111$dataset <- 'ALS111'
ALS112$dataset <- 'ALS112'
ALS114$dataset <- 'ALS114'
ALS115$dataset <- 'ALS115'
ALS116$dataset <- 'ALS116'
ALS117$dataset <- 'ALS117'
ALS118$dataset <- 'ALS118'
ALS120$dataset <- 'ALS120'
ALS122$dataset <- 'ALS122'
ALS124$dataset <- 'ALS124'
ALS126$dataset <- 'ALS126'
ALS127$dataset <- 'ALS127'
ALS128$dataset <- 'ALS128'
ALS129$dataset <- 'ALS129'
ALS130$dataset <- 'ALS130'



# Merging all samples
unintegrated_ALS_PN <- merge(
  x = PN301, 
  y = c(PN303, PN304, PN306, PN307, PN308, PN309, PN317, PN318, PN319, PN322, PN323, PN324, PN325, PN328, ALS110, ALS111, ALS112, ALS114, ALS115, ALS116, ALS117, ALS118, ALS120, ALS122, ALS124, ALS126, ALS127, ALS128, ALS129, ALS130),
  add.cell.ids = c("PN301", "PN303", "PN304", "PN306", "PN307", "PN308", "PN309", "PN317", "PN318", "PN319", "PN322", "PN323", "PN324", "PN325", "PN328", "ALS110", "ALS111", "ALS112", "ALS114", "ALS115", "ALS116", "ALS117", "ALS118", "ALS120", "ALS122", "ALS124", "ALS126", "ALS127", "ALS128", "ALS129", "ALS130")
)



# Cleaning memory
rm(PN301, PN303, PN304, PN306, PN307, PN308, PN309, PN317, PN318, PN319, PN322, PN323, PN324, PN325, PN328, ALS110, ALS111, ALS112, ALS114, ALS115, ALS116, ALS117, ALS118, ALS120, ALS122, ALS124, ALS126, ALS127, ALS128, ALS129, ALS130)



#Saving merged sample
saveRDS(unintegrated_ALS_PN, file = "unintegrated_ALS_PN.rds")



#Computing QC Metrics

#Explanation for each of these steps in the vignette
#Each one of these adds their values into the Seurat object

# compute nucleosome signal score per cell
unintegrated_ALS_PN <- NucleosomeSignal(object = unintegrated_ALS_PN)

# compute TSS enrichment score per cell
unintegrated_ALS_PN <- TSSEnrichment(object = unintegrated_ALS_PN)

# add fraction of reads in peaks
unintegrated_ALS_PN$pct_reads_in_peaks <- unintegrated_ALS_PN$peak_region_fragments / unintegrated_ALS_PN$passed_filters * 100

# add blacklist ratio
unintegrated_ALS_PN$blacklist_ratio <- FractionCountsInRegion(
  object = unintegrated_ALS_PN, 
  assay = 'peaks',
  regions = blacklist_hg38_unified
)



#Saving merged sample with QC
saveRDS(unintegrated_ALS_PN, file = "unintegrated_ALS_QC_PN.rds")



#Visualizing the relationship between these different variables in the metadata

DensityScatter <- DensityScatter(unintegrated_ALS_PN, x = 'nCount_peaks', y = 'TSS.enrichment', log_x = TRUE, quantiles = TRUE)
ggsave("DensityScatter.png", plot = DensityScatter, width = 12, height = 6, dpi = 300)
DensityScatter



#Filtering out cells that are outliers for the mononucleosomal/nucleosome-free ratio

unintegrated_ALS_PN$nucleosome_group <- ifelse(unintegrated_ALS_PN$nucleosome_signal > 3, 'NS > 3', 'NS < 3')
FragmentHistogram(object = unintegrated_ALS_PN, group.by = 'nucleosome_group')

#Note the y-axes on these plots
#This plot looks very different from the tutorial



#Plotting the distribution of each QC metric separately with violin plots

unintegrated_ALS_PN_QC_ViolinPlots <- VlnPlot(
  object = unintegrated_ALS_PN,
  features = c('nCount_peaks', 'TSS.enrichment', 'blacklist_ratio', 'nucleosome_signal', 'pct_reads_in_peaks'),
  pt.size = 0.1,
  ncol = 5
)

ggsave("unintegrated_ALS_PN_QC_ViolinPlots.png", plot = unintegrated_ALS_PN_QC_ViolinPlots, width = 12, height = 6, dpi = 300)

unintegrated_ALS_PN_QC_ViolinPlots

#The resolution on these is terrible in this viewer, but if you save the image as a .png they display much better.

#This plot also looks very different from the tutorial data.



#Remove cells based on outliers for these QC metrics.
#The exact thresholds for these can and should all be tailored based on your dataset

unintegrated_ALS_PN <- subset(
  x = unintegrated_ALS_PN,
  subset = nCount_peaks > 9000 &
    nCount_peaks < 100000 &
    pct_reads_in_peaks > 40 &
    blacklist_ratio < 0.01 &
    nucleosome_signal < 3 &
    TSS.enrichment > 3
)
unintegrated_ALS_PN

#Saving filted object
saveRDS(unintegrated_ALS_PN, file = "unintegrated_ALS_PN_QC_Filtered.rds")



#Normalization and Linear Dimensional Reduction
#Details for both of these can be found in the vignette

#NOTE: One thing that's important to consider is the cutoff threshold for features.
#In this example, they are NOT excluding any features!
#However, in your own samples, you may want to set this threshold higher (ex. 'q75' to use only the top 25% of peaks), as it will run faster.

unintegrated_ALS_PN <- RunTFIDF(unintegrated_ALS_PN) #This is the normalization step: Term frequency-inverse document frequency normalization
unintegrated_ALS_PN <- FindTopFeatures(unintegrated_ALS_PN, min.cutoff = 'q0') #This is the feature selection step: min.cutoff sets the percentile threshold
unintegrated_ALS_PN <- RunSVD(unintegrated_ALS_PN) #This is the dimension reduction step: singular value decomposition

#The combination of TF-IDF and SVD is known as "latent semantic indexing" (LSI). After running these steps, you should see that appear in the Seurat obj
unintegrated_ALS_PN

#The output from all of these steps is *conceptually* similar to what you get after PCA in scRNAseq analysis

#Saving object
saveRDS(unintegrated_ALS_PN, file = "unintegrated_ALS_PN_LinDimRed.rds")



#The first LSI component frequently captures sequencing depth (technical variation) rather that biological variation.
#If this is the case, it should be removed from downstream analysis.
#We can assess the correlation between each LSI component and sequencing depth with the following function:

DepthCor <- DepthCor(unintegrated_ALS_PN)
ggsave("DepthCor.png", plot = DepthCor, width = 12, height = 6, dpi = 300)
DepthCor



#Looking at the above plot for the tutorial data, it is clear that the first LSI component is basically just capturing sequencing depth,
#thus it should not be used in downstream analyses.



#Non-linear dimension reduction and clustering

#These steps perform nonlinear and graph-based clustering on the now-relatively-low-dimensional data.
#These steps all come from the Seurat pacakge and are the same as what are used in scRNAseq analysis

unintegrated_ALS_PN <- RunUMAP(object = unintegrated_ALS_PN, reduction = 'lsi', dims = 2:30)
unintegrated_ALS_PN <- FindNeighbors(object = unintegrated_ALS_PN, reduction = 'lsi', dims = 2:30)
unintegrated_ALS_PN <- FindClusters(object = unintegrated_ALS_PN, verbose = FALSE, algorithm = 3)
unintegrated_ALS_PN_clusters <- DimPlot(object = unintegrated_ALS_PN, label = TRUE) + NoLegend()
unintegrated_ALS_PN_overlap <- DimPlot(unintegrated_ALS_PN, group.by = 'dataset', pt.size = 0.1)

ggsave("unintegrated_ALS_PN_clusters.png", plot = unintegrated_ALS_PN_clusters, width = 12, height = 6, dpi = 300)
ggsave("unintegrated_ALS_PN_overlap.png", plot = unintegrated_ALS_PN_overlap, width = 12, height = 6, dpi = 300)

unintegrated_ALS_PN_clusters
unintegrated_ALS_PN_overlap

#Saving rds
saveRDS(unintegrated_ALS_PN, file = "unintegrated_ALS_PN_NonLimDimRed.rds")


#Nelow step only required if starting from an intermediary .rds
#unintegrated_ALS_PN <- readRDS("/home/guidubjl/Signac/results_integrated_ALS_PN/unintegrated_ALS_PN_NonLimDimRed.rds")


#Batch correct using Harmony
hm_integrated_ALS_PN <- RunHarmony(object = unintegrated_ALS_PN, group.by.vars = 'dataset', reduction.use = 'lsi', assay.use = 'peaks', project.dim = FALSE)
hm_integrated_ALS_PN <- RunUMAP(hm_integrated_ALS_PN, dims = 2:30, reduction = 'harmony')
hm_integrated_ALS_PN_UMAP_Plot <- DimPlot(hm_integrated_ALS_PN, group.by = 'dataset', pt.size = 0.1)

rm(unintegrated_ALS_PN)

ggsave("hm_integrated_ALS_PN_UMAP_Plot.png", plot = hm_integrated_ALS_PN_UMAP_Plot, width = 12, height = 6, dpi = 300)

hm_integrated_ALS_PN_UMAP_Plot



sessionInfo <- sessionInfo()
saveRDS(sessionInfo, file = "sessionInfo.rds")
saveRDS(hm_integrated_ALS_PN, file = "hm_integrated_ALS_PN.rds")

#Below step only required if starting from an intermediary .rds
#hm_integrated_ALS_PN <- readRDS("/home/guidubjl/Signac/results_integrated_ALS_PN/hm_integrated_ALS_PN.rds")


gene.activities <- GeneActivity(hm_integrated_ALS_PN)


# add the gene activity matrix to the Seurat object as a new assay and normalize it
hm_integrated_ALS_PN[['RNA']] <- CreateAssayObject(counts = gene.activities)
hm_integrated_ALS_PN <- NormalizeData(
  object = hm_integrated_ALS_PN,
  assay = 'RNA',
  normalization.method = 'LogNormalize',
  scale.factor = median(hm_integrated_ALS_PN$nCount_RNA)
)


DefaultAssay(hm_integrated_ALS_PN) <- 'RNA'

glymph_genes_clusters <- FeaturePlot(
  object = hm_integrated_ALS_PN,
  features = c('GFAP', 'AQP4', 'VEGFA', 'VEGFB', 'SHTN1'),
  pt.size = 0.1,
  max.cutoff = 'q95',
  ncol = 3
)

ggsave("glymph_genes_clusters.png", plot = glymph_genes_clusters, width = 12, height = 6, dpi = 300)
glymph_genes_clusters

saveRDS(hm_integrated_ALS_PN, file = "hm_integrated_ALS_PN_GeneActivity.rds")

## Integrating with scRNA-seq data


# Load the pre-processed scRNA-seq data
library(dplyr)
library(purrr)
library(SummarizedExperiment)


ALS_rna <- readRDS("/data/l3_varc/Data/snATAC+RNAseq/2021_MIT_D002_snRNA-seq_PFC-MCX/se_agg_BA4_fil/se_agg_BA4_SALS_fil.rds")
PN_rna <- readRDS("/data/l3_varc/Data/snATAC+RNAseq/2021_MIT_D002_snRNA-seq_PFC-MCX/se_agg_BA4_fil/se_agg_BA4_PN_fil.rds")
merged_rna_data <- cbind(ALS_rna, PN_rna)

merged_rna_data_counts <- SummarizedExperiment::assay(merged_rna_data, "counts")
rownames(merged_rna_data_counts) <- rowRanges(merged_rna_data)@elementMetadata$Gene
merged_rna_data_counts <- merged_rna_data_counts[rownames(merged_rna_data) != "N/A", ]
merged_rna_data_counts <- merged_rna_data_counts[match(unique(rownames(merged_rna_data_counts)), rownames(merged_rna_data_counts)),]
merged_rna_data_metadata <- as.data.frame(SummarizedExperiment::colData(merged_rna_data))

merged_rna_data_seurat <- CreateSeuratObject(counts = merged_rna_data_counts, meta.data = merged_rna_data_metadata, assay = "RNA")

merged_rna_data_seurat_updated <- UpdateSeuratObject(merged_rna_data_seurat)

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
# packageVersion('SeuratData')
# library(SeuratData)
# #remotes::install_github("satijalab/azimuth", ref = 'master')
# packageVersion('Azimuth')
# library(Azimuth)

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

merged_rna_data_complete <- merged_rna_data_seurat_updated

saveRDS(merged_rna_data_complete, file = "merged_rna_data_complete.rds")

rm(merged_rna_data)
rm(merged_rna_data_seurat)
rm(merged_rna_data_seurat_updated)

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

combined_cluster_plots <- plot1 + plot2

ggsave("snRNAseq_clusters.png", plot = plot1, width = 12, height = 6, dpi = 300)
ggsave("snATACseq_clusters.png", plot = plot2, width = 12, height = 6, dpi = 300)
ggsave("combined_cluster_plots.png", plot = combined_cluster_plots, width = 12, height = 6, dpi = 300)

VlnPlot(hm_integrated_ALS_PN, 'prediction.score.max', group.by = 'predicted.id')

pred_score_max_plot <- VlnPlot(hm_integrated_ALS_PN, 'prediction.score.max', group.by = 'predicted.id')
ggsave("pred_score_max_plot.png", plot = pred_score_max_plot, width = 12, height = 6, dpi = 300)


# # Identify the metadata columns that start with "prediction.score."
# metadata_attributes <- colnames(hm_integrated_ALS_PN[[]])
# prediction_score_attributes <- grep("^prediction.score.", metadata_attributes, value = TRUE)
# prediction_score_attributes <- setdiff(prediction_score_attributes, "prediction.score.max")
# 
# # Extract the prediction score attributes for these cells
# predicted_platelets <- which(hm_integrated_ALS_PN$predicted.id == "Platelet")
# platelet_scores <- hm_integrated_ALS_PN[[]][predicted_platelets, prediction_score_attributes]
# 
# # Order the columns by their average values in descending order
# ordered_columns <- names(sort(colMeans(platelet_scores, na.rm = TRUE), decreasing = TRUE))
# ordered_platelet_scores_df <- platelet_scores[, ordered_columns]
# 
# print(ordered_platelet_scores_df)


predicted_id_counts <- table(hm_integrated_ALS_PN$predicted.id)

# Identify the predicted.id values that have more than 20 cells
major_predicted_ids <- names(predicted_id_counts[predicted_id_counts > 20])
hm_integrated_ALS_PN <- hm_integrated_ALS_PN[, hm_integrated_ALS_PN$predicted.id %in% major_predicted_ids]


# change cell identities to the per-cell predicted labels
Idents(hm_integrated_ALS_PN) <- hm_integrated_ALS_PN$predicted.id


# # replace each cluster label with its most likely predicted label
# for(i in levels(pbmc)) {
#   cells_to_reid <- WhichCells(pbmc, idents = i)
#   newid <- names(which.max(table(pbmc$predicted.id[cells_to_reid])))
#   Idents(pbmc, cells = cells_to_reid) <- newid
# }


## Find differentially accessible peaks between cell types


# if (!requireNamespace("remotes", quietly = TRUE))
#   install.packages('remotes')
# remotes::install_github('immunogenomics/presto')

DefaultAssay(hm_integrated_ALS_PN) <- 'peaks'


#AFTER THIS POINT, I AM ACTUALLY SWITCHING TO RUNNING THE CODE THAT I TAILORED ON THE LINUX MACHINE, AS I THINK
#IT GIVES US MORE MEANINGFUL AND INTERPRETABLE RESULTS WITH OUR DATA


hm_integrated_ALS_PN_da_peaks_All <- FindAllMarkers(hm_integrated_ALS_PN)
hm_integrated_ALS_PN_da_peaks_All_sliced <- hm_integrated_ALS_PN_da_peaks_All
hm_integrated_ALS_PN_da_peaks_All_sliced <- hm_integrated_ALS_PN_da_peaks_All_sliced |> group_by(cluster) %>% slice_max(n = 2, order_by = avg_log2FC)
#closest_genes_all <- ClosestFeature(hm_integrated_ALS, regions = granges(hm_integrated_ALS_PN_da_peaks_All))


# wilcox is the default option for test.use
hm_integrated_ALS_PN_da_peaks_Astro <- FindMarkers(
  object = hm_integrated_ALS_PN,
  ident.1 = "Astro",
  test.use = 'wilcox',
  min.pct = 0.1
)

#Top DA peaks for Astro
top.da.peak_Astro <- rownames(hm_integrated_ALS_PN_da_peaks_Astro[hm_integrated_ALS_PN_da_peaks_Astro$p_val < 0.005 & hm_integrated_ALS_PN_da_peaks_Astro$pct.1 > 0.2, ])
head(top.da.peak_Astro)


hm_integrated_ALS_PN_da_peaks_Endo <- FindMarkers(
  object = hm_integrated_ALS_PN,
  ident.1 = "Endo",
  test.use = 'wilcox',
  min.pct = 0.1
)

#Top DA peaks for Endo
top.da.peak_Endo <- rownames(hm_integrated_ALS_PN_da_peaks_Endo[hm_integrated_ALS_PN_da_peaks_Endo$p_val < 0.005 & hm_integrated_ALS_PN_da_peaks_Endo$pct.1 > 0.2, ])
head(top.da.peak_Endo)

head(hm_integrated_ALS_PN_da_peaks_Astro)
head(hm_integrated_ALS_PN_da_peaks_Endo)

#Basically translating the top DA peaks into genes
closest_genes_Astro <- ClosestFeature(hm_integrated_ALS_PN, regions = top.da.peak_Astro)
closest_genes_Endo <- ClosestFeature(hm_integrated_ALS_PN, regions = top.da.peak_Endo)

#Annotation(hm_integrated_ALS_PN) <- annotations
#closest_genes_all <- ClosestFeature(hm_integrated_ALS_PN, regions = hm_integrated_ALS_PN_da_peaks_All)

#ChatGPT Correction to solve the ClosestFeature issues:
library(GenomicRanges)
library(tidyr)
library(dplyr)

# If your peak coordinates are row names
df <- hm_integrated_ALS_PN_da_peaks_All %>%
  tibble::rownames_to_column(var = "peak")

# Split "chr-start-end" format into separate columns
regions_df <- df %>%
  separate(peak, into = c("seqnames", "start", "end"), sep = "-", convert = TRUE)

# Create a GRanges object
regions_gr <- GRanges(
  seqnames = regions_df$seqnames,
  ranges = IRanges(start = regions_df$start, end = regions_df$end)
)

# # Add metadata columns (optional but helpful)
# mcols(regions_gr) <- regions_df %>%
#   select(-seqnames, -start, -end)

# Now run ClosestFeature
closest_genes_all <- ClosestFeature(hm_integrated_ALS_PN, regions = regions_gr)
head(closest_genes_all)


library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)

closest_genes_all_GO <- enrichGO(gene = closest_genes_all$gene_id,
                                 keyType = "ENSEMBL",
                                 OrgDb = org.Hs.eg.db,
                                 ont = "BP",
                                 pAdjustMethod = "BH",
                                 pvalueCutoff = 0.05,
                                 qvalueCutoff = 0.05,
                                 readable = TRUE)

closest_genes_all_GO_plot <- barplot(closest_genes_all_GO,showCategory = 20)

ggsave("closest_genes_all_GO_plot.png", plot = closest_genes_all_GO_plot, width = 12, height = 6, dpi = 300)




library(Signac)
library(GenomicRanges)
library(dplyr)

# Ensure cluster is character
hm_integrated_ALS_PN_da_peaks_All$cluster <- as.character(hm_integrated_ALS_PN_da_peaks_All$cluster)

# --- Split into list of data frames by cluster ---
da_peaks_by_celltype <- split(hm_integrated_ALS_PN_da_peaks_All, hm_integrated_ALS_PN_da_peaks_All$cluster)

# --- For each cell type, parse region strings and run ClosestFeature ---
closest_genes_by_celltype <- lapply(da_peaks_by_celltype, function(df) {
  # Parse peak strings to GRanges
  coords <- strsplit(df$gene, "-")
  coords_df <- do.call(rbind, coords) %>% as.data.frame()
  colnames(coords_df) <- c("chr", "start", "end")
  coords_df$start <- as.numeric(as.character(coords_df$start))
  coords_df$end <- as.numeric(as.character(coords_df$end))
  
  gr <- GRanges(
    seqnames = coords_df$chr,
    ranges = IRanges(start = coords_df$start, end = coords_df$end)
  )
  
  # Run ClosestFeature
  ClosestFeature(hm_integrated_ALS_PN, regions = gr)
})

#Example for one cell type
head(closest_genes_by_celltype$Astro)  # example for Astrocytes


#Running GO for Astrocytes:
# Step 1: Extract the Astro genes
astro_genes_df <- closest_genes_by_celltype[["Astro"]]

# Step 2: Get unique ENSEMBL gene IDs
#Looking back, I'm not sure if this is a good step to take (reducing things down to just the unique IDs)
#However, it doesn't run if you try to run the non-unique object, so maybe it is necessary.
astro_gene_ids <- unique(astro_genes_df$gene_id)

# Step 3: Run GO enrichment
closest_genes_astro_GO <- enrichGO(
  gene = astro_gene_ids,
  keyType = "ENSEMBL",
  OrgDb = org.Hs.eg.db,
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05,
  readable = TRUE
)

# Step 4: Create barplot
closest_genes_astro_GO_plot <- barplot(closest_genes_astro_GO, showCategory = 20)

# Step 5: Save plot
ggsave("closest_genes_astro_GO_plot.png",
       plot = closest_genes_astro_GO_plot,
       width = 12, height = 10, dpi = 300)

closest_genes_astro_GO_plot


#Same thing for Endo
#Running GO for Endothelial cells:
# Step 1: Extract the Endo genes
endo_genes_df <- closest_genes_by_celltype[["Endo"]]

# Step 2: Get unique ENSEMBL gene IDs
endo_gene_ids <- unique(endo_genes_df$gene_id)

# Step 3: Run GO enrichment
closest_genes_endo_GO <- enrichGO(
  gene = endo_gene_ids,
  keyType = "ENSEMBL",
  OrgDb = org.Hs.eg.db,
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05,
  readable = TRUE
)

# Step 4: Create barplot
closest_genes_endo_GO_plot <- barplot(closest_genes_endo_GO, showCategory = 20)

# Step 5: Save plot
ggsave("closest_genes_endo_GO_plot.png",
       plot = closest_genes_endo_GO_plot,
       width = 12, height = 10, dpi = 300)

regions_gr

#Coverage Plot for GFAP

hm_integrated_ALS_PN <- SortIdents(hm_integrated_ALS_PN)

# find DA peaks overlapping gene of interest
#This first line is extra, provided by ChatGPT to match earlier objects:
astro_regions <- astro_genes_df$query_region
regions_highlight <- subsetByOverlaps(StringToGRanges(astro_regions), LookupGeneCoords(hm_integrated_ALS_PN, "GFAP"))

#The following plot code won't run if Azimuth is installed
GFAP_Coverage_Plot <- CoveragePlot(
  object = hm_integrated_ALS_PN,
  region = "GFAP",
  region.highlight = regions_highlight,
  extend.upstream = 1000,
  extend.downstream = 1000
)


ggsave("GFAP_Coverage_Plot.png",
       plot = GFAP_Coverage_Plot,
       width = 12, height = 10, dpi = 300)
GFAP_Coverage_Plot


# 
# #Setting up the interactive coverage viewer
# Coverage_Browser <- CoverageBrowser(object = hm_integrated_ALS_PN,
#                                     region = "GFAP")

sessionInfo()
