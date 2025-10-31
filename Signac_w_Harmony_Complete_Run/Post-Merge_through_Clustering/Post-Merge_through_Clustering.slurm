#!/bin/bash
#SBATCH --mail-user=joseph.l.guidubaldi@vanderbilt.edu
#SBATCH --mail-type=ALL
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=1-00:00:00
#SBATCH --mem=128G
#SBATCH --output=Post-Merge_through_Clustering.slurm.out

cd /home/guidubjl/Signac/results_integrated_ALS_PN/Signac_w_Harmony_Complete_Run/Post-Merge_through_Clustering

module load gcc openmpi r/4.5.0 hdf5

Rscript --no-save 'Post-Merge_through_Clustering.R'
