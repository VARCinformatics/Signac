#!/bin/bash
#SBATCH --mail-user=joseph.l.guidubaldi@vanderbilt.edu
#SBATCH --mail-type=ALL
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=1-00:00:00
#SBATCH --mem=64G
#SBATCH --output=Generating_TSSEnrichmnet5.slurm.out

cd /home/guidubjl/Signac/results_integrated_ALS_PN/Signac_IndividualSample_QC/Testing_Different_TSSEnrichment_Cutoffs

module load gcc openmpi r/4.5.0 hdf5

Rscript --no-save 'Generating_TSSEnrichmnet5.R'
