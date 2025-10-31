#!/bin/bash
#SBATCH --mail-user=joseph.l.guidubaldi@vanderbilt.edu
#SBATCH --mail-type=ALL
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=1-00:00:00
#SBATCH --mem=64G
#SBATCH --output=Signac_ALS_PN_IndividualSampleQCPlots.slurm.out

cd /home/guidubjl/Signac/results_integrated_ALS_PN/Signac_IndividualSample_QC/QC_Plots_by_Sample

module load gcc openmpi r/4.5.0 hdf5

Rscript --no-save 'Signac_ALS_PN_IndividualSampleQCPlots.R'
