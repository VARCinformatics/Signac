#!/bin/bash
#SBATCH --mail-user=joseph.l.guidubaldi@vanderbilt.edu
#SBATCH --mail-type=ALL
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=2-00:00:00
#SBATCH --mem=64G
#SBATCH --output=TSSEnrichment5_MoreStrict_Cutoff_test.slurm.out

cd /home/guidubjl/Signac/results_integrated_ALS_PN/Signac_IndividualSample_QC/Testing_Different_TSSEnrichment_Cutoffs/TSSEnrichment5_test/TSSEnrichment5_MoreStrict

module load gcc openmpi r/4.5.0 hdf5

Rscript --no-save 'TSSEnrichment5_MoreStrict_Cutoff_test.R'
