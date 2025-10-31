#!/bin/bash
#SBATCH --mail-user=joseph.l.guidubaldi@vanderbilt.edu
#SBATCH --mail-type=ALL
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=7-00:00:00
#SBATCH --mem=128G
#SBATCH --output=ACCRE_Signac_w_Harmony_ALS_PN_Complete.slurm.out

cd /home/guidubjl/Signac/results_integrated_ALS_PN/Signac_w_Harmony_Complete_Run

module load gcc openmpi r/4.5.0 hdf5

Rscript --no-save 'ACCRE_Signac_w_Harmony_ALS_PN_Complete.R'
