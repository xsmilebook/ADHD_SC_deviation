#!/bin/bash
#SBATCH --job-name=ABCD
#SBATCH --output=/ibmgpfs/cuizaixu_lab/xuhaoshu/ADHD_SC_deviation/reports/log/Normative_Model/NM.out
#SBATCH --error=/ibmgpfs/cuizaixu_lab/xuhaoshu/ADHD_SC_deviation/reports/log/Normative_Model/NM.err
#SBATCH --partition=q_fat_l
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=72

module load R/4.2.2
Rscript /ibmgpfs/cuizaixu_lab/xuhaoshu/ADHD_SC_deviation/src/S2_normativemodeling/S4_constructNM_forDeviation_ABCD.R
