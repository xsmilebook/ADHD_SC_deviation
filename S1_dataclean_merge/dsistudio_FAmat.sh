#!/bin/bash
#SBATCH -J DSIstudio
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=4
#SBATCH -p q_cn

module load dsi-studio/2022.08.03
module load fsl
module load singularity
module load mrtrix3

# set path
subj=$1 # sub-NDARINVWF7C1DEL
siteid=$2 # site8
event=$3 # baselineYear1Arm1

processedpath=/ibmgpfs/cuizaixu_lab/xuxiaoyu/ABCD/processed/qsiPrep/${event}/SIEMENS/$siteid/derivatives/qsiprep
predwipath=$processedpath/$subj/qsiprep/$subj/ses-${event}/dwi
reconpath=$processedpath/$subj/qsirecon/$subj/ses-${event}/dwi
predwipath2=/GPFS/cuizaixu_lab_permanent/xuxiaoyu/ABCD/processed/qsiPrep/${event}/SIEMENS/$siteid/derivatives/qsiprep/$subj/qsiprep/$subj/ses-${event}/dwi

mkdir -p $predwipath2

# step1 flip bvec
dwigradcheck ${predwipath}/${subj}_ses-${event}_space-T1w_desc-preproc_dwi.nii.gz -grad ${predwipath}/*_dwi.b -export_grad_fsl ${predwipath}/${subj}_ses-${event}_space-T1w_desc-preproc_dwi_c.bvec ${predwipath}/${subj}_ses-${event}_space-T1w_desc-preproc_dwi_c.bval -force -info -scratch /GPFS/cuizaixu_lab_permanent/xuxiaoyu/ABCD/processed/qsiPrep/${subj}_${event}_$(date +%Y%m%d_%H%M%)

# step2 Convert 4D nifti to SRC files
module unload mrtrix3
dsi_studio --action=src --source=${predwipath}/${subj}_ses-${event}_space-T1w_desc-preproc_dwi.nii.gz --bval=${predwipath}/${subj}_ses-${event}_space-T1w_desc-preproc_dwi_c.bval --bvec=${predwipath}/${subj}_ses-${event}_space-T1w_desc-preproc_dwi_c.bvec --overwrite=1 --output=${predwipath2}/${subj}_ses-${event}_space-T1w_desc-preproc_dwi.src.gz


# Step3 reconstruction
# Use the GQI model
dsi_studio --action=rec --source=${predwipath2}/${subj}_ses-${event}_space-T1w_desc-preproc_dwi.src.gz --method=4 --mask=${predwipath}/${subj}_ses-${event}_space-T1w_desc-brain_mask.nii.gz --check_btable=1 --align_acpc=0 --output=${predwipath2}/


# Step4 tractography
dsi_studio --action=trk --source=${predwipath2}/${subj}_ses-${event}_space-T1w_desc-preproc_dwi.src.gz.gqi.1.25.fib.gz --seed_count=10000000 --threshold_index=qa --turning_angle=0 --step_size=0 --min_length=30 --max_length=300 --output=${predwipath2}/${subj}_ses-${event}_space-T1w_desc-preproc_dwi.tt.gz

# Step5 export metrics
dsi_studio --action=exp --source=${predwipath2}/${subj}_ses-${event}_space-T1w_desc-preproc_dwi.src.gz.gqi.1.25.fib.gz --export=qa,iso,dti_fa,rd,ad

# Step6 extract connectivity
cp $predwipath2/${subj}_ses-${event}_space-T1w_desc-preproc_dwi.src.gz.gqi.1.25.fib.gz.dti_fa.nii.gz ${predwipath2}/native_fa.nii.gz

# qa matrix
dsi_studio --action=ana --source=${predwipath2}/${subj}_ses-${event}_space-T1w_desc-preproc_dwi.src.gz.gqi.1.25.fib.gz --tract=${predwipath2}/${subj}_ses-${event}_space-T1w_desc-preproc_dwi.tt.gz --connectivity=$reconpath/Yeo17_schaefer376_merge.nii.gz --connectivity_type=pass --connectivity_value=qa --connectivity_output=matrix,connectogram,measure

# fa matrix
dsi_studio --action=ana --source=${predwipath2}/${subj}_ses-${event}_space-T1w_desc-preproc_dwi.src.gz.gqi.1.25.fib.gz --tract=${predwipath2}/${subj}_ses-${event}_space-T1w_desc-preproc_dwi.tt.gz --connectivity=$reconpath/Yeo17_schaefer376_merge.nii.gz --other_slices=${predwipath2}/native_fa.nii.gz --connectivity_type=pass --connectivity_value=native_fa --connectivity_output=matrix,connectogram,measure


