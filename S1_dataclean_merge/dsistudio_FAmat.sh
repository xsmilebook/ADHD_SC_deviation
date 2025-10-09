#!/bin/bash
#SBATCH -J DSIstudio
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=2
#SBATCH -p batch
# module load fsl
source ~/.bashrc
source /share/apps/anaconda3/etc/profile.d/conda.sh
conda activate mrtrix3_new
# set path
subj=$1
predwipath=/share/home/xuxiaoyu/PKUADHD/results/qsiprep/${subj}/qsiprep/${subj}/dwi
reconpath=/share/home/xuxiaoyu/PKUADHD/results/qsirecon/${subj}/qsiprep/${subj}/dwi
# step1 flip bvec - removed event parameter
dwigradcheck ${predwipath}/${subj}_space-T1w_desc-preproc_dwi.nii.gz -grad ${predwipath}/*_dwi.b -export_grad_fsl ${predwipath}/${subj}_space-T1w_desc-preproc_dwi_c.bvec ${predwipath}/${subj}_space-T1w_desc-preproc_dwi_c.bval -force -info -scratch /share/home/xuxiaoyu/PKUADHD/wd/dsistudio/${subj}_$(date +%Y%m%d_%H%M%)
# step2 Convert 4D nifti to SRC files - removed event parameter
dsi_studio --action=src --source=${predwipath}/${subj}_space-T1w_desc-preproc_dwi.nii.gz --bval=${predwipath}/${subj}_space-T1w_desc-preproc_dwi_c.bval --bvec=${predwipath}/${subj}_space-T1w_desc-preproc_dwi_c.bvec --overwrite=1 --output=${predwipath}/${subj}_space-T1w_desc-preproc_dwi.src.gz
# Step3 reconstruction - removed event parameter
# Use the GQI model
dsi_studio --action=rec --source=${predwipath}/${subj}_space-T1w_desc-preproc_dwi.src.gz --method=4 --mask=${predwipath}/${subj}_space-T1w_desc-brain_mask.nii.gz --check_btable=1 --align_acpc=0 --output=${predwipath}/
# Step4 tractography - removed event parameter
dsi_studio --action=trk --source=${predwipath}/${subj}_space-T1w_desc-preproc_dwi.src.gz.gqi.1.25.fib.gz --seed_count=10000000 --threshold_index=qa --turning_angle=0 --step_size=0 --min_length=30 --max_length=300 --output=${predwipath}/${subj}_space-T1w_desc-preproc_dwi.tt.gz
# Step5 export metrics - removed event parameter
dsi_studio --action=exp --source=${predwipath}/${subj}_space-T1w_desc-preproc_dwi.src.gz.gqi.1.25.fib.gz --export=qa,iso,dti_fa,rd,ad
# Step6 extract connectivity - removed event parameter
cp $predwipath/${subj}_space-T1w_desc-preproc_dwi.src.gz.gqi.1.25.fib.gz.dti_fa.nii.gz ${predwipath}/native_fa.nii.gz
# qa matrix - removed event parameter
dsi_studio --action=ana --source=${predwipath}/${subj}_space-T1w_desc-preproc_dwi.src.gz.gqi.1.25.fib.gz --tract=${predwipath}/${subj}_space-T1w_desc-preproc_dwi.tt.gz --connectivity=$reconpath/Yeo17_schaefer376_merge.nii.gz --connectivity_type=pass --connectivity_value=qa --connectivity_output=matrix,connectogram,measure
# fa matrix - removed event parameter
dsi_studio --action=ana --source=${predwipath}/${subj}_space-T1w_desc-preproc_dwi.src.gz.gqi.1.25.fib.gz --tract=${predwipath}/${subj}_space-T1w_desc-preproc_dwi.tt.gz --connectivity=$reconpath/Yeo17_schaefer376_merge.nii.gz --other_slices=${predwipath}/native_fa.nii.gz --connectivity_type=pass --connectivity_value=native_fa --connectivity_output=matrix,connectogram,measure