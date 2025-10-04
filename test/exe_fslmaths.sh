#!/bin/bash
#SBATCH -J test_fsl
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH -p q_fat_c

module load fsl/6.3.0


atlas_dir=/ibmgpfs/cuizaixu_lab/xuxiaoyu/ABCD/processed/qsiPrep/baselineYear1Arm1/SIEMENS/site5/derivatives/qsiprep/sub-NDARINVHAP0JZTR/qsirecon/sub-NDARINVHAP0JZTR/ses-baselineYear1Arm1/dwi
atlas_path=${atlas_dir}/sub-NDARINVHAP0JZTR_ses-baselineYear1Arm1_space-T1w_desc-preproc_desc-schaefer400_atlas.nii.gz
result_path=/ibmgpfs/cuizaixu_lab/xuhaoshu/ADHD_SC_deviation/data/test/extract_mask/fsl
mkdir -p ${result_path}/schaefer376
mkdir -p ${result_path}/Yeo17_sep

for i in `cat /ibmgpfs/cuizaixu_lab/xuxiaoyu/ABCD/batch_code/SIEMENS/dMRI/schaefer376_index_Yeo17.csv`
do
index=$(echo $i | awk -F "," '{print $2}')
echo index${index}
fslmaths $atlas_path -thr $index -uthr $index -bin ${result_path}/schaefer376/schaeferDelLM_${index}.nii.gz
done

for k in $( seq 1 15 )
do
fslmaths $atlas_path -thr 500 -bin ${result_path}/Yeo17_sep/Yeo17_${k}.nii.gz

done

# 初始化后检查是否为空
fslstats ${result_path}/Yeo17_sep/Yeo17_1.nii.gz -V
# 输出应为 "0 0"（非零体素数=0，总体积=0）

dos2unix /ibmgpfs/cuizaixu_lab/xuxiaoyu/ABCD/batch_code/SIEMENS/dMRI/schaefer376_index_Yeo17.csv


for i in `cat /ibmgpfs/cuizaixu_lab/xuxiaoyu/ABCD/batch_code/SIEMENS/dMRI/schaefer376_index_Yeo17.csv`
do
index=$(echo $i | awk -F "," '{print $2}')
name=$(echo $i | awk -F "," '{print $16}')

for j in $( seq 1 15 )
do
if [ "$name" -eq $j ]; then

echo Schaefer${index} network${name}
fslmaths ${result_path}/Yeo17_sep/Yeo17_${j}.nii.gz -add ${result_path}/schaefer376/schaeferDelLM_${index}.nii.gz -bin ${result_path}/Yeo17_sep/Yeo17_${j}.nii.gz

fi

done
done

rm -rf ${result_path}/schaefer376/

for j in $( seq 2 15 )
do

fslmaths ${result_path}/Yeo17_sep/Yeo17_${j}.nii.gz -mul $j ${result_path}/Yeo17_sep/Yeo17_${j}.nii.gz
fslmaths ${result_path}/Yeo17_sep/Yeo17_1.nii.gz -add ${result_path}/Yeo17_sep/Yeo17_${j}.nii.gz ${result_path}/Yeo17_sep/Yeo17_1.nii.gz

done

cp ${result_path}/Yeo17_sep/Yeo17_1.nii.gz ${result_path}/Yeo17_schaefer376_merge.nii.gz
