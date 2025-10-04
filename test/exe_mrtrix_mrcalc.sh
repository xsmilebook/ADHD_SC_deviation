#!/bin/bash
#SBATCH -J test_mrtrix_nii
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH -p q_fat_c

module load mrtrix3

atlas_dir=/ibmgpfs/cuizaixu_lab/xuxiaoyu/ABCD/processed/qsiPrep/baselineYear1Arm1/SIEMENS/site5/derivatives/qsiprep/sub-NDARINVHAP0JZTR/qsirecon/sub-NDARINVHAP0JZTR/ses-baselineYear1Arm1/dwi
atlas_path=${atlas_dir}/sub-NDARINVHAP0JZTR_ses-baselineYear1Arm1_space-T1w_desc-preproc_desc-schaefer400_atlas.nii.gz
result_path=/ibmgpfs/cuizaixu_lab/xuhaoshu/ADHD_SC_deviation/data/test/extract_mask/mrtrix
mkdir -p ${result_path}/schaefer376
mkdir -p ${result_path}/Yeo17_sep

# Step 1: Extract Schaefer regions → .nii.gz
for i in `cat /ibmgpfs/cuizaixu_lab/xuxiaoyu/ABCD/batch_code/SIEMENS/dMRI/schaefer376_index_Yeo17.csv`
do
    index=$(echo $i | awk -F "," '{print $2}')
    echo "index${index}"
    mrcalc $atlas_path $index -eq ${result_path}/schaefer376/schaeferDelLM_${index}.nii.gz
done

# Step 2: Initialize empty masks → .nii.gz
for k in $(seq 1 15)
do
    mrcalc $atlas_path 0 -mult ${result_path}/Yeo17_sep/Yeo17_${k}.nii.gz
done

# Step 3: Merge using -max (logical OR)
dos2unix /ibmgpfs/cuizaixu_lab/xuxiaoyu/ABCD/batch_code/SIEMENS/dMRI/schaefer376_index_Yeo17.csv
for i in `cat /ibmgpfs/cuizaixu_lab/xuxiaoyu/ABCD/batch_code/SIEMENS/dMRI/schaefer376_index_Yeo17.csv`; do
    index=$(echo $i | awk -F "," '{print $2}')
    name=$(echo $i | awk -F "," '{print $16}')
    for j in $(seq 1 15); do
        if [ "$name" -eq $j ]; then
            echo "Schaefer${index} network${name}"
            tmp_out=$(mktemp --suffix=.nii.gz)
            mrcalc ${result_path}/Yeo17_sep/Yeo17_${j}.nii.gz \
                   ${result_path}/schaefer376/schaeferDelLM_${index}.nii.gz \
                   -max $tmp_out -force
            mv $tmp_out ${result_path}/Yeo17_sep/Yeo17_${j}.nii.gz
        fi
    done
done

# Step 4: Combine into label map
cp ${result_path}/Yeo17_sep/Yeo17_1.nii.gz ${result_path}/Yeo17_sep/merged.nii.gz

for j in $(seq 2 15)
do
    # 将 Yeo_j 赋标签值 j
    tmp_label=$(mktemp --suffix=.nii.gz)
    mrcalc ${result_path}/Yeo17_sep/Yeo17_${j}.nii.gz $j -mult $tmp_label -force

    # 累加到 merged
    tmp_merged=$(mktemp --suffix=.nii.gz)
    mrcalc ${result_path}/Yeo17_sep/merged.nii.gz $tmp_label -add $tmp_merged -force
    mv $tmp_merged ${result_path}/Yeo17_sep/merged.nii.gz

    # 清理临时文件
    rm -f $tmp_label
done

# 最终结果
cp ${result_path}/Yeo17_sep/merged.nii.gz ${result_path}/Yeo17_schaefer376_merge.nii.gz
