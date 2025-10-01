#!/bin/bash
#SBATCH -J MRtrix
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=2
#SBATCH --mem-per-cpu 4000
#SBATCH -p q_cn

# PKU6 construct matrix(17*17)
conda activate mrtrix3
module load fsl/6.3.0

subj=$1

datapath=/share/home/xuxiaoyu/PKUADHD/results/qsiprep/${subj}/qsiprep/${subj}/dwi
atlasdir=/share/home/xuxiaoyu/PKUADHD/results/qsiprep/${subj}/qsirecon//${subj}/dwi
atlaspath=${atlasdir}/${subj}_space-T1w_desc-preproc_desc-schaefer400_atlas.nii.gz
atlasdir2=/share/home/xuxiaoyu/PKUADHD/results/qsiprep/${subj}/qsirecon/${subj}/dwi

mkdir -p ${atlasdir}/Yeo17_sep
mkdir -p ${atlasdir2}/schaefer376

for i in `cat /ibmgpfs/cuizaixu_lab/xuxiaoyu/ABCD/batch_code/SIEMENS/dMRI/schaefer376_index_Yeo17.csv`
do
index=$(echo $i | awk -F "," '{print $2}')
echo index${index}
fslmaths $atlaspath -thr $index -uthr $index -bin ${atlasdir2}/schaefer376/schaeferDelLM_${index}.nii.gz

done

for k in $( seq 1 15 )
do

fslmaths $atlasdir/${scanID}_space-T1w_desc-preproc_desc-schaefer400_atlas.nii.gz -thr 500 -bin ${atlasdir}/Yeo17_sep/Yeo17_${k}.nii.gz

done


dos2unix /ibmgpfs/cuizaixu_lab/xuxiaoyu/ABCD/batch_code/SIEMENS/dMRI/schaefer376_index_Yeo17.csv

for i in `cat /ibmgpfs/cuizaixu_lab/xuxiaoyu/ABCD/batch_code/SIEMENS/dMRI/schaefer376_index_Yeo17.csv`
do
index=$(echo $i | awk -F "," '{print $2}')
name=$(echo $i | awk -F "," '{print $16}')

for j in $( seq 1 15 )
do
if [ "$name" -eq $j ]; then

echo Schaefer${index} network${name}
fslmaths ${atlasdir}/Yeo17_sep/Yeo17_${j}.nii.gz -add $atlasdir2/schaefer376/schaeferDelLM_${index}.nii.gz -bin ${atlasdir}/Yeo17_sep/Yeo17_${j}.nii.gz

fi

done

done

rm -rf ${atlasdir2}/schaefer376/

for j in $( seq 2 15 )
do

fslmaths ${atlasdir}/Yeo17_sep/Yeo17_${j}.nii.gz -mul $j ${atlasdir}/Yeo17_sep/Yeo17_${j}.nii.gz
fslmaths ${atlasdir}/Yeo17_sep/Yeo17_1.nii.gz -add ${atlasdir}/Yeo17_sep/Yeo17_${j}.nii.gz ${atlasdir}/Yeo17_sep/Yeo17_1.nii.gz

done

cp ${atlasdir}/Yeo17_sep/Yeo17_1.nii.gz ${atlasdir}/Yeo17_schaefer376_merge.nii.gz
rm -rf ${atlasdir}/Yeo17_sep


