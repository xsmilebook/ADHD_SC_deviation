#!/bin/bash
#SBATCH -J MRtrix
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=2
#SBATCH -p batch

# PKU6 construct matrix(17*17)
conda activate mrtrix3

subj=$1

datapath=/share/home/xuxiaoyu/PKUADHD/results/qsiprep/${subj}/qsiprep/${subj}/dwi
atlasdir=/share/home/xuxiaoyu/PKUADHD/results/qsiprep/${subj}/qsirecon/${subj}/dwi
atlaspath=${atlasdir}/${subj}_space-T1w_desc-preproc_desc-schaefer400_atlas.nii.gz

mkdir -p ${atlasdir}/Yeo17_sep
mkdir -p ${atlasdir}/schaefer376

for i in `cat /share/home/xuxiaoyu/PKUADHD/code/schaefer376_index_Yeo17.csv`
do
index=$(echo $i | awk -F "," '{print $2}')
echo index${index}
# fslmaths $atlaspath -thr $index -uthr $index -bin ${atlasdir}/schaefer376/schaeferDelLM_${index}.nii.gz
mrcalc $atlaspath $index -eq ${atlasdir}/schaefer376/schaeferDelLM_${index}.nii.gz -datatype bit

done

for k in $( seq 1 15 )
do

# fslmaths $atlasdir/${subj}_space-T1w_desc-preproc_desc-schaefer400_atlas.nii.gz -thr 500 -bin ${atlasdir}/Yeo17_sep/Yeo17_${k}.nii.gz
mrcalc $atlaspath 0 -mult ${atlasdir}/Yeo17_sep/Yeo17_${k}.nii.gz -datatype bit

done


dos2unix /share/home/xuxiaoyu/PKUADHD/code/schaefer376_index_Yeo17.csv

for i in `cat /share/home/xuxiaoyu/PKUADHD/code/schaefer376_index_Yeo17.csv`
do
index=$(echo $i | awk -F "," '{print $2}')
name=$(echo $i | awk -F "," '{print $16}')

for j in $( seq 1 15 )
do
if [ "$name" -eq $j ]; then

echo Schaefer${index} network${name}
# fslmaths ${atlasdir}/Yeo17_sep/Yeo17_${j}.nii.gz -add $atlasdir/schaefer376/schaeferDelLM_${index}.nii.gz -bin ${atlasdir}/Yeo17_sep/Yeo17_${j}.nii.gz
mrcalc ${atlasdir}/Yeo17_sep/Yeo17_${j}.nii.gz \
        ${atlasdir}/schaefer376/schaeferDelLM_${index}.nii.gz \
        -add 0 -gt \
        ${atlasdir}/Yeo17_sep/Yeo17_${j}.nii.gz -datatype bit
fi

done

done

rm -rf ${atlasdir}/schaefer376/

for j in $( seq 2 15 )
do

# fslmaths ${atlasdir}/Yeo17_sep/Yeo17_${j}.nii.gz -mul $j ${atlasdir}/Yeo17_sep/Yeo17_${j}.nii.gz
# fslmaths ${atlasdir}/Yeo17_sep/Yeo17_1.nii.gz -add ${atlasdir}/Yeo17_sep/Yeo17_${j}.nii.gz ${atlasdir}/Yeo17_sep/Yeo17_1.nii.gz
mrcalc ${atlasdir}/Yeo17_sep/Yeo17_${j}.nii.gz $j -mult ${atlasdir}/Yeo17_sep/Yeo17_${j}.nii.gz
mrcalc ${atlasdir}/Yeo17_sep/Yeo17_1.nii.gz \
        ${atlasdir}/Yeo17_sep/Yeo17_${j}.nii.gz \
        -add ${atlasdir}/Yeo17_sep/Yeo17_1.nii.gz
done

cp ${atlasdir}/Yeo17_sep/Yeo17_1.nii.gz ${atlasdir}/Yeo17_schaefer376_merge.nii.gz
rm -rf ${atlasdir}/Yeo17_sep


