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

# Step1
for i in `cat /share/home/xuxiaoyu/PKUADHD/code/schaefer376_index_Yeo17.csv`
do
    index=$(echo $i | awk -F "," '{print $2}')
    echo index${index}
    mrcalc $atlaspath $index -eq ${atlasdir}/schaefer376/schaeferDelLM_${index}.nii.gz
done

# Step2: initialize empty mask
for k in $( seq 1 15 )
do
    mrcalc $atlaspath 0 -mult ${atlasdir}/Yeo17_sep/Yeo17_${k}.nii.gz
done

# Step3: 
# dos2unix /share/home/xuxiaoyu/PKUADHD/code/schaefer376_index_Yeo17.csv
for i in `cat /share/home/xuxiaoyu/PKUADHD/code/schaefer376_index_Yeo17.csv`; do
    index=$(echo $i | awk -F "," '{print $2}')
    name=$(echo $i | awk -F "," '{print $16}')

    for j in $( seq 1 15 ); do
        if [ "$name" -eq $j ]; then
                echo "Schaefer${index} network${name}"
                tmp_out=$(mktemp --suffix=.nii.gz)
                mrcalc ${atlasdir}/Yeo17_sep/Yeo17_${j}.nii.gz \
                        ${atlasdir}/schaefer376/schaeferDelLM_${index}.nii.gz \
                        -max $tmp_out -force
                mv $tmp_out ${atlasdir}/Yeo17_sep/Yeo17_${j}.nii.gz
        fi
    done

done

rm -rf ${atlasdir}/schaefer376/

# Step4
cp ${atlasdir}/Yeo17_sep/Yeo17_1.nii.gz ${atlasdir}/Yeo17_sep/merged.nii.gz
for j in $( seq 2 15 )
do
    tmp_label=$(mktemp --suffix=.nii.gz)
    mrcalc ${atlasdir}/Yeo17_sep/Yeo17_${j}.nii.gz $j -mult $tmp_label -force

    tmp_merged=$(mktemp --suffix=.nii.gz)
    mrcalc ${atlasdir}/Yeo17_sep/merged.nii.gz $tmp_label -add $tmp_merged -force
    
    mv $tmp_merged ${atlasdir}/Yeo17_sep/merged.nii.gz

    rm -f $tmp_label
done

cp ${atlasdir}/Yeo17_sep/merged.nii.gz ${atlasdir}/Yeo17_schaefer376_merge.nii.gz
rm -rf ${atlasdir}/Yeo17_sep


