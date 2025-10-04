#!/bin/bash

# 比较 FSL 和 MRtrix3 生成的两个图谱是否完全一致

FSL_RESULT=$1
MRTRIX_RESULT=$2

if [ $# -ne 2 ]; then
    echo "Usage: $0 <fsl_result.nii.gz> <mrtrix_result.nii.gz>"
    exit 1
fi

if [ ! -f "$FSL_RESULT" ]; then
    echo "Error: FSL result file not found: $FSL_RESULT"
    exit 1
fi

if [ ! -f "$MRTRIX_RESULT" ]; then
    echo "Error: MRtrix3 result file not found: $MRTRIX_RESULT"
    exit 1
fi

echo "Comparing:"
echo "  FSL      : $FSL_RESULT"
echo "  MRtrix3  : $MRTRIX_RESULT"
echo

# Step 1: 检查维度是否一致（防止形状不同）
dim1=$(fslinfo "$FSL_RESULT" | grep dim1 | awk '{print $2}')
dim2=$(fslinfo "$FSL_RESULT" | grep dim2 | awk '{print $2}')
dim3=$(fslinfo "$FSL_RESULT" | grep dim3 | awk '{print $2}')

dim1_m=$(fslinfo "$MRTRIX_RESULT" | grep dim1 | awk '{print $2}')
dim2_m=$(fslinfo "$MRTRIX_RESULT" | grep dim2 | awk '{print $2}')
dim3_m=$(fslinfo "$MRTRIX_RESULT" | grep dim3 | awk '{print $2}')

if [ "$dim1" != "$dim1_m" ] || [ "$dim2" != "$dim2_m" ] || [ "$dim3" != "$dim3_m" ]; then
    echo "❌ ERROR: Image dimensions do not match!"
    echo "FSL:     $dim1 x $dim2 x $dim3"
    echo "MRtrix3: $dim1_m x $dim2_m x $dim3_m"
    exit 1
fi

# Step 2: 计算两个图像的差值（绝对值）
diff_img=$(mktemp --suffix=.nii.gz)
fslmaths "$FSL_RESULT" -sub "$MRTRIX_RESULT" -abs "$diff_img"

# Step 3: 检查最大差值
max_diff=$(fslstats "$diff_img" -R | awk '{print $2}')

# Step 4: 判断是否完全一致
tolerance=1e-6  # 允许极小浮点误差（但你的图谱是整数，应为0）

if (( $(echo "$max_diff <= $tolerance" | bc -l) )); then
    echo "✅ SUCCESS: The two images are identical (max difference = $max_diff)"
    
    # 额外验证：非零体素数是否相同（适用于整数标签图）
    vox_fsl=$(fslstats "$FSL_RESULT" -V | awk '{print $1}')
    vox_mrt=$(fslstats "$MRTRIX_RESULT" -V | awk '{print $1}')
    if [ "$vox_fsl" -eq "$vox_mrt" ]; then
        echo "   - Non-zero voxel count matches: $vox_fsl"
    else
        echo "   ⚠️ Warning: Non-zero voxel counts differ! ($vox_fsl vs $vox_mrt)"
    fi
    
    # 检查标签范围是否一致
    range_fsl=$(fslstats "$FSL_RESULT" -R)
    range_mrt=$(fslstats "$MRTRIX_RESULT" -R)
    if [ "$range_fsl" = "$range_mrt" ]; then
        echo "   - Label range matches: $range_fsl"
    else
        echo "   ⚠️ Warning: Label ranges differ!"
        echo "     FSL:     $range_fsl"
        echo "     MRtrix3: $range_mrt"
    fi
else
    echo "❌ FAILURE: Images differ! Maximum absolute difference = $max_diff"
    echo "   (Expected 0 for integer label maps)"
    
    # 显示差异位置统计
    diff_vox=$(fslstats "$diff_img" -V | awk '{print $1}')
    echo "   Number of differing voxels: $diff_vox"
fi

# 清理临时文件
rm -f "$diff_img"