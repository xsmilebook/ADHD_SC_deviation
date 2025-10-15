#!/bin/bash

TOTAL_JOBS=1200

SBATCH_DIR="/ibmgpfs/cuizaixu_lab/xuhaoshu/ADHD_SC_deviation/data/bootstrap_code"
mkdir -p $SBATCH_DIR

for i in $(seq 1 $TOTAL_JOBS); do
    cat > "$SBATCH_DIR/bootstrap_job_${i}.slurm" << EOF
#!/bin/bash
#SBATCH --job-name=bootstrap_${i}
#SBATCH --output=/ibmgpfs/cuizaixu_lab/xuhaoshu/ADHD_SC_deviation/reports/log/bootstrap_${i}_%j.out
#SBATCH --error=/ibmgpfs/cuizaixu_lab/xuhaoshu/ADHD_SC_deviation/reports/log/bootstrap_${i}_%j.err
#SBATCH --partition=q_cn
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4

module load R/4.2.2

cd /ibmgpfs/cuizaixu_lab/xuhaoshu/ADHD_SC_deviation/src/S2_normativemodeling/

Rscript S2_bootstrap_ABCD_exe.R $i

echo "Bootstrap job $i completed at \$(date)"
EOF

    sbatch "$SBATCH_DIR/bootstrap_job_${i}.slurm"
    echo "Submitted job ${i}"
done

echo "All ${TOTAL_JOBS} jobs submitted!"