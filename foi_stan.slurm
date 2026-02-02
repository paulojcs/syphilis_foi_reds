#!/bin/bash
#SBATCH --job-name=syphilis_foi_stan
#SBATCH --partition=SP2
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=48:00:00
#SBATCH --output=foi_stan_%j.out
#SBATCH --error=foi_stan_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=your_email@usp.br

# =============================================================================
# SYPHILIS FOI ANALYSIS - CUSTOM STAN MODELS
# =============================================================================
# Submit with: sbatch foi_stan.slurm
# =============================================================================

echo "=========================================="
echo "Job started: $(date)"
echo "Job ID: $SLURM_JOB_ID"
echo "Node: $SLURMD_NODENAME"
echo "CPUs: $SLURM_CPUS_PER_TASK"
echo "=========================================="

# Load modules
module purge
module load R/4.5.1

# OpenSSL workaround (if needed)
export PATH=/temporario2/apps/gnu/openssl-3.5.1/bin:$PATH
export LD_LIBRARY_PATH=/temporario2/apps/gnu/openssl-3.5.1/lib64:$LD_LIBRARY_PATH

# R library path
export R_LIBS_USER=/scratch/7631403/R_packages

# Navigate to working directory
cd /scratch/7631403/reds_data

# Create Stan models directory and copy models
mkdir -p stan_models
cp /scratch/7631403/reds_data/stan_models/*.stan stan_models/ 2>/dev/null || echo "Stan files already in place"

# Verify Stan model files exist
echo ""
echo "Checking Stan model files..."
for f in foi_constant.stan foi_piecewise.stan foi_random_walk.stan foi_random_walk_v2.stan; do
    if [ -f "stan_models/$f" ]; then
        echo "  ✓ $f"
    else
        echo "  ✗ MISSING: $f"
    fi
done

# Run the analysis
echo ""
echo "Starting R analysis..."
Rscript syphilis_foi_stan_analysis.R

# Check completion
if [ $? -eq 0 ]; then
    echo ""
    echo "=========================================="
    echo "Job completed successfully: $(date)"
    echo "=========================================="
else
    echo ""
    echo "=========================================="
    echo "Job failed with error code: $?"
    echo "=========================================="
fi
