#!/bin/bash
#SBATCH --job-name=foi_simulation
#SBATCH --partition=SP2
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=24:00:00
#SBATCH --output=foi_simulation_%j.out
#SBATCH --error=foi_simulation_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=your_email@usp.br

# =============================================================================
# FOI SIMULATION STUDY - RANDOM WALK MODEL VALIDATION
# =============================================================================
# Submit with: sbatch foi_simulation.slurm
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

# Verify Stan model file exists
echo ""
echo "Checking Stan model file..."
if [ -f "stan_models/foi_random_walk_v2.stan" ]; then
    echo "  OK: foi_random_walk_v2.stan"
else
    echo "  MISSING: stan_models/foi_random_walk_v2.stan"
    exit 1
fi

# Run the simulation study
echo ""
echo "Starting simulation study..."
Rscript foi_simulation_study.R

# Check completion
if [ $? -eq 0 ]; then
    echo ""
    echo "=========================================="
    echo "Simulation completed successfully: $(date)"
    echo "=========================================="
    echo ""
    echo "Output files:"
    ls -la foi_simulation_*.png foi_simulation_*.csv foi_simulation_*.txt 2>/dev/null
else
    echo ""
    echo "=========================================="
    echo "Simulation failed with error code: $?"
    echo "=========================================="
fi
