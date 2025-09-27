#!/bin/bash
#SBATCH --job-name=bm_cluster_python
#SBATCH --array=1-135
#SBATCH --partition=campus-new
#SBATCH --time=6:00:00
#SBATCH --mem=32G
#SBATCH --output=/fh/fast/setty_m/user/dotto/benchmarkDA_private/SlurmLog/benchmark_cluster_python_%A_%a.out
#SBATCH --error=/fh/fast/setty_m/user/dotto/benchmarkDA_private/SlurmLog/benchmark_cluster_python_%A_%a.err

# Benchmark execution
# Generated automatically

set -e

# Load environment
eval "$(micromamba shell hook --shell bash 2>/dev/null)" || echo "micromamba not available"

# Navigate to project root
cd /fh/fast/setty_m/user/dotto/benchmarkDA_private

# Execute direct benchmark with SLURM array task ID
if [ -n "${MAMBA_EXE}" ]; then
    ${MAMBA_EXE} run -n benchmarkda python bin/direct_benchmark.py \
        --dataset cluster \
        --method_type python \
        --n_dm 10 \
        --balance No \
        --iteration_num 0 \
        --methods meld_default
else
    python bin/direct_benchmark.py \
        --dataset cluster \
        --method_type python \
        --n_dm 10 \
        --balance No \
        --iteration_num 0 \
        --methods meld_default
fi

echo "Benchmark completed for task $SLURM_ARRAY_TASK_ID"
