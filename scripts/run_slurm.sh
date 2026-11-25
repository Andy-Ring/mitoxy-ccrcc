#!/bin/bash
#SBATCH --job-name=mitoxy
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=24:00:00
#SBATCH --output=mitoxy-%j.log
set -euo pipefail
cd "${SLURM_SUBMIT_DIR:-$PWD}"
source .venv/bin/activate
export MITOXY_THREADS="${SLURM_CPUS_PER_TASK:-1}"
export OPENBLAS_NUM_THREADS=1
export OMP_NUM_THREADS=1
if (( $# == 0 )); then set -- figure4; fi
python scripts/reproduce.py --download --stages "$@"
