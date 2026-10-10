#!/bin/bash -l
#SBATCH --job-name=absCSS_B2
#SBATCH --partition=cpu_parti
#SBATCH --cpus-per-task=1
#SBATCH --mem=6G
#SBATCH --time=04:00:00
#SBATCH --output=absCSS_B2_%j.out
#SBATCH --error=absCSS_B2_%j.err

module purge
module load anaconda/2022.10
set -euo pipefail
cd "$(dirname "$0")"
export PYTHONPATH="$PWD/code:${PYTHONPATH:-}"

# M5 production grid: all three interfaces have common half-cell phase.
for p in 1 1.25 1.5 1.75 2; do
  for m in 63 126 188; do
    python code/B2_dipole_solver.py --model M5 --p "$p" --m "$m" --T 1 \
      --interface-phase 0.5 --refresh 2 --tol 1e-9 \
      --out "M5_p${p}_m${m}.json"
  done
done
