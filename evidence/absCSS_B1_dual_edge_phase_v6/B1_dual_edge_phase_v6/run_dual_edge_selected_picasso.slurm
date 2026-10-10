#!/bin/bash -l
#SBATCH --job-name=absCSS_B1_dual
#SBATCH --partition=cpu_parti
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=06:00:00
#SBATCH --array=0-11%12

set -euo pipefail
module purge
module load anaconda/2022.10

CASES=(
  "63 1.75 .75 .25 M1" "126 1.75 .75 .25 M1" "188 1.75 .75 .25 M1" "251 1.75 .75 .25 M1"
  "63 2.0 .75 .25 M1" "126 2.0 .75 .25 M1" "188 2.0 .75 .25 M1" "251 2.0 .75 .25 M1"
  "63 2.0 .75 .25 M5" "126 2.0 .75 .25 M5" "188 2.0 .75 .25 M5" "251 2.0 .75 .25 M5"
)
read -r M P PL PR METHOD <<< "${CASES[$SLURM_ARRAY_TASK_ID]}"
mkdir -p results
python checkpoint_case.py \
  --m "$M" --p "$P" --psiL "$PL" --psiR "$PR" --method "$METHOD" \
  --chunk 1000000 --state "results/${METHOD}_p${P}_L${PL}_R${PR}_m${M}.npz"
