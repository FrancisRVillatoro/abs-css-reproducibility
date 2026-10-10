#!/bin/bash -l
#SBATCH --job-name=absCSS_B3
#SBATCH --partition=cpu_parti
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=04:00:00
#SBATCH --array=0-5%2
set -euo pipefail
module purge
module load anaconda/2022.10
ROOT="$SLURM_SUBMIT_DIR"
cd "$ROOT"
# Six clean epsilon-limit cases p=3/2, L=8, r=1,2, eps=.1,.06,.03
case "$SLURM_ARRAY_TASK_ID" in
0) E=.1; R=1;; 1) E=.1; R=2;; 2) E=.06; R=1;; 3) E=.06; R=2;; 4) E=.03; R=1;; 5) E=.03; R=2;; esac
python code/B3_regularized_solver.py --initial sine --p 1.5 --r "$R" --epsilon "$E" --points 16 --L 8 --T .01 --dt-ratio .1 --out "repro_p15_e${E}_r${R}.json"
