#!/bin/bash -l
#SBATCH --job-name=absCSS_B5
#SBATCH --partition=cpu_parti
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=06:00:00
#SBATCH --output=logs/B5_%j.out
#SBATCH --error=logs/B5_%j.err

module purge
module load anaconda/2022.10
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
mkdir -p logs out
python3 code/B5_css_rh.py profile --kind lobe --model CSS --p 1.5 --m 125 --T 1 --dt-ratio 0.05 --out out/lobe_CSS.json
python3 code/B5_css_rh.py profile --kind lobe --model RH  --p 1.5 --m 125 --T 1 --dt-ratio 0.05 --out out/lobe_RH.json
python3 code/B5_css_rh.py chase --model CSS --p 1.5 --m 125 --T 1.2 --dt-ratio 0.05 --eta 0.001 --cfast 3 --cslow 1 --separation 8 --out out/chase_CSS.json
python3 code/B5_css_rh.py chase --model RH  --p 1.5 --m 125 --T 1.2 --dt-ratio 0.05 --eta 0.001 --cfast 3 --cslow 1 --separation 8 --out out/chase_RH.json
