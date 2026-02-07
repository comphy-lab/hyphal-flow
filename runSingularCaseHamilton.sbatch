#!/bin/bash
#SBATCH --job-name=hyphaSingleMPI
#SBATCH -N 1
#SBATCH -n 128
#SBATCH -c 1
#SBATCH --mem=240G
#SBATCH --time=3-00:00:00
#SBATCH --gres=tmp:100G
#SBATCH -p multi
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=your.email@durham.ac.uk
#SBATCH --output=slurm-single-%j.out
#SBATCH --error=slurm-single-%j.err

# Compatibility entrypoint for the singular naming variant.
SCRIPT_DIR="${SLURM_SUBMIT_DIR:-$(pwd)}"
exec bash "${SCRIPT_DIR}/runSingleCaseHamilton.sbatch" "$@"
