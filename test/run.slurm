#!/bin/bash
#SBATCH --job-name=job1
#SBATCH --account=b1042
#SBATCH --partition=genomics
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --time=1:00:00
#SBATCH --mem=50G

cd $SLURM_SUBMIT_DIR
source /home/qgn1237/2_software/mambaforge/etc/profile.d/conda.sh
source /home/qgn1237/2_software/mambaforge/etc/profile.d/mamba.sh

mamba activate directclean 

poetry run python ./test.py
