#!/bin/bash
#SBATCH --job-name=nih_2024
#SBATCH --nodes=1
#SBATCH --mail-type=END,FAIL          # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=david.hodgson@lshtm.ac.uk 
#SBATCH --ntasks-per-node=9
#SBATCH --mem=8G              # Request 8 GB of RAM
#SBATCH --time=48:00:00       # Set runtime to 48 hours
#SBATCH --output=nih_2024_h3_%j.log 
#SBATCH --array=1-6

# Load any necessary modules
#module load R/4.1.2 
module load gnu9/9.4.0
module load boost/1.76.0

source activate ~/miniconda3/envs/R

# Change to the directory where your job script is located
~/miniconda3/envs/R/bin/Rscript nih_2024_run.R $SLURM_ARRAY_TASK_ID
