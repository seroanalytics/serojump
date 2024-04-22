#!/bin/bash
#SBATCH --job-name=transvir_w2
#SBATCH --nodes=1
#SBATCH --mail-type=END,FAIL          # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=david.hodgson@lshtm.ac.uk 
#SBATCH --ntasks=4
#SBATCH --mem=8G              # Request 8 GB of RAM
#SBATCH --time=48:00:00       # Set runtime to 12 hours
#SBATCH --output=transvir_w2_%j.log 

# Load any necessary modules
#module load R/4.1.2 
module load gnu9/9.4.0 
module load boost/1.76.0

source activate ~/miniconda3/envs/R

# Change to the directory where your job script is located
~/miniconda3/envs/R/bin/Rscript transvir_w2_run.R 
