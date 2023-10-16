#!/bin/bash
#SBATCH -J 18
#SBATCH -o cerebellum_tutorial_outputs/seed18/out.18
#SBATCH -e cerebellum_tutorial_outputs/seed18/err.18
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH -t 0-01:00:00
#SBATCH --mem-per-cpu=5000
source /n/fs/ragr-data/users/uchitra/miniconda3/bin/activate base
conda activate gaston-package
gaston -i cerebellum_data/cerebellum_coords_mat.npy -o cerebellum_data/F_glmpca_penalty_10_rep1.npy --epochs 10000 -d cerebellum_tutorial_outputs --hidden_spatial 20 20 --hidden_expression 20 20 --optimizer adam --seed 18 -c 500
