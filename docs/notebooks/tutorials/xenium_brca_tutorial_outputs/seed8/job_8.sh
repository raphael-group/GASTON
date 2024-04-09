#!/bin/bash
#SBATCH -J 8
#SBATCH -o xenium_brca_tutorial_outputs/seed8/out.8
#SBATCH -e xenium_brca_tutorial_outputs/seed8/err.8
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH -t 0-01:00:00
#SBATCH --mem-per-cpu=5000
source /n/fs/ragr-data/users/uchitra/miniconda3/bin/activate base
conda activate gaston-package
gaston -i xenium_tumor_data/coords_mat_bounded.npy -o xenium_tumor_data/glmpca_bounded.npy --epochs 100000 -d xenium_brca_tutorial_outputs --hidden_spatial 20 20 --hidden_expression 20 20 --optimizer adam --seed 8 -c 500
