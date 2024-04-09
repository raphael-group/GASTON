#!/bin/bash
#SBATCH -J 10
#SBATCH -o xenium_brca_tutorial_outputs/seed10/out.10
#SBATCH -e xenium_brca_tutorial_outputs/seed10/err.10
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH -t 0-01:00:00
#SBATCH --mem-per-cpu=5000
source /n/fs/ragr-data/users/uchitra/miniconda3/bin/activate base
conda activate gaston-package
gaston -i xenium_tumor_data/coords_mat_bounded.npy -o xenium_tumor_data/glmpca_bounded.npy --epochs 100000 -d xenium_brca_tutorial_outputs --hidden_spatial 20 20 --hidden_expression 20 20 --optimizer adam --seed 10 -c 500
