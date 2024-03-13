#!/bin/bash
#SBATCH -J 15
#SBATCH -o motor_cortex_tutorial_outputs/seed15/out.15
#SBATCH -e motor_cortex_tutorial_outputs/seed15/err.15
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH -t 0-01:00:00
#SBATCH --mem-per-cpu=5000
source /n/fs/ragr-data/users/uchitra/miniconda3/bin/activate base
conda activate gaston-package
gaston -i motor_cortex_data/coords.npy -o motor_cortex_data/glmpca.npy --epochs 10000 -d motor_cortex_tutorial_outputs --hidden_spatial 20 20 --hidden_expression 20 20 --optimizer adam --seed 15 -c 500
