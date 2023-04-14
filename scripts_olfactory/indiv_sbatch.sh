#!/usr/bin/env bash
#SBATCH -N 1                  ## Node count
#SBATCH --ntasks-per-node=1   ## Processors per node
#SBATCH -t 7:59:59            ## Walltime
#SBATCH --mem=100GB            ## Memory
#SBATCH -A raphael            ## partition

scripts_dir=/n/fs/ragr-research/projects/network-mutations/manifold-alignment/olfactory_glmpca/scripts_NN
source /n/fs/ragr-data/users/uchitra/miniconda3/bin/activate base
conda activate belayer2

echo $1
echo $2
echo $3
echo $4

echo $scripts_dir/run_opt.py
python $scripts_dir/run_opt.py \
    -s $1 \
    -o $2 \
    -p $3 \
    -t $4