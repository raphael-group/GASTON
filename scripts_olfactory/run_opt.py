#################################################################
############################ IMPORTS ############################
#################################################################
import torch
import torch.nn as nn
import torch.nn.functional as F
import torch.utils
import torch.distributions

device = 'cuda' if torch.cuda.is_available() else 'cpu'

import numpy as np
import matplotlib.pyplot as plt

from sklearn.metrics import adjusted_rand_score
from sklearn import preprocessing

import seaborn as sns

import anndata
import scanpy

import pandas as pd
import copy

import argparse
import csv
import os
import sys


########################################################################
# MAIN
########################################################################

def get_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument('-s', '--hiddenS', type=str, required=True)
    parser.add_argument('-o', '--optimizer', type=str, required=True)
    parser.add_argument('-p', '--partition', type=str, required=True)
    parser.add_argument('-t', '--trial', type=str, required=True)
    return parser

# Run script.
def run(args):
    nhS=int(args.hiddenS)
    optimizer=str(args.optimizer)
    partition=str(args.partition)
    trial=int(args.trial)

    ###########################
    # STEP 1: LOAD DATA
    ###########################

    folder = '/n/fs/ragr-data/datasets/olfactory_bulb_stereoseq/Dataset1_LiuLongQi_MouseOlfactoryBulb/'

    counts_file=folder+'RNA_counts.tsv'
    coor_file=folder+'position.tsv'

    counts = pd.read_csv(counts_file, sep='\t', index_col=0)
    coor_df = pd.read_csv(coor_file, sep='\t')
    spot_labels=list(coor_df['label'])

    counts_mat=counts.to_numpy().T
    gene_labels=np.array(counts.index)

    coords=coor_df[['x', 'y']].to_numpy()

    N,G=counts_mat.shape

    ########### load GLM-PCA ############

    penalty=10
    glmpca_path=f'/n/fs/ragr-research/projects/network-mutations/manifold-alignment/olfactory_glmpca/F_glmpca_2d_trial_0_penalty_{penalty}.npy'
    F_glmpca=np.load(glmpca_path)

    ###########################
    # STEP 2: CROP TO LOWER BULB
    ###########################

    used_spot_labels = pd.read_csv('/n/fs/ragr-research/projects/network-mutations/manifold-alignment/olfactory_glmpca/used_barcodes.txt')
    used_spot_labels = [int(x[5:]) for x in list(used_spot_labels['Spot_1'])]
    used_spots=np.array( [spot_labels.index(i) for i in used_spot_labels] )

    if partition=='lower':
        used_spots=used_spots[ np.where(coords[used_spots,1] < 9200)[0] ]
    elif partition=='upper':
        used_spots=used_spots[ np.where(coords[used_spots,1] >= 9200)[0] ]

    S=coords[used_spots,:]
    A=F_glmpca[used_spots,:]


    # next, normalize A, S

    scaler = preprocessing.StandardScaler().fit(A)
    A_scaled = scaler.transform(A)

    scaler = preprocessing.StandardScaler().fit(S)
    S_scaled = scaler.transform(S)

    #####

    S_torch=torch.tensor(S_scaled,dtype=torch.float32)
    A_torch=torch.tensor(A_scaled,dtype=torch.float32)


    # S_torch=torch.tensor(S,dtype=torch.float32)
    # A_torch=torch.tensor(A,dtype=torch.float32)


    ###########################
    # STEP 3: RUN MODEL
    ###########################

    # nhS_"$nhS"_optimizer_"$optimizer"_partition_"$partition"_trial_"$trial"

    SAVE_PATH=f'/n/fs/ragr-research/projects/network-mutations/manifold-alignment/olfactory_glmpca/intermediate_NN_v2/'
    SAVE_PATH += f'nhS_{nhS}_optimizer_{optimizer}_partition_{partition}_trial_{trial}/'


    # first load model
    sys.path.append('/n/fs/ragr-research/projects/network-mutations/manifold-alignment/SpatialNN_code')
    from SpatialNN import SpatialNN, train

    mod, loss_list = train(S_torch, A_torch,
            S_hidden_list=[nhS], A_hidden_list=[10], epochs=50000,
            checkpoint=500, SAVE_PATH=SAVE_PATH, optim=optimizer, seed=trial)



    torch.save(mod, SAVE_PATH + f'final_model.pt')
    np.save(SAVE_PATH+'loss_list.npy', loss_list)

    torch.save(S_torch, SAVE_PATH+'Storch.pt')
    torch.save(A_torch, SAVE_PATH+'Atorch.pt')



if __name__=='__main__':
    run(get_parser().parse_args(sys.argv[1:]))