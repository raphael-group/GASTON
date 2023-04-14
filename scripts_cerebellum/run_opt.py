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

    rep=1
    folder='/n/fs/ragr-research/projects/network-mutations/manifold-alignment/slideseq_cerebellum/'

    counts_mat=np.load(folder+f'counts_mat_rep{rep}.npy')
    coords_mat=np.load(folder+f'coords_mat_rep{rep}.npy')
    gene_labels=np.load(folder+f'gene_labels_rep{rep}.npy',allow_pickle=True)
    F_glmpca=np.load(folder+f'F_glmpca_penalty_10_rep{rep}.npy')
    cell_types=np.load(folder+'cell_types_ALL.npy',allow_pickle=True)

    ###########################
    # STEP 2: load/normalize A, S
    ###########################

    if partition!='all':
        inds=np.load(folder+f'{partition}_cerebellum_inds.npy')
        S=coords_mat[inds,:]
        A=F_glmpca[inds,:]
    else:
        S=coords_mat
        A=F_glmpca

    scaler = preprocessing.StandardScaler().fit(A)
    A_scaled = scaler.transform(A)

    scaler = preprocessing.StandardScaler().fit(S)
    S_scaled = scaler.transform(S)

    S_torch=torch.tensor(S_scaled,dtype=torch.float32)
    A_torch=torch.tensor(A_scaled,dtype=torch.float32)


    ###########################
    # STEP 3: RUN MODEL
    ###########################

    SAVE_PATH=f'/n/fs/ragr-research/users/bjarnold/spatial_transcriptomics/test_dir/slideseq_cerebellum/intermediate_NN_v2/'
    SAVE_PATH += f'nhS_{nhS}_optimizer_{optimizer}_partition_{partition}_trial_{trial}/'


    # first load model
    sys.path.append('/n/fs/ragr-research/users/bjarnold/spatial_transcriptomics/SpatialNN')
    from SpatialNN import SpatialNN, train

    mod, loss_list = train(S_torch, A_torch,
            S_hidden_list=[nhS], A_hidden_list=[10], epochs=50000,
            checkpoint=500, SAVE_PATH=SAVE_PATH, optim=optimizer, seed=trial)

    # sm_init=SpatialNN(A_torch.shape[1], [nhS], [10])

    # mod, loss_list = train(sm_init, S_torch, A_torch, 
    #                                 epochs=20000, batch_size=10, checkpoint=100,
    #                                 SAVE_PATH=SAVE_PATH, optim=optimizer, seed=trial, 
    #                                 momentum=0.9, weight_decay=0)


    torch.save(mod, SAVE_PATH + f'final_model.pt')
    np.save(SAVE_PATH+'loss_list.npy', loss_list)

    torch.save(S_torch, SAVE_PATH+'Storch.pt')
    torch.save(A_torch, SAVE_PATH+'Atorch.pt')


if __name__=='__main__':
    run(get_parser().parse_args(sys.argv[1:]))