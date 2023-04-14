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

import seaborn as sns

import anndata
import scanpy

import argparse
import sys



########################################################################
# MAIN
########################################################################

def get_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument('-s', '--sample', type=str, required=True)
    parser.add_argument('-p', '--hiddenS', type=str, required=True)
    parser.add_argument('-a', '--hiddenA', type=str, required=True)
    parser.add_argument('-o', '--optimizer', type=str, required=True)
    parser.add_argument('-t', '--trial', type=str, required=True)
    return parser

# Run script.
def run(args):
    sample=int(args.sample)
    nhS=int(args.hiddenS)
    nhA=int(args.hiddenA)
    optimizer=str(args.optimizer)
    trial=int(args.trial)

    ###########################
    # STEP 1: LOAD DATA
    ###########################

    adata = anndata.read_h5ad("/n/fs/ragr-data/users/congma/Codes/spatial_embedding/data/{}/sample_{}.h5ad".format(sample,sample))

    data=adata.X.T # n_genes x n_layers
    int_data=adata.layers['count'].T
    coords=adata.obsm["X_pos"]

    # load GLM-PCs
    base_folder='/n/fs/ragr-research/projects/network-mutations/spatial_trans/F_glmpca_2d/'
    F_glmpca=np.load(base_folder+'F_glmpca_2d_{}_penalty_10_trial_0.npy'.format(sample))

    N,G=F_glmpca.shape


    # next, z-score normalize features/spots
    S=coords
    A=F_glmpca

    from sklearn import preprocessing

    scaler = preprocessing.StandardScaler().fit(A)
    A_scaled = scaler.transform(A)

    scaler = preprocessing.StandardScaler().fit(S)
    S_scaled = scaler.transform(S)

    #####

    S_torch=torch.tensor(S_scaled,dtype=torch.float32)
    A_torch=torch.tensor(A_scaled,dtype=torch.float32)


    ###########################
    # STEP 2: RUN MODEL
    ###########################

    SAVE_PATH=f'/n/fs/ragr-research/projects/network-mutations/manifold-alignment/NN_intermediate_ALL_v2/'
    SAVE_PATH += f'sample_{sample}_nhS_{nhS}_nhA_{nhA}_optimizer_{optimizer}_trial_{trial}/'

    # first load model
    sys.path.append('../SpatialNN_code')
    from SpatialNN import SpatialNN, train

    mod, loss_list = train(S_torch, A_torch,
            S_hidden_list=[nhS], A_hidden_list=[nhA], epochs=30000,
            checkpoint=500, SAVE_PATH=SAVE_PATH, optim=optimizer, seed=trial)


    torch.save(mod, SAVE_PATH + f'final_model.pt')
    np.save(SAVE_PATH+'loss_list.npy', loss_list)

    torch.save(S_torch, SAVE_PATH+'Storch.pt')
    torch.save(A_torch, SAVE_PATH+'Atorch.pt')

    np.save(SAVE_PATH+'layer_labels.npy', adata.obs['layer'])

    if sample > 151600:
        # also save harmonic coordinates
        adata = anndata.read_h5ad("/n/fs/ragr-data/users/congma/Codes/spatial_embedding/data/{}/sample_{}_noninit.h5ad".format(sample,sample))
        harmonic_coords=adata.obs["harmonic noninit"]
        np.save(SAVE_PATH+'harmonic_coords.npy', harmonic_coords)


if __name__=='__main__':
    run(get_parser().parse_args(sys.argv[1:]))