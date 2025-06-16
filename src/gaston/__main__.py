#!/usr/bin/python -tt

"""
"""

import os,sys
import torch
import torch.nn as nn
import torch.nn.functional as F
import torch.utils
import torch.distributions
import random
import numpy as np
from sklearn import preprocessing
# from prettytable import PrettyTable

# from gaston.parse_args import args

import argparse
from gaston.neural_net import load_rescale_input_data,train
from gaston import dp_related


def main():
    # PARSE ARGS
    description="GASTON: interpretable deep learning model for mapping tissue topography"
    parser = argparse.ArgumentParser(description=description)

    parser.add_argument('-o', '--output_layer', type=str, required=True, help="filename for N x G numpy array of GLM-PC values (column) for each spatial coordinate (row)")
    parser.add_argument('-i', '--input_layer', type=str, required=True, help="filename for N x 2 numpy array of position (x,y) of each spatial location in sample")
    parser.add_argument('-d', '--output_dir', type=str, required=False, default="./", help="The directory to save the output files")
    
    parser.add_argument('-e', '--epochs', type=int, required=False, default=10000, help="number of epochs to train the neural network")
    parser.add_argument('-c', '--checkpoint', type=int, required=False, default=500, help="save model every checkpoint epochs")

    # parser.add_argument('-u', '--hidden_units_spatial', type=int, required=False, default=50, help="number of hidden units encoding spatial coordinates to depth")
    # parser.add_argument('--hidden_layers_spatial', type=int, required=False, default=1, help="number of hidden layers encoding spatial coordinates to depth")
    parser.add_argument('-p', '--hidden_spatial', nargs='+', type=int, required=True, help="architecture of fully connected NN transforming spatial coordinates (x,y) to isodepth")

    # parser.add_argument('-x', '--hidden_units_expression', type=int, required=False, default=10, help="number of hidden units encoding depth to expression")
    # parser.add_argument('--hidden_layers_expression', type=int, required=False, default=1, help="number of hidden layers encoding depth to expression")
    parser.add_argument('-x', '--hidden_expression', nargs='+', type=int, required=True, help="architecture of fully connected NN transforming isodepth to expression GLM-PCs")

    # parser.add_argument('-t', '--positional_encoding', action='store_true', help="positional encoding option")
    parser.add_argument('-b', '--embedding_size', type=int, required=False, default=4, help="positional encoding embedding size")
    parser.add_argument('-g', '--sigma', type=float, required=False, default=0.2, help="positional encoding sigma hyperparameter")

    parser.add_argument('-z', '--optimizer', type=str, required=False, default="adam", help="optimizer for fitting the neural network")
    parser.add_argument('-s', '--seed', type=int, required=False, default=0, help="Set random seed for reproducibility")

    args=parser.parse_args(sys.argv[1:])

    ###### RUN NN
    out_dir_seed=f"{args.output_dir}/seed{args.seed}" # save in rep{seed}
    os.makedirs(out_dir_seed, exist_ok=True) 
    # Load input data and rescale, S_torch and A_torch are torch tensors representing the spatial coordinates and expression data, respectively
    S=np.load(args.input_layer)
    A=np.load(args.output_layer)
    S_torch, A_torch = load_rescale_input_data(S,A)

    # train neural net
    mod, loss_list = train(S_torch, A_torch,
                          S_hidden_list=args.hidden_spatial, A_hidden_list=args.hidden_expression, 
                          epochs=args.epochs, checkpoint=args.checkpoint, 
                          save_dir=out_dir_seed, optim=args.optimizer, seed=args.seed, save_final=True, 
                          sigma=args.sigma) 



if __name__ == '__main__':
    main()




# def count_parameters(model):
#     table = PrettyTable(["Modules", "Parameters"])
#     total_params = 0
#     for name, parameter in model.named_parameters():
#         if not parameter.requires_grad: continue
#             param = parameter.numel()
#             table.add_row([name, param])
#             total_params+=param
#   print(table)
#   print(f"Total Trainable Params: {total_params}")
#   return total_params