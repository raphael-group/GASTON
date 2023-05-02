#!/usr/bin/python -tt

"""
"""

import os
import torch
import torch.nn as nn
import torch.nn.functional as F
import torch.utils
import torch.distributions
import random
import numpy as np
from sklearn import preprocessing
from prettytable import PrettyTable

from spatialnn.parse_args import args
from spatialnn.neural_net import train
from spatialnn import dp_related



def load_rescale_input_data(args):
  coords = np.loadtxt(args.input_layer)
  expression = np.loadtxt(args.output_layer)
  assert coords.shape[0] == expression.shape[0], 'Input and output files do not have same number of rows! Some spots are missing or do not have expression PC values!'
  
  scaler = preprocessing.StandardScaler().fit(expression)
  A_scaled = scaler.transform(expression)
  A_torch = torch.tensor(A_scaled,dtype=torch.float32)

  scaler = preprocessing.StandardScaler().fit(coords)
  S_scaled = scaler.transform(coords)
  S_torch = torch.tensor(S_scaled,dtype=torch.float32)

  return S_torch, A_torch

def count_parameters(model):
  table = PrettyTable(["Modules", "Parameters"])
  total_params = 0
  for name, parameter in model.named_parameters():
    if not parameter.requires_grad: continue
    param = parameter.numel()
    table.add_row([name, param])
    total_params+=param
  print(table)
  print(f"Total Trainable Params: {total_params}")
  return total_params

def main():

  os.makedirs(args.output_dir, exist_ok=True) 

  # Load input data and rescale, S_torch and A_torch are torch tensors representing the spatial coordinates and expression data, respectively
  S_torch, A_torch = load_rescale_input_data(args)
  
  # rescale input data
  """scaler = preprocessing.StandardScaler().fit(expression)
  A_scaled = scaler.transform(expression)
  A_torch=torch.tensor(A_scaled,dtype=torch.float32)

  scaler = preprocessing.StandardScaler().fit(coords)
  S_scaled = scaler.transform(coords)
  S_torch=torch.tensor(S_scaled,dtype=torch.float32)"""

  # train neural network
  mod, loss_list = train(S_torch, A_torch,
                          S_hidden_list=[args.hidden_units_spatial], A_hidden_list=[args.hidden_units_expression], 
                          epochs=args.epochs, checkpoint=args.checkpoint, 
                          save_dir=args.output_dir, optim=args.optimizer, seed=args.seed) 
  count_parameters(mod)

  # save final model and losses per training checkpoint
  torch.save(mod, f'{args.output_dir}/final_model.pt')
  np.savetxt(f'{args.output_dir}/loss_list.txt', loss_list)
  with open(f'{args.output_dir}/min_loss.txt', 'w') as f:
    f.write(str(min(loss_list)) + "\n")
  torch.save(S_torch, f'{args.output_dir}/Storch.pt')
  torch.save(A_torch, f'{args.output_dir}/Atorch.pt')


  A = A_torch.numpy()
  S = S_torch.detach().numpy()
  for max_layers in range(2, args.max_num_layers+1):
    belayer_depth, belayer_labels = dp_related.get_depth_labels(mod, A, S, max_layers)
    dp_related.plot_clusters(belayer_labels, A, S, args.output_dir, max_layers, figsize=(6,6))
    dp_related.plot_vector_field(mod, S, belayer_labels, max_layers, args.output_dir)
  dp_related.plot_depth(belayer_depth, S, args.output_dir)


if __name__ == '__main__':
  main()


