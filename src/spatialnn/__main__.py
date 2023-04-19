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

from spatialnn.parse_args import args
from spatialnn.neural_net import train
from spatialnn import dp_related

def set_seeds(seed):
    torch.manual_seed(seed)
    random.seed(seed)
    os.environ['PYTHONHASHSEED'] = str(seed)
    torch.cuda.manual_seed(seed)
    torch.cuda.manual_seed_all(seed) # if you are using multi-GPU.
    np.random.seed(seed) # Numpy module.
    torch.manual_seed(seed)
    torch.backends.cudnn.benchmark = False
    torch.backends.cudnn.deterministic = True

def main():

  set_seeds(args.seed)
  os.makedirs(args.output_dir, exist_ok=True) 

  # Load input data
  coords = np.loadtxt(args.input_layer, delimiter='\t')
  expression = np.loadtxt(args.output_layer, delimiter='\t')
  assert coords.shape[0] == expression.shape[0], 'Input and output files do not have same number of rows! Some spots are missing or do not have expression PC values!'

  # rescale input data
  scaler = preprocessing.StandardScaler().fit(expression)
  A_scaled = scaler.transform(expression)
  A_torch=torch.tensor(A_scaled,dtype=torch.float32)

  scaler = preprocessing.StandardScaler().fit(coords)
  S_scaled = scaler.transform(coords)
  S_torch=torch.tensor(S_scaled,dtype=torch.float32)

  # train neural network
  mod, loss_list = train(S_torch, A_torch,
                          S_hidden_list=[args.hidden_units_spatial], A_hidden_list=[10], 
                          epochs=args.epochs, checkpoint=1000, 
                          save_dir=args.output_dir, optim=args.optimizer, seed=args.seed) 

  A = A_torch.numpy()
  S = S_torch.detach().numpy()
  belayer_depth, belayer_labels = dp_related.get_depth_labels(mod, A, S, args.num_layers)
  dp_related.plot_clusters(belayer_labels, A, S, args.output_dir, figsize=(6,6))

  # save final model and losses per training checkpoint
  torch.save(mod, f'{args.output_dir}/final_model.pt')
  np.savetxt(f'{args.output_dir}/loss_list.npy', loss_list)
  torch.save(S_torch, f'{args.output_dir}/Storch.pt')
  torch.save(A_torch, f'{args.output_dir}/Atorch.pt')

if __name__ == '__main__':
  main()


