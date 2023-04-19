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

from spatialnn.parse_args import args

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
  print(args.seed)
  set_seeds(args.seed)

if __name__ == '__main__':
  main()


