import sys
import os
import argparse
import multiprocessing as mp
import importlib.resources

import random
import numpy as np

def parse_args():
    description = "SpatialNN."
    parser = argparse.ArgumentParser(description=description)

    parser.add_argument('-u', '--hidden_units_spatial', type=int, required=False, default=50, help="number of hidden units in the layer")
    parser.add_argument('-o', '--optimizer', type=str, required=False, default="adam", help="the optimizer to use for fitting the neural network")
    parser.add_argument('-p', '--partition', type=str, required=False, help="the portion of the slide to anlayze")
    parser.add_argument('-s', '--seed', type=int, required=False, default=random.randrange(10000), help="Set random seed for reproducibility")


    """
    # some other kind of options copied-pasted from elsewhere we can use
    parser.add_argument("-p","--purityfile", type=str, required=True, help="File with purity of each sample (TSV file in two columns`SAMPLE PURITY`)")
    parser.add_argument("--betabinomial", required=False, default=False, action='store_true', help="Use betabinomial likelihood to cluster mutations (default: binomial)")
    # snpfile and segfile are conditionally required if --betabinomial specified
    parser.add_argument("-i","--snpfile", type=str, required='--betabinomial' in sys.argv, default=None, help="File with precisions for betabinomial fit (default: binomial likelihood)")
    """

    """if not os.path.isfile(args.INPUT):
        raise ValueError("INPUT file does not exist!")"""
    """
    args = parser.parse_args()
    return {
        "hiddenS" : args.hiddenS,
        "optimizer" : args.optimizer,
        "partition" : args.partition,
        "seed" : args.seed
    }"""
    return parser.parse_args()

args = parse_args()




