import torch
import numpy as np

import os
import torch
import torch.nn as nn

from gaston.neural_net import get_loss

def process_files(output_folder):
    smallest_loss = np.Inf
    best_model_folder_path = None
    best_mod=None

    for folder_name in os.listdir(output_folder):
        folder_path = os.path.join(output_folder, folder_name)

        if os.path.isdir(folder_path):
            model_path = os.path.join(folder_path, 'final_model.pt')

            if os.path.exists(model_path):
                try:
                    mod =  torch.load(model_path)
                    St=torch.load(os.path.join(folder_path, 'Storch.pt'))
                    At=torch.load(os.path.join(folder_path, 'Atorch.pt'))
                    loss = get_loss(mod,St,At)

                    if loss < smallest_loss:
                        smallest_loss = loss
                        best_model_folder_path = folder_path
                        best_mod=mod
                except Exception as e:
                    raise Exception(f"Error loading model from {model_path}: {str(e)}")
    print(f'best model: {best_model_folder_path}')
    if best_model_folder_path:
        # folder_name = os.path.basename(os.path.dirname(best_model_path))
        storch_path = os.path.join(best_model_folder_path, 'Storch.pt')
        atorch_path = os.path.join(best_model_folder_path, 'Atorch.pt')

        if os.path.exists(storch_path) and os.path.exists(atorch_path):
            A_torch = torch.load(atorch_path)
            S_torch = torch.load(storch_path)
            
            A = A_torch.detach().numpy()
            S = S_torch.detach().numpy()

    else:
        raise Exception("No 'final_model.pt' found in any folder.")

    return best_mod, A, S
