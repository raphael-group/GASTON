import numpy as np
import matplotlib.pyplot as plt

from sklearn.preprocessing import normalize

from mpl_toolkits.axes_grid1.inset_locator import mark_inset, inset_axes

import torch
import torch.nn as nn
import torch.nn.functional as F
import torch.utils
import torch.distributions

def plot_clusters(belayer_labels, A, S, fig=None, ax=None, figsize=(5,8), colors=None, color_palette=plt.cm.Dark2, s=20,labels=None):
    
    if fig is None or ax is None:
        fig,ax=plt.subplots(figsize=figsize)

    if colors is None:
        colors=np.array([color_palette(i) for i in range(len(np.unique(belayer_labels)))])

    for i,t in enumerate(np.unique(belayer_labels)):
        pts_t=np.where(belayer_labels==t)[0]
        l=None
        if labels is not None:
            l=labels[i]
        ax.scatter(S[pts_t,0], S[pts_t,1],s=s, color=colors[int(t)], label=l)
    plt.axis('off')
    
def plot_clusters_boundary(belayer_labels, A, S, boundary_locs, belayer_depth=None, 
                           fig=None, ax=None, figsize=(5,8), colors=None, 
                           color_palette=plt.cm.Dark2, color_boundary='lime', alpha_boundary = 1, s=20):
    
    if fig is None or ax is None:
        fig,ax=plt.subplots(figsize=figsize)
    
    if colors is None:
        colors=np.array([color_palette(i) for i in range(len(np.unique(belayer_labels)))])
    
    for t in np.unique(belayer_labels):
        pts_t=np.where(belayer_labels==t)[0]
        ax.scatter(S[pts_t,0],S[pts_t,1],s=s,c=colors[int(t)], alpha=1)

    ax.scatter(S[boundary_locs,0], S[boundary_locs,1],s=s, alpha=alpha_boundary, c=color_boundary)
    
    if belayer_depth is not None:
        ax.tricontour(S[:,0], S[:,1], belayer_depth, 
                      levels=[np.min(belayer_depth[belayer_labels==i]) for i in range(5)[1:]], 
                      linewidths=2, colors='k', linestyles='-')

    plt.axis('off')
    plt.tight_layout()
                           


def plot_vector_field(model, belayer_labels, S, transform_mat=None, figsize=(5,8), colors=None, color_palette=plt.cm.Dark2, normalize_grads=True, scale=10, width=1.5e-3):
    x=torch.tensor(S,requires_grad=True).float()
    grads=torch.autograd.grad(outputs=model.spatial_embedding(x).flatten(),inputs=x, grad_outputs=torch.ones_like(x[:,0]))[0]

    N=S.shape[0]
    L=len(np.unique(belayer_labels))
    
    if transform_mat is None:
        transform_mat=np.eye(2)
    
    S=(transform_mat @ S.T).T
    grads=(transform_mat @ grads.detach().numpy().T).T

    if colors is None:
        colors=np.array([color_palette(i) for i in range(L)])
    spot_colors=colors[[int(t) for t in belayer_labels]]

    fig,ax=plt.subplots(figsize=figsize)

    if normalize_grads:
        grads=normalize(grads,axis=1,norm='l2')

    im1=ax.quiver(S[:,0],S[:,1], grads[:,0], grads[:,1], 
                  scale=scale, scale_units='inches', width=width,color=spot_colors)
    plt.axis('off')
    
def plot_depth(belayer_depth, S, figsize=(5,8), contours=True, contour_levels=4, contour_lw=1, contour_fs=10, colorbar=False,s=20,cbar_fs=10, axis_off=True):
    fig,ax=plt.subplots(figsize=figsize)

    im1=ax.scatter(S[:,0], S[:,1], c=belayer_depth, cmap='Reds', s=s)
    if axis_off:
        plt.axis('off')

    if contours:
        CS=ax.tricontour(S[:,0], S[:,1], belayer_depth, levels=contour_levels, linewidths=contour_lw, colors='k', linestyles='solid')
        ax.clabel(CS, CS.levels, inline=True, fontsize=contour_fs)
    if colorbar:
        cbar=plt.colorbar(im1)
        cbar.ax.tick_params(labelsize=cbar_fs)