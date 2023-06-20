import numpy as np
from sklearn.preprocessing import normalize
import matplotlib.pyplot as plt
from scipy.stats import linregress
from matplotlib.collections import LineCollection

layer_boundary_label='Layer boundaries \n from expression or \n cell type proportion'

# example of ct_colors: {'Oligodendrocytes': 'C6', 'Granule': 'mediumseagreen', 'Purkinje': 'red', 'Bergmann': 'C4', 'MLI1': 'gold', 'MLI2': 'goldenrod',  'Astrocytes': 'C0', 'Golgi': 'C9', 'Fibroblast': 'C5'}

# example of ct_pseudocounts: {3: 1} -- ie layer 4 needs CT pseudocount of 1

# def plot_ct_props(binning_output, ct_colors=None, color_palette=plt.cm.tab20, 
#                   ct_pseudocounts=None, linewidth=8, figsize=(15,6), exclude_ct=[],
#                  layer_ct_threshold=0.6, ticksize=25, layer_boundary_label=layer_boundary_label,
#                  include_lgd=True, lgd_fontsize=25, lgd_bbox=(1.75,1), output_LLR_pvals=True):
#     unique_binned_depths=binning_output['unique_binned_depths']
#     binned_labels=binning_output['binned_labels']
#     ct_count_mat=binning_output['binned_cell_type_mat'].T # len(unique_cell_types) x binned_labels
#     cell_type_names=binning_output['cell_type_names']
    
#     L=len(np.unique(binned_labels))
    
#     left_bps=[]
#     right_bps=[]

#     for i in range(len(binned_labels)-1):
#         if binned_labels[i] != binned_labels[i+1]:
#             left_bps.append(unique_binned_depths[i])
#             right_bps.append(unique_binned_depths[i+1])
    
#     pc_mat=np.zeros( ct_count_mat.shape )
#     if ct_pseudocounts is None:
#         ct_pseudocounts={}
    
#     for l in ct_pseudocounts:
#         pc=ct_pseudocounts[l]
#         pc_mat[:,np.where(binned_labels==l)[0]]=pc
#     ct_count_prop=normalize(ct_count_mat+pc_mat,axis=0,norm='l1')

#     fig,ax=plt.subplots(figsize=figsize)

#     # get biggest cell types in layer
#     layer_ct_markers={}
#     belayer_labels=binning_output['belayer_labels'] # unbinned labels
#     cell_type_mat=binning_output['cell_type_mat'] # unbinned cell types

#     for t in range(L):
#         pts_t=np.where(belayer_labels==t)[0]
#         ct_counts_t=np.sum(cell_type_mat[pts_t,:],0)

#         argsort_cts=np.argsort(ct_counts_t)[::-1]
#         i=0
#         while np.sum(ct_counts_t[argsort_cts[:i]]) / np.sum(ct_counts_t) < layer_ct_threshold:
#             i+=1
#         layer_t_cts=cell_type_names[argsort_cts[:i]]
#         layer_ct_markers[t]=cell_type_names[argsort_cts[:i]]
    
    
#     if ct_colors is not None:
#         ct_plotted=[]
#         for t in layer_ct_markers:
#             ct_list=layer_ct_markers[t]
#             pts_t=np.where(binned_labels==t)[0]
#             for ct in ct_list:
#                 if ct in exclude_ct:
#                     continue
#                 ct_ind=np.where(cell_type_names==ct)[0][0]
#                 if ct not in ct_plotted:
#                     ax.plot(unique_binned_depths[pts_t], ct_count_prop[ct_ind,pts_t], label=ct, color=ct_colors[ct], linewidth=linewidth)
#                 else:
#                     ax.plot(unique_binned_depths[pts_t], ct_count_prop[ct_ind,pts_t], color=ct_colors[ct], linewidth=linewidth)
#                 ct_plotted.append(ct)
#     else:
#         colors_list = iter([color_palette(i) for i in range(20)])
#         ct_color_assignment={}
        
#         for t in layer_ct_markers:
#             ct_list=layer_ct_markers[t]
#             pts_t=np.where(binned_labels==t)[0]
#             for ct in ct_list:
#                 if ct in exclude_ct:
#                     continue
#                 ct_ind=np.where(cell_type_names==ct)[0][0]
#                 if ct not in ct_color_assignment:
#                     ct_color=next(colors_list)
#                     ax.plot(unique_binned_depths[pts_t], ct_count_prop[ct_ind,pts_t], label=ct, color=ct_color, linewidth=linewidth)
#                     ct_color_assignment[ct]=ct_color
#                 else:
#                     ax.plot(unique_binned_depths[pts_t], ct_count_prop[ct_ind,pts_t], color=ct_color_assignment[ct], linewidth=linewidth)
        
        
#     # AXES tick size
#     ax.tick_params(axis='both', which='major', labelsize=ticksize)

#     for i in range(L-1):
#         if i==0:
#             plt.axvline((left_bps[i]+right_bps[i])*0.5, color='black', ls='--', linewidth=3, label=layer_boundary_label)
#         else:
#             plt.axvline((left_bps[i]+right_bps[i])*0.5, color='black', ls='--', linewidth=3)

#     if include_lgd:
#         plt.legend(fontsize=lgd_fontsize, bbox_to_anchor=lgd_bbox)
                            

#     ##################

#     if output_LLR_pvals:
#         for t in np.unique(binned_labels):
#             pts_t=np.where(binned_labels==t)[0]
#             # pts_t=pts_t[3:-3]
#             for ct in layer_ct_markers[int(t)]:
#                 ct_ind=np.where(cell_type_names==ct)[0][0]
#                 ct_prop_in_layer_t=ct_count_prop[ct_ind,pts_t]
#                 # print(ct_prop_in_layer_t)

#                 _,_,_,pv,_=linregress(unique_binned_depths[pts_t], y=ct_prop_in_layer_t)
#                 print(f'layer {int(t)} cell type: {ct}, p-val: {pv}')
#             print('')

############################################################################################################

def get_layer_cts(binning_output, layer_ct_threshold):
    layer_ct_markers={}
    belayer_labels=binning_output['belayer_labels'] # unbinned labels
    cell_type_mat=binning_output['cell_type_mat'] # unbinned cell types
    cell_type_names=binning_output['cell_type_names']
    L=len(np.unique(belayer_labels))

    for t in range(L):
        pts_t=np.where(belayer_labels==t)[0]
        ct_counts_t=np.sum(cell_type_mat[pts_t,:],0)

        argsort_cts=np.argsort(ct_counts_t)[::-1]
        i=0
        while np.sum(ct_counts_t[argsort_cts[:i]]) / np.sum(ct_counts_t) < layer_ct_threshold:
            i+=1
        layer_t_cts=cell_type_names[argsort_cts[:i]]
        layer_ct_markers[t]=cell_type_names[argsort_cts[:i]]
    return layer_ct_markers
            
def plot_ct_props2(binning_output, ct_colors=None, color_palette=plt.cm.tab20, 
                   ct_pseudocounts=None, linewidth=8, figsize=(15,6),
                   layer_ct_threshold=0.6, ticksize=25, layer_boundary_label=layer_boundary_label, exclude_ct=[],
                   alpha1=1, alpha2=0.2,
                   include_lgd=True, lgd_fontsize=25, lgd_bbox=(1.75,1), output_LLR_pvals=True):
    unique_binned_depths=binning_output['unique_binned_depths']
    binned_labels=binning_output['binned_labels']
    ct_count_mat=binning_output['binned_cell_type_mat'].T # len(unique_cell_types) x binned_labels
    cell_type_names=binning_output['cell_type_names']
    
    L=len(np.unique(binned_labels))
    
    left_bps=[]
    right_bps=[]

    for i in range(len(binned_labels)-1):
        if binned_labels[i] != binned_labels[i+1]:
            left_bps.append(unique_binned_depths[i])
            right_bps.append(unique_binned_depths[i+1])
    
    pc_mat=np.zeros( ct_count_mat.shape )
    if ct_pseudocounts is None:
        ct_pseudocounts={}
    
    for l in ct_pseudocounts:
        pc=ct_pseudocounts[l]
        pc_mat[:,np.where(binned_labels==l)[0]]=pc
    ct_count_prop=normalize(ct_count_mat+pc_mat,axis=0,norm='l1')

    fig,ax=plt.subplots(figsize=figsize)

    # get biggest cell types in layer
    layer_ct_markers=get_layer_cts(binning_output, layer_ct_threshold)
        
    if ct_colors is None:
        ct_colors_list=[color_palette(i) for i in range(20)]
    c=0
    # ct_colors_used=[]
    for i,ct in enumerate(cell_type_names):
        if ct in exclude_ct:
            continue
        alphas=np.ones(len(unique_binned_depths))*alpha2
        for l in range(L):
            if ct in layer_ct_markers[l]:
                pts_l=np.where(binned_labels==l)[0]
                alphas[pts_l]=alpha1
        for s in range(len(alphas)-1):
            if alphas[s]==alpha1 and alphas[s+1] < alpha1:
                alphas[s]=alpha2
                          
        if alpha1 in alphas:
            x,y=unique_binned_depths, ct_count_prop[i,:]
            points = np.vstack((x, y)).T.reshape(-1, 1, 2)
            segments = np.hstack((points[:-1], points[1:]))
            if ct_colors is None:
                lc = LineCollection(segments, alpha=alphas, color=ct_colors_list[c], lw=linewidth, label=ct)
            else:
                lc = LineCollection(segments, alpha=alphas, color=ct_colors[ct], lw=linewidth, label=ct)
            line = ax.add_collection(lc)
            c+=1
            # ct_colors_used.append(ct_colors[i])
        # ax.plot(unique_binned_depths, ct_count_prop[i,:], label=ct, color=ct_colors[i], linewidth=linewidth)
        
    # AXES tick size
    ax.tick_params(axis='both', which='major', labelsize=ticksize)

    for i in range(L-1):
        if i==0:
            plt.axvline((left_bps[i]+right_bps[i])*0.5, color='black', ls='--', linewidth=3, label=layer_boundary_label)
        else:
            plt.axvline((left_bps[i]+right_bps[i])*0.5, color='black', ls='--', linewidth=3)

    if include_lgd:
        lgd=plt.legend(fontsize=lgd_fontsize, bbox_to_anchor=lgd_bbox, labelcolor='linecolor')
    
        for lh in lgd.legendHandles: 
            lh.set_alpha(1)
        print(list(lgd.get_texts()))
        for text in lgd.get_texts():
            text.set_alpha(1)
        print(list(lgd.get_texts()))

    ##################

    if output_LLR_pvals:
        for t in np.unique(binned_labels):
            pts_t=np.where(binned_labels==t)[0]
            # pts_t=pts_t[3:-3]
            for ct in layer_ct_markers[int(t)]:
                ct_ind=np.where(cell_type_names==ct)[0][0]
                ct_prop_in_layer_t=ct_count_prop[ct_ind,pts_t]
                # print(ct_prop_in_layer_t)

                _,_,_,pv,_=linregress(unique_binned_depths[pts_t], y=ct_prop_in_layer_t)
                print(f'layer {int(t)} cell type: {ct}, p-val: {pv}')
            print('')
