import numpy as np
from sklearn.preprocessing import normalize
import matplotlib.pyplot as plt
from scipy.stats import linregress
from matplotlib.collections import LineCollection

# layer_boundary_label='Layer boundaries \n from expression or \n cell type proportion'
layer_boundary_label='Layer boundaries'

# example of ct_colors: {'Oligodendrocytes': 'C6', 'Granule': 'mediumseagreen', 'Purkinje': 'red', 'Bergmann': 'C4', 'MLI1': 'gold', 'MLI2': 'goldenrod',  'Astrocytes': 'C0', 'Golgi': 'C9', 'Fibroblast': 'C5'}

# example of ct_pseudocounts: {3: 1} -- ie layer 4 needs CT pseudocount of 1

############################################################################################################

def get_layer_cts(binning_output, layer_ct_threshold, exclude_ct=[]):
    layer_ct_markers={}
    belayer_labels=binning_output['belayer_labels'] # unbinned labels
    cell_type_mat=binning_output['cell_type_mat'] # unbinned cell types
    cell_type_names=binning_output['cell_type_names']
    L=len(np.unique(belayer_labels))
    
    ct_names_not_excluded=np.array([ct for ct in cell_type_names if ct not in exclude_ct])
    ct_inds_not_excluded=np.array([i for i,ct in enumerate(cell_type_names) if ct not in exclude_ct])
    ct_mat_not_excluded=cell_type_mat[:, ct_inds_not_excluded]
    
    for t in range(L):
        pts_t=np.where(belayer_labels==t)[0]
        
        ct_counts_t=np.sum(ct_mat_not_excluded[pts_t,:],0)

        argsort_cts=np.argsort(ct_counts_t)[::-1]
        i=0
        while np.sum(ct_counts_t[argsort_cts[:i]]) / np.sum(ct_counts_t) < layer_ct_threshold:
            i+=1
        layer_t_cts=ct_names_not_excluded[argsort_cts[:i]]
        layer_ct_markers[t]=layer_t_cts
    return layer_ct_markers
            
def plot_ct_props2(binning_output, ct_list=None, ct_colors=None, color_palette=plt.cm.tab20, 
                   ct_pseudocounts=None, linewidth=8, figsize=(15,6),
                   layer_ct_threshold=0.6, ticksize=25, layer_boundary_label=layer_boundary_label, 
                   exclude_ct=[],
                   # alpha1=1, alpha2=0.2,
                   width1=8, width2=4,
                   include_lgd=True, lgd_fontsize=25, lgd_bbox=(1.75,1), lgd_frameon=True, output_LLR_pvals=False):
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
            
    if ct_list is None:
        ct_inds=np.array([i for i,ct in enumerate(cell_type_names) if ct not in exclude_ct])
        ct_list=cell_type_names[ct_inds]
    else:
        ct_inds=np.array([np.where(cell_type_names==i)[0][0] for i in ct_list if i not in exclude_ct])
        ct_list=cell_type_names[ct_inds]
    
    pc_mat=np.zeros( ct_count_mat.shape )
    if ct_pseudocounts is None:
        ct_pseudocounts={}
    
    for l in ct_pseudocounts:
        pc=ct_pseudocounts[l]
        pc_mat[:,np.where(binned_labels==l)[0]]=pc
    
    ct_count_mat=(ct_count_mat+pc_mat)[ct_inds]
    ct_count_prop=normalize(ct_count_mat,axis=0,norm='l1')

    fig,ax=plt.subplots(figsize=figsize)

    # get biggest cell types in layer
    layer_ct_markers=get_layer_cts(binning_output, layer_ct_threshold, exclude_ct=exclude_ct)
        
    if ct_colors is None:
        ct_colors_list=[color_palette(i) for i in range(20)]
    c=0
    
    for i,ct in enumerate(ct_list):
        # ct=ct_list[i]
        # alphas=np.ones(len(unique_binned_depths))*alpha2
        widths=np.ones(len(unique_binned_depths))*width2
        for l in range(L):
            if ct in layer_ct_markers[l]:
                pts_l=np.where(binned_labels==l)[0]
                widths[pts_l]=width1
        for s in range(len(widths)-1):
            if widths[s]==width1 and widths[s+1] < width1:
                widths[s]=width2
        
        x,y=unique_binned_depths, ct_count_prop[i,:]
        points = np.vstack((x, y)).T.reshape(-1, 1, 2)
        segments = np.hstack((points[:-1], points[1:]))
        if ct_colors is None:
            lc = LineCollection(segments, alpha=1, color=ct_colors_list[c], lw=widths, label=ct)
        else:
            lc = LineCollection(segments, alpha=1, color=ct_colors[ct], lw=widths, label=ct)
        line = ax.add_collection(lc)
        c+=1
        
    # AXES tick size
    ax.tick_params(axis='both', which='major', labelsize=ticksize)

    for i in range(L-1):
        if i==0:
            plt.axvline((left_bps[i]+right_bps[i])*0.5, color='black', ls='--', linewidth=3, label=layer_boundary_label)
        else:
            plt.axvline((left_bps[i]+right_bps[i])*0.5, color='black', ls='--', linewidth=3)

    if include_lgd:
        lgd=plt.legend(fontsize=lgd_fontsize, bbox_to_anchor=lgd_bbox, labelcolor='linecolor', frameon=lgd_frameon)
    
        for lh in lgd.legendHandles: 
            lh.set_alpha(1)
        for text in lgd.get_texts():
            text.set_alpha(1)
        for legobj in lgd.legendHandles:
            legobj.set_linewidth(width2)

    ##################
    # print(layer_ct_markers)
    if output_LLR_pvals:
        for t in np.unique(binned_labels):
            pts_t=np.where(binned_labels==t)[0]
            # pts_t=pts_t[3:-3]
            for ct in layer_ct_markers[int(t)]:
                ct_ind=np.where(ct_list==ct)[0][0]
                ct_prop_in_layer_t=ct_count_prop[ct_ind,pts_t]
                # print(ct_prop_in_layer_t)

                _,_,_,pv,_=linregress(unique_binned_depths[pts_t], y=ct_prop_in_layer_t)
                print(f'layer {int(t)} cell type: {ct}, p-val: {pv}')
            print('')
