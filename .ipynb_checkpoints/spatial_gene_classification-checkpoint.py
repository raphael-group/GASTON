import numpy as np

def get_discont_genes(pw_fit_dict, binning_output, q=0.95, cell_type=None):
    
    if cell_type is None:
        cell_type = 'all_cell_types'
    
    _,_,discont_mat,_=pw_fit_dict[cell_type]
    # print(discont_mat)
    
    discont_q=np.tile( np.quantile(np.abs(discont_mat), q,0), (discont_mat.shape[0],1))
    discont_genes=list( np.where(np.sum(np.abs(discont_mat) > discont_q,1))[0] )
    
    gene_labels_idx=binning_output['gene_labels_idx']
    return gene_labels_idx[discont_genes]

def get_SV_genes(pw_fit_dict, binning_output, q=0.95, cell_type=None):
    
    if cell_type is None:
        cell_type = 'all_cell_types'
    
    slope_mat,_,_,_=pw_fit_dict[cell_type]
    
    slope_q=np.tile( np.quantile(np.abs(slope_mat), q,0), (slope_mat.shape[0],1))
    svgs=list( np.where(np.sum(np.abs(slope_mat) > slope_q,1))[0] )
    
    gene_labels_idx=binning_output['gene_labels_idx']
    
    return gene_labels_idx[svgs]

def get_layer_specific_genes(pw_fit_dict, l, binning_output, q=0.95, cell_type=None):
    
    L=binning_output['L']
    gene_labels_idx=binning_output['gene_labels_idx']
    
    if cell_type is None:
        cell_type = 'all_cell_types'
    
    slope_mat,_,discont_mat,_=pw_fit_dict[cell_type]
    
    slope_mat_l=slope_mat[:,l]
    large_slope_genes=np.where(np.abs(slope_mat_l) > np.quantile(np.abs(slope_mat_l), q,0))[0]
    
    if l>0:
        discont_left=discont_mat[:,l-1]
        large_left_genes=np.where(np.abs(discont_left) > np.quantile(np.abs(discont_left), q,0))[0]
    else:
        large_left_genes=[]
    
    if l < L-1:
        discont_right=discont_mat[:,l]
        large_right_genes=np.where(np.abs(discont_right) > np.quantile(np.abs(discont_right), q,0))[0]
    else:
        large_right_genes=[]
    
    all_genes=[g for g in range(len(gene_labels_idx)) if g in large_slope_genes or g in large_left_genes or g in large_right_genes]
    
    return gene_labels_idx[all_genes]


# EXAMPLE of layer_cts: {0: ['Oligodendrocytes'], 1: ['Granule'], 2: ['Purkinje', 'Bergmann'], 3: ['MLI1', 'MLI2']}
def get_SV_genes_cell_type(pw_fit_dict, layer_cts, binning_output, perc_diff_threshold=0.5, q=0.95):
    
    slope_mat_all=pw_fit_dict['all_cell_types'][0]
    G,L=slope_mat_all.shape
    
    slope_q=np.tile( np.quantile(np.abs(slope_mat_all), q,0), (slope_mat_all.shape[0],1))
    
    # genes whose layer l slope decreases by at least 50%
    svg_ct_dep_per_layer = {l: [] for l in range(L)} # initialize svg_ct_dep_per_layer with empty lists for each layer
    
    for l in range(L):
        svgs_layer_l = np.where((np.abs(slope_mat_all) > slope_q)[:,l])[0]
        # print(l,np.where(svgs_layer_l==1895)[0])
        cts = layer_cts[l]
        pc_tot = np.zeros(len(svgs_layer_l))
        
        for ct in cts:
            slope_mat_ct,_,_,_ = pw_fit_dict[ct]
            orig_slopes = slope_mat_all[svgs_layer_l,l]
            new_slopes = slope_mat_ct[svgs_layer_l,l]
            pc_tot += np.abs( (orig_slopes-new_slopes)/(orig_slopes) )
        
        pc_tot = pc_tot / len(cts)
        svg_ct_dep_per_layer[l] = list(svgs_layer_l[np.where(pc_tot > perc_diff_threshold)[0]])
    
    for l in range(4):
        svg_ct_dep_per_layer[l] = np.array(list(set(svg_ct_dep_per_layer[l])))
        # print(l,len(svg_ct_dep_per_layer[l]))
    
    svg_ct_dep = []
    for l in range(4):
        svg_ct_dep += list(svg_ct_dep_per_layer[l])
    svg_ct_dep = list(set(svg_ct_dep))
    
    # get genes whose layer specific slope does not decrease by at least 50%
    svgs=list( np.where(np.sum(np.abs(slope_mat_all) > slope_q,1))[0] )
    svg_ct_ind = [g for g in svgs if g not in svg_ct_dep]
    
    gene_labels_idx=binning_output['gene_labels_idx']
    
    return gene_labels_idx[svg_ct_dep], gene_labels_idx[svg_ct_ind]
