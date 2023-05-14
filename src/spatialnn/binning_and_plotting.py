from scipy.stats import mode
import numpy as np
import matplotlib.pyplot as plt

# counts mat has to be of shape G x N
def bin_data(counts_mat, belayer_labels, belayer_depth, 
              cell_type_df, gene_labels, num_bins=70, 
             idx_kept=None, umi_threshold=500):
    
    if idx_kept is None:
        idx_kept=np.where(np.sum(counts_mat,1) > umi_threshold)[0]
    gene_labels_idx=gene_labels[idx_kept]
    
    pseudo_counts_mat=counts_mat+1
    exposure=np.sum(pseudo_counts_mat,axis=0)
    cmat=pseudo_counts_mat[idx_kept,:]
    
    cell_type_mat=cell_type_df.to_numpy()
    cell_type_names=np.array(cell_type_df.columns)

    G,N=cmat.shape

    # BINNING
    depth_min, depth_max=np.floor(np.min(belayer_depth))-0.5, np.ceil(np.max(belayer_depth))+0.5
    bins=np.linspace(depth_min, depth_max, num=num_bins+1)

    unique_binned_depths=np.array( [0.5*(bins[i]+bins[i+1]) for i in range(len(bins)-1)] )
    binned_depth_inds=np.digitize(belayer_depth, bins)-1 #ie [1,0,3,15,...]
    binned_depths=unique_binned_depths[binned_depth_inds]
    
    # remove bins not used
    unique_binned_depths=np.delete(unique_binned_depths,
                                   [np.where(unique_binned_depths==t)[0][0] for t in unique_binned_depths if t not in binned_depths])

    N_1d=len(unique_binned_depths)
    binned_count=np.zeros( (G, N_1d) )
    binned_exposure=np.zeros( N_1d )
    binned_labels=np.zeros(N_1d)
    binned_cell_type_mat=np.zeros((N_1d, len(cell_type_names)))
    
    binned_count_per_ct={ct: np.zeros( (G, N_1d) ) for ct in cell_type_names}
    binned_exposure_per_ct={ct: np.zeros( N_1d ) for ct in cell_type_names}

    map_1d_bins_to_2d={} # map b -> [list of cells in bin b]
    for ind, b in enumerate(unique_binned_depths):
        bin_pts=np.where(binned_depths==b)[0]

        binned_count[:,ind]=np.sum(cmat[:,bin_pts],axis=1)
        binned_exposure[ind]=np.sum(exposure[bin_pts])
        binned_labels[ind]= int(mode( belayer_labels[bin_pts],keepdims=False).mode)
        binned_cell_type_mat[ind,:] = np.sum( cell_type_mat[bin_pts,:], axis=0)
        map_1d_bins_to_2d[b]=bin_pts
        
        for ct_ind, ct in enumerate(cell_type_names):
            
            ct_spots=np.where(cell_type_mat[:,ct_ind] > 0)[0]
            ct_spots_bin = [t for t in ct_spots if t in bin_pts]
            ct_spots_bin_proportions=cell_type_mat[ct_spots_bin,ct_ind]
            
            if len(ct_spots_bin)>0:
                binned_count_per_ct[ct][:,ind]=np.sum(cmat[:,ct_spots_bin] * np.tile(ct_spots_bin_proportions,(G,1)), axis=1)
                binned_exposure_per_ct[ct][ind]=np.sum(exposure[ct_spots_bin] * ct_spots_bin_proportions)
            
    L=len(np.unique(belayer_labels))
    segs=[np.where(binned_labels==i)[0] for i in range(L)]

    to_return={}
    
    to_return['L']=len(np.unique(belayer_labels))
    to_return['umi_threshold']=umi_threshold
    to_return['belayer_labels']=belayer_labels
    to_return['pseudo_counts_mat_idx']=cmat
    to_return['cell_type_mat']=cell_type_mat
    to_return['cell_type_names']=cell_type_names
    to_return['gene_labels_idx']=gene_labels_idx
    
    to_return['binned_depths']=binned_depths
    to_return['unique_binned_depths']=unique_binned_depths
    to_return['binned_count']=binned_count
    to_return['binned_exposure']=binned_exposure
    to_return['binned_labels']=binned_labels
    to_return['binned_cell_type_mat']=binned_cell_type_mat
    
    to_return['binned_count_per_ct']=binned_count_per_ct
    to_return['binned_exposure_per_ct']=binned_exposure_per_ct
    
    to_return['map_1d_bins_to_2d']=map_1d_bins_to_2d
    to_return['segs']=segs

    return to_return

def plot_gene_pwlinear(gene_name, pw_fit_dict, belayer_labels, belayer_depth, binning_output, 
                       cell_type=None, spot_threshold=0.25, pt_size=10, 
                       colors=None, color_palette=plt.cm.Dark2, linear_fit=True, layer_list=None, ticksize=20, figsize=(7,3)):
    
    # idx_kept=np.where(np.sum(binning_output['counts_mat'],1) > umi_threshold)[0]
    # gene_labels_idx=gene_labels[idx_kept]
    gene_labels_idx=binning_output['gene_labels_idx']
    if gene_name in gene_labels_idx:
        gene=np.where(gene_labels_idx==gene_name)[0]
    else:
        umi_threshold=binning_output['umi_threshold']
        raise ValueError(f'gene does not have UMI count above threshold {umi_threshold}')
    
    unique_binned_depths=binning_output['unique_binned_depths']
    binned_labels=binning_output['binned_labels']
    
    if cell_type is None:
        binned_count=binning_output['binned_count']
        binned_exposure=binning_output['binned_exposure']
        
    else:
        binned_count=binning_output['binned_count_per_ct'][cell_type]
        binned_exposure=binning_output['binned_exposure_per_ct'][cell_type]
        binned_cell_type_mat=binning_output['binned_cell_type_mat']
        
        ct_ind=np.where(binning_output['cell_type_names']==cell_type)[0][0]
    
    segs=binning_output['segs']
    L=len(segs)

    fig,ax=plt.subplots(figsize=figsize)

    if layer_list is None:
        layer_list=range(L)
        
    for seg in layer_list:
        pts_seg=np.where(binned_labels==seg)[0]
        if cell_type is not None:
            pts_seg=[p for p in pts_seg if binned_cell_type_mat[p,ct_ind] / binned_cell_type_mat[p,:].sum() > spot_threshold]
        
        if colors is not None:
            plt.scatter(unique_binned_depths[pts_seg],
                       np.log( (binned_count[gene,pts_seg]) / binned_exposure[pts_seg] ), color=colors[seg], s=pt_size)
        else:
            c = np.array([color_palette(i) for i in range(len(np.unique(belayer_labels)))])
            plt.scatter(unique_binned_depths[pts_seg], 
                        np.log( (binned_count[gene,pts_seg]) / binned_exposure[pts_seg] ), color=c[seg], s=pt_size)
        if linear_fit:
            if cell_type is None:
                slope_mat, intercept_mat, _, _ = pw_fit_dict['all_cell_types']
            else:
                slope_mat, intercept_mat, _, _ = pw_fit_dict[cell_type]
            
            slope=slope_mat[gene,seg]
            intercept=intercept_mat[gene,seg]
            # print(slope,intercept)

            plt.plot( unique_binned_depths[pts_seg], intercept + slope*unique_binned_depths[pts_seg], color='black', alpha=1, linewidth=2 )
    plt.title(gene_name, fontsize=20)
    plt.xticks(fontsize=ticksize)
    plt.yticks(fontsize=ticksize)
    plt.savefig(f'./figures/{gene_name}_pwlinear.png', bbox_inches='tight', dpi=300)
    plt.close()
    
    
