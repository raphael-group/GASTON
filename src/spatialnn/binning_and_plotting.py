from scipy.stats import mode
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import os

from spatialnn import segmented_fit

# counts mat has to be of shape G x N
def bin_data(counts_mat, belayer_labels, belayer_depth, 
              cell_type_df, gene_labels, num_bins=70, num_bins_per_layer=None,
             idx_kept=None, umi_threshold=500, pc=1, extra_data=[]):
    
    if idx_kept is None:
        idx_kept=np.where(np.sum(counts_mat,1) > umi_threshold)[0]
    gene_labels_idx=gene_labels[idx_kept]
    
    # pseudo_counts_mat=counts_mat+1
    # exposure=np.sum(pseudo_counts_mat,axis=0)
    pseudo_counts_mat, exposure = segmented_fit.add_pc(counts_mat, pc=pc)
    
    cmat=pseudo_counts_mat[idx_kept,:]
    
    if cell_type_df is not None:
        cell_type_mat=cell_type_df.to_numpy()
        cell_type_names=np.array(cell_type_df.columns)
    else:
        N=len(exposure)
        cell_type_mat=np.ones((N,1))
        cell_type_names=['All']

    G,N=cmat.shape

    # BINNING
    if num_bins_per_layer is not None:
        bins=np.array([])
        L=len(np.unique(belayer_labels))
        
        for l in range(L):
            depth_l=belayer_depth[np.where(belayer_labels==l)[0]]
            
            if l>0:
                depth_lm1=belayer_depth[np.where(belayer_labels==l-1)[0]]
                depth_left=0.5*(np.min(depth_l) + np.max(depth_lm1))
            else:
                depth_left=np.min(depth_l)-0.01
                
            if l<L-1:
                depth_lp1=belayer_depth[np.where(belayer_labels==l+1)[0]]
                depth_right=0.5*(np.max(depth_l) + np.min(depth_lp1))
            else:
                depth_right=np.max(depth_l)+0.01
            
            bins_l=np.linspace(depth_left, depth_right, num=num_bins_per_layer[l]+1)
            if l!=0:
                bins_l=bins_l[1:]
            bins=np.concatenate((bins, bins_l))
    else:
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
    binned_number_spots=np.zeros(N_1d)

    binned_count_per_ct={ct: np.zeros( (G, N_1d) ) for ct in cell_type_names}
    binned_exposure_per_ct={ct: np.zeros( N_1d ) for ct in cell_type_names}
    binned_extra_data=[np.zeros(N_1d) for i in range(len(extra_data))]
    map_1d_bins_to_2d={} # map b -> [list of cells in bin b]
    for ind, b in enumerate(unique_binned_depths):
        bin_pts=np.where(binned_depths==b)[0]

        binned_count[:,ind]=np.sum(cmat[:,bin_pts],axis=1)
        binned_exposure[ind]=np.sum(exposure[bin_pts])
        binned_labels[ind]= int(mode( belayer_labels[bin_pts],keepdims=False).mode)
        binned_cell_type_mat[ind,:] = np.sum( cell_type_mat[bin_pts,:], axis=0)
        binned_number_spots[ind]=len(bin_pts)
        map_1d_bins_to_2d[b]=bin_pts

        for i, eb in enumerate(extra_data):
            binned_extra_data[i][ind]=np.mean(extra_data[i][bin_pts])
                
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
    to_return['binned_number_spots']=binned_number_spots
    
    to_return['binned_count_per_ct']=binned_count_per_ct
    to_return['binned_exposure_per_ct']=binned_exposure_per_ct
    to_return['binned_extra_data']=binned_extra_data
    
    to_return['map_1d_bins_to_2d']=map_1d_bins_to_2d
    to_return['segs']=segs

    return to_return

"""
def plot_gene_pwlinear(gene_name, pw_fit_dict, belayer_labels, belayer_depth, binning_output, 
                       cell_type=None, spot_threshold=0.25, pt_size=10, 
                       colors=None, linear_fit=True, lw=2, layer_list=None, ticksize=20, figsize=(7,3),
                      offset=1, save=False, save_dir="./"):
    
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
        if colors is None:
            c=None
        else:
            c=colors[seg]
        
        xax=unique_binned_depths[pts_seg]
        yax=np.log( (binned_count[gene,pts_seg] / binned_exposure[pts_seg]) * offset )
        
        plt.scatter(xax, yax, color=c, s=pt_size)
        
        if linear_fit:
            if cell_type is None:
                slope_mat, intercept_mat, _, _ = pw_fit_dict['all_cell_types']
            else:
                slope_mat, intercept_mat, _, _ = pw_fit_dict[cell_type]

            slope=slope_mat[gene,seg]
            intercept=intercept_mat[gene,seg]
            plt.plot( unique_binned_depths[pts_seg], np.log(offset) + intercept + slope*unique_binned_depths[pts_seg], color='grey', alpha=1, lw=lw )
                
    plt.xticks(fontsize=ticksize)
    plt.yticks(fontsize=ticksize)
    plt.title(f"{gene_name}")
    if save:
        plt.savefig(f"{save_dir}/{gene_name}_pwlinear.png", bbox_inches="tight")
        plt.close()
"""
    
    

def plot_gene_pwlinear(gene_name, pw_fit_dict, belayer_labels, belayer_depth, binning_output, 
                       cell_type_list=None, ct_colors=None, spot_threshold=0.25, pt_size=10, 
                       colors=None, linear_fit=True, lw=2, layer_list=None, ticksize=20, figsize=(7,3),
                      offset=1, xticks=None, yticks=None, alpha=1, layer_boundary_plotting=False, 
                      save=False, save_dir="./"):
    
    gene_labels_idx=binning_output['gene_labels_idx']
    if gene_name in gene_labels_idx:
        gene=np.where(gene_labels_idx==gene_name)[0]
    else:
        umi_threshold=binning_output['umi_threshold']
        raise ValueError(f'gene does not have UMI count above threshold {umi_threshold}')
    
    unique_binned_depths=binning_output['unique_binned_depths']
    binned_labels=binning_output['binned_labels']
    
    binned_count_list=[]
    binned_exposure_list=[]
    ct_ind_list=[]
    
    if cell_type_list is None:
        binned_count_list.append(binning_output['binned_count'])
        binned_exposure_list.append(binning_output['binned_exposure'])
        
    else:
        for ct in cell_type_list:
            binned_count_list.append(binning_output['binned_count_per_ct'][ct])
            binned_exposure_list.append(binning_output['binned_exposure_per_ct'][ct])
            ct_ind_list.append( np.where(binning_output['cell_type_names']==ct)[0][0] )
    
    segs=binning_output['segs']
    L=len(segs)

    fig,ax=plt.subplots(figsize=figsize)

    if layer_list is None:
        layer_list=range(L)
        
    for seg in layer_list:
        pts_seg=np.where(binned_labels==seg)[0]
        
        for i in range(len(binned_count_list)):
            binned_count=binned_count_list[i]
            binned_exposure=binned_exposure_list[i]
            ct=None
            if cell_type_list is not None:
                ct=cell_type_list[i]
                # if restricting cell types, then restrict spots also
                binned_cell_type_mat=binning_output['binned_cell_type_mat']
                ct_ind=ct_ind_list[i]
                pts_seg=[p for p in pts_seg if binned_cell_type_mat[p,ct_ind] / binned_cell_type_mat[p,:].sum() > spot_threshold]
                
                # set colors for cell types
                if ct_colors is None:
                    c=None
                else:
                    c=ct_colors[ct]
            else:
                # set colors for layers
                if colors is None:
                    c=None
                else:
                    c=colors[seg]
                

            xax=unique_binned_depths[pts_seg]
            yax=np.log( (binned_count[gene,pts_seg] / binned_exposure[pts_seg]) * offset )
            plt.scatter(xax, yax, color=c, s=pt_size, alpha=alpha)

            if linear_fit:
                if ct is None:
                    slope_mat, intercept_mat, _, _ = pw_fit_dict['all_cell_types']
                else:
                    slope_mat, intercept_mat, _, _ = pw_fit_dict[ct]

                slope=slope_mat[gene,seg]
                intercept=intercept_mat[gene,seg]
                plt.plot( unique_binned_depths[pts_seg], np.log(offset) + intercept + slope*unique_binned_depths[pts_seg], color='grey', alpha=1, lw=lw )
    
    if xticks is None:
        plt.xticks(fontsize=ticksize)
    else:
        plt.xticks(xticks,fontsize=ticksize)
        
    if yticks is None:
        plt.yticks(fontsize=ticksize)
    else:
        plt.yticks(yticks,fontsize=ticksize)
        
    if layer_boundary_plotting and len(layer_list)>1:
        binned_labels=binning_output['binned_labels']
        
        left_bps=[]
        right_bps=[]

        for i in range(len(binned_labels)-1):
            if binned_labels[i] != binned_labels[i+1]:
                left_bps.append(unique_binned_depths[i])
                right_bps.append(unique_binned_depths[i+1])
        
        for i in layer_list[:-1]:
            plt.axvline((left_bps[i]+right_bps[i])*0.5, color='black', ls='--', linewidth=1.5, alpha=0.2)

    sns.despine()
    if save:
        os.makedirs(save_dir, exist_ok=True)
        plt.savefig(f"{save_dir}/{gene_name}_pwlinear.pdf", bbox_inches="tight")
        plt.close()
