import numpy as np
import scanpy as sc
import squidpy as sq
import pandas as pd

def get_gaston_input_adata(data_folder, get_rgb=False, spot_umi_threshold=50):
    adata = sc.read_10x_h5(f'{data_folder}/filtered_feature_bc_matrix.h5')

    df_pos = pd.read_csv(f'{data_folder}/spatial/tissue_positions_list.csv', sep=",", header=None, names=["barcode", "in_tissue", "array_row", "array_col", "pxl_row_in_fullres", "pxl_col_in_fullres"])
    df_pos = df_pos[df_pos.in_tissue == True]
    assert set(list(df_pos.barcode)) == set(list(adata.obs.index))
    df_pos.barcode = pd.Categorical(df_pos.barcode, categories=list(adata.obs.index), ordered=True)
    df_pos.sort_values(by="barcode", inplace=True)
    adata.obsm["coords"] = np.vstack([df_pos.array_row, df_pos.array_col]).T
    
    gene_labels = adata.var.index.to_numpy()
    counts_mat = adata.X
    counts_mat = np.array(counts_mat.todense())
    coords_mat = adata.obsm["coords"]
    
    # filter spots with low UMI count, based on prior knowledge
    print(f'removing spots with UMI count < {spot_umi_threshold}')
    nonzero_spots = np.where(np.sum(counts_mat,1) >= spot_umi_threshold)[0]
    counts_mat = counts_mat[nonzero_spots,:]
    coords_mat = coords_mat[nonzero_spots,:]
    df_pos = df_pos.iloc[nonzero_spots,:]

    if not get_rgb:
        return counts_mat, coords_mat, gene_labels

    # get mean RGB values per spot
    adata2 = sc.read_visium(data_folder,
                        count_file = f"filtered_feature_bc_matrix.h5",
                        source_image_path = f"spatial/tissue_hires_image.png")

    library_id = list(adata2.uns['spatial'].keys())[0] # adata2.uns['spatial'] should have only one key
    scale=adata2.uns['spatial'][library_id]['scalefactors']['tissue_hires_scalef']
    img = sq.im.ImageContainer(adata2.uns['spatial'][library_id]['images']['hires'],
                               scale=scale,
                               layer="img1")
    print('calculating RGB')
    sq.im.calculate_image_features(adata2, img, features="summary", key_added="features")
    columns = ['summary_ch-0_mean', 'summary_ch-1_mean', 'summary_ch-2_mean']
    RGB_mean = adata2.obsm["features"][columns]
    
    RGB_mean = RGB_mean[RGB_mean.index.isin(df_pos['barcode'])]

    return counts_mat, coords_mat, gene_labels, RGB_mean.to_numpy()

def get_top_pearson_residuals(num_pcs, counts_mat, coords_mat, gene_labels=None, n_top_genes=5000, clip=0.01):
    df=pd.DataFrame(counts_mat, columns=gene_labels)
    adata=sc.AnnData(df)
    adata.obsm["coords"] = coords_mat
    
    sc.experimental.pp.highly_variable_genes(
            adata, flavor="pearson_residuals", n_top_genes=n_top_genes
    )
    
    adata = adata[:, adata.var["highly_variable"]]
    adata.layers["raw"] = adata.X.copy()
    adata.layers["sqrt_norm"] = np.sqrt(
        sc.pp.normalize_total(adata, inplace=False)["X"]
    )
    
    theta=np.Inf
    sc.experimental.pp.normalize_pearson_residuals(adata, clip=clip, theta=theta)
    sc.pp.pca(adata, n_comps=num_pcs)
    return adata.obsm['X_pca']

def get_gaston_input_xenium(folder,filter_zero_cells=False):
    adata = sc.read_10x_h5(filename=folder+'cell_feature_matrix.h5')

    df = pd.read_csv(folder+"cells.csv")
    df.set_index(adata.obs_names, inplace=True)
    adata.obs = df.copy()

    adata.obsm["spatial"] = adata.obs[["x_centroid", "y_centroid"]].copy().to_numpy()

    counts_mat=np.array(adata.X.todense())
    coords_mat=adata.obs[["x_centroid", "y_centroid"]].to_numpy()

    cts=pd.read_csv(folder+'Cell_Barcode_Type_Matrices.csv')
    cts=np.array(cts['Cluster'])

    gene_labels=np.array(adata.var.index)

    if filter_zero_cells:
        inds_to_keep=np.where(np.sum(counts_mat,1)>0)[0]
        counts_mat=counts_mat[inds_to_keep,:]
        coords_mat=coords_mat[inds_to_keep,:]
        cts=cts[inds_to_keep]

    return counts_mat, coords_mat, gene_labels, cts