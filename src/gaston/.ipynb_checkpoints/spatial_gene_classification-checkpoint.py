import numpy as np
from collections import defaultdict

def get_discont_genes(pw_fit_dict, binning_output, q=0.95):
    
    _,_,discont_mat,_=pw_fit_dict['all_cell_types']
    gene_labels_idx=binning_output['gene_labels_idx']
    
    discont_genes=defaultdict(list) # gene -> [domain boundary p], ie bounary between R_p and R_{p+1}
    
    discont_q=np.quantile(np.abs(discont_mat), q,0)
    K=len(discont_q)
    for i,g in enumerate(gene_labels_idx):
        for l in range(K):
            if np.abs(discont_mat[i,l]) > discont_q[l]:
                #if g not in discont_genes:
                #    discont_genes[g]=[l]
                #else:
                discont_genes[g].append(l)
    
    # discont_genes=list( np.where(np.sum(np.abs(discont_mat) > discont_q,1))[0] )    

    return discont_genes

def get_cont_genes(pw_fit_dict, binning_output, q=0.95, ct_attributable=False, domain_cts=None, ct_perc=0.6):
    
    cont_genes=defaultdict(list) # dict of gene -> [list of domains]
    gene_labels_idx=binning_output['gene_labels_idx']
    

    slope_mat_all,_,_,_=pw_fit_dict['all_cell_types']
    slope_q=np.quantile(np.abs(slope_mat_all), q,0)
    
    L=len(slope_q)
    for i,g in enumerate(gene_labels_idx):
        for l in range(L):
            if np.abs(slope_mat_all[i,l]) > slope_q[l]:
                #if g not in cont_genes:
                #    cont_genes[g]=[l]
                #else:
                cont_genes[g].append(l)
    
    if not ct_attributable:
        return cont_genes
    
    cont_genes_domain_ct={g: [] for g in cont_genes} # dict gene -> [(domain,ct)]

    for g in cont_genes:
        for l in cont_genes[g]:
            other=True
            for ct in domain_cts[l]:
                if np.abs( pw_fit_dict[ct][0][gene_labels_idx==g,l] ) / np.abs(pw_fit_dict['all_cell_types'][0][gene_labels_idx==g,l]) > ct_perc:
                    other=False
                    cont_genes_domain_ct[g].append( (l,ct) )
                
            if other:
                cont_genes_domain_ct[g].append( (l, 'Other') )
                
    return cont_genes_domain_ct


# EXAMPLE of domain_cts: {0: ['Oligodendrocytes'], 1: ['Granule'], 2: ['Purkinje', 'Bergmann'], 3: ['MLI1', 'MLI2']}
# if other=True, then get cont genes NOT attributed to cell type
def get_cell_type_cont_genes(pw_fit_dict, domain_cts, binning_output, perc_threshold=0.6, q=0.95, other=False):
    
    # first get SVGs, as well as the domain they vary in
    cont_genes=set([]) # list of genes
    cont_genes_domain=[] # list of (genes, domain)
    gene_labels_idx=binning_output['gene_labels_idx']
    

    slope_mat_all,_,_,_=pw_fit_dict['all_cell_types']
    slope_q=np.quantile(np.abs(slope_mat_all), q,0)
    
    L=len(slope_q)
    for i,g in enumerate(gene_labels_idx):
        for l in range(L):
            if np.abs(slope_mat_all[i,l]) > slope_q[l]:
                cont_genes.add(g)
                cont_genes_domain.append((g,l))
    
    cont_genes_ct=set([])
    cont_genes_domain_ct=[] # dict (gene, domain): ct

    for g,l in cont_genes_domain:
        for ct in domain_cts[l]:
            if np.abs( pw_fit_dict[ct][0][gene_labels_idx==g,l] ) / np.abs(pw_fit_dict['all_cell_types'][0][gene_labels_idx==g,l]) > perc_threshold:
                cont_genes_ct.add(g)
                cont_genes_ct_domain[(g,l)] = ct
                
    cont_genes_other=[g for g in cont_genes if g not in cont_genes_ct]
    cont_genes_other_domain=[(g,l) for g,l in cont_genes_domain if (g,l) not in cont_genes_ct_domain]
    
    return cont_genes_ct,cont_genes_ct_domain, cont_genes_other, cont_genes_other_domain