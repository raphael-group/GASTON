import numpy as np

import torch
import torch.nn as nn
import torch.nn.functional as F
import torch.utils
import torch.distributions

import matplotlib.pyplot as plt


def rotate_by_theta(coords, theta, rotate_about=np.array([0,0])):
    """
    Rotate coordinate array by angle theta.
    :param coords: np.array of shape (N,2) points (rows).
    :param theta: angle theta (in radians) to rotate by
    :param rotate_about: (OPTIONAL) point to rotate about
    :return: np.array of shape (N,2) where each point is rotated by theta
    """
    coordsT=coords.T
    
    c,s=np.cos(theta), np.sin(theta)
    rotation_matrix=np.array(((c, -s), (s, c)))
    
    return (rotation_matrix @ coordsT).T
    

def opt_linear(y, xcoords):
    """
    Fit each row y_g of Y as y_g = xcoords * theta + const using linear regression
    :param y: np.array of shape (G,N)
    :param xcoords: np.array of shape (N,) of x-coordinates.
    :return: (1) theta, const for each gene g and (2) Total loss \sum_g ||y_g - yhat_g||^2 
    """
    g,n = y.shape
    if n==1:
        return np.array([]),0
    
    X=np.ones((n,2))
    X[:,0]=xcoords

    try:
        theta=np.linalg.inv(X.T @ X) @ X.T @ y.T
        error=np.linalg.norm(X @ theta - y.T)**2

        # linreg coefficients theta.T have size g x 2
        # each row is the linreg coefficients for gene g
        return theta.T,error
    except: # in case no linreg fit
        placeholder=np.ones((2,g))-2
        return placeholder, np.Inf
    
def dp_raw(data, Lmax, xcoords, opt_function=opt_linear):
    """
    Dynamic programming algorithm for segmented regression
    :param data: np.array of shape (G,N), either integer counts or GLM-PCs
    :param Lmax: maximum number of domains
    :param xcoords: np.array of shape (N,) of x-coordinates of each spot
    :param opt_function: function to fit each segment (default is least squares for linear regression)
    :return: (1) error_mat: DP table of size N x L
             where error_mat[n,l] = error from using linreg to fit first n+1 
             spots (after sorting by x-coordinate) using (l+1)-piece segmented regression,
             (2) segment_map: pointers for DP table
             (n,l) -> (n', l-1) where you use l-1 pieces for spots 1, ..., n' 
             and one piece for spots n'+1, ..., N
    """

    G=data.shape[0]
    N=data.shape[1]
    
    sorted_xcoords=np.sort(xcoords)
    sorted_xcoords_inds=np.argsort(xcoords)

    # dp table error_mat (T x L)
    # where error_mat[t,l] = error from fitting first t+1 coords using (l+1)-piece segmented regression
    error_mat=np.zeros((N,Lmax))

    # map (t,p) -> (t', p-1) where you use p-1 segments for times 1,...,t' and a segment for t'+1,...,t
    segment_map={}

    # save previous calls to opt_function
    saved_opt_functions=np.zeros((N+1,N+1)) - 1

    # fill out first column of matrix (0th row is just 0)
    for n in range(N):
        xc=sorted_xcoords[:n+1]
        xc_inds=sorted_xcoords_inds[:n+1]
        
        _,err=opt_function(data[:,xc_inds],xc)
        error_mat[n,0]=err

    # fill out each subsequent column l
    for l in range(1,Lmax):

        # for each column, go from top to bottom [ignoring first row]
        for n in range(N):
            best_nprime=-1
            best_error=np.Inf

            for nprime in range(n):
                if saved_opt_functions[nprime+1,n+1] >= 0:
                    nprime_fit=saved_opt_functions[nprime+1,n+1]
                else:
                    xc=sorted_xcoords[nprime+1:n+1]
                    xc_inds=sorted_xcoords_inds[nprime+1:n+1]
                    
                    nprime_fit=opt_function(data[:,xc_inds], xc)[1]
                    saved_opt_functions[nprime+1,n+1]=nprime_fit
                cur_error=error_mat[nprime,l-1] + nprime_fit

                if cur_error < best_error:
                    best_error=cur_error
                    best_nprime=nprime
            error_mat[n,l] = best_error
            segment_map[(n,l)] = (best_nprime,l-1)

    return error_mat, segment_map


def dp_bucketized(data, bucket_endpoints, Lmax, xcoords, opt_function=opt_linear):
    """
    Dynamic programming algorithm for segmented regression
    where the x-coordinates are partitioned into B equally spaced buckets
    :param data: np.array of shape (G,N), either integer counts or GLM-PCs
    :param bucket_endpoints: np.array of shape (B+1), the endpoints of the B buckets
    :param Lmax: maximum number of domains
    :param xcoords: np.array of shape (N,) of x-coordinates for each spot
    :param opt_function: function to fit each segment (default is least squares for linear regression)
    :return: (1) error_mat: DP table of size B x L
             where error_mat[b,l] = error from using linreg to fit first 
             spots in first b+1 buckets using (l+1)-piece segmented regression,
             (2) segment_map: pointers for DP table
             (b,l) -> (b', l-1) where you use l-1 pieces for buckets 1,...,b' 
             and one piece for buckets b'+1,...,t
    """
    G=data.shape[0]
    N=data.shape[1]
    
    B=len(bucket_endpoints)-1
    buckets=np.digitize(xcoords, bucket_endpoints) - 1 # for some reason bucket labels are 1-indexed
    
    # dp on matrix error_mat of size B x Lmax
    
    # where error_mat[b,p] = error from using linreg to fit 
    # first b+1 buckets using (p+1)-part segmented regression
    error_mat=np.zeros((B,Lmax))

    # map (b,p) -> (b', p-1) where you use p-1 segments for buckets 1,...,b' and a segment for b'+1,...,b
    segment_map={}

    # saved opt_function
    saved_opt_functions=np.zeros((B+1,B+1)) - 1

    # fill out first column of matrix (0th row is just 0)
    for b in range(B):
        inds_upto_bucket_b = np.where(buckets <= b)[0]
        
        xc=xcoords[inds_upto_bucket_b]
        _,err=opt_function(data[:,inds_upto_bucket_b],xc)
        error_mat[b,0]=err

    # fill out each subsequent column p
    for p in range(1,Lmax):

        # for each column, go from top to bottom [ignoring first row]
        for b in range(B):
            # fill out entry t,p

            best_bprime=-1
            best_error=np.Inf
            for bprime in range(b):
                
                if saved_opt_functions[bprime+1,b+1] >= 0:
                    bprime_fit=saved_opt_functions[bprime+1,b+1]
                else:
                    inds_between_bprime_b=np.where( (buckets > bprime) & (buckets <= b) )[0]
                    
                    xc=xcoords[inds_between_bprime_b]                    
                    
                    bprime_fit=opt_function(data[:,inds_between_bprime_b], xc)[1]
                    saved_opt_functions[bprime+1,b+1]=bprime_fit
                    
                cur_error=error_mat[bprime,p-1] + bprime_fit
                if cur_error < best_error:
                    best_error=cur_error
                    best_bprime=bprime
            error_mat[b,p] = best_error
            segment_map[(b,p)] = (best_bprime,p-1)
    
    return error_mat, segment_map


def dp(data, isodepth, Lmax, use_buckets=True, num_buckets=150, opt_function=opt_linear):
    G, N = data.shape
    
    if use_buckets:
        bin_endpoints=np.linspace(np.min(isodepth),np.max(isodepth)+0.01,num_buckets+1)
        print('running DP with buckets')
        error_mat,seg_map=dp_bucketized(data, bin_endpoints, Lmax, 
                                        isodepth, opt_function=opt_linear)

        # get labels, save to res
        losses,labels=get_losses_labels(error_mat,seg_map, isodepth, 
                                            use_buckets=True, 
                                            bin_endpoints=bin_endpoints)
        return losses, labels

    else:
        print('running DP without buckets')
        error_mat, seg_map = dp_raw(data, Lmax, isodepth)

        # get labels, save to res
        losses,labels=get_losses_labels(error_mat,seg_map, isodepth, use_buckets=False)
        return losses, labels
    


def rotation_dp(data, coords, Lmax=8, rotation_angle_list=[0,5,10,15,17.5,20], 
                use_buckets=True, num_buckets=150, opt_function=opt_linear):
    """
    Rotation mode of gaston.
    Rotates the tissue slice by different angles and runs the above DP functions.
    :param data: np.array of shape (G,N), either integer counts or GLM-PCs
    :param coords: np.array of shape (N,2), i-th row is coordinate of spot i
    :param Lmax: maximum number of domains
    :param rotation_angle_list: list of angles (in degrees) to rotate tissue by
    :param use_buckets: boolean indicating whether to divide x-coordinates into 
                        equally spaced buckets or not. 
                        (Lower spatial resolution but faster run-time)
    :param num_buckets: (if use_buckets=True) number of buckets to divide 
                        x-coordinate range into
    :param opt_function: function to fit each segment (default is least squares for linear regression)
    :return: (1) loss_array: np.array of shape (number of rotation angles, Lmax)
            where losses[a,l] = DP loss from running segmented regression with l pieces
            on tissue rotated by angle rotation_angle_list[a]
            (2) label_dict: dict where [a,l] -> DP labels  from running segmented regression 
            with l pieces on tissue rotated by angle rotation_angle_list[a]
    """
    G = data.shape[0]
    N = data.shape[1]
    
    # loop over all angles in rotation_angle_list
    theta_list=[deg*np.pi/180 for deg in rotation_angle_list]
    
    # each row is [loss at domains 1, ..., Lmax] for each rotation angle
    losses=np.zeros( (len(rotation_angle_list), Lmax) )
    labels={}
    
    for ind_t,theta in enumerate(theta_list):
        print('\n angle: {}'.format(rotation_angle_list[ind_t]))
        coords_rotated=rotate_by_theta(coords,theta)
        xcoords_rotated=coords_rotated[:,0]
        if use_buckets:
            bin_endpoints=np.linspace(np.min(xcoords_rotated),np.max(xcoords_rotated)+0.01,num_buckets+1)
            print('running DP with buckets')
            error_mat,seg_map=dp_bucketized(data, bin_endpoints, Lmax, 
                                            xcoords_rotated, opt_function=opt_linear)

            # get labels, save to res
            losses_t,labels_t=get_losses_labels(error_mat,seg_map, xcoords_rotated, 
                                                use_buckets=True, 
                                                bin_endpoints=bin_endpoints)
            
            losses[ind_t,:]=losses_t
            for l in range(1,Lmax+1):
                labels[(ind_t,l)]=labels_t[l]
            
        else:
            print('running DP without buckets')
            error_mat, seg_map = dp_raw(data, Lmax, xcoords_rotated)
            
            # get labels, save to res
            losses_t,labels_t=get_losses_labels(error_mat,seg_map, xcoords_rotated, use_buckets=False)
            
            losses[ind_t,:]=losses_t
            for l in range(1,Lmax+1):
                labels[(ind_t,l)]=labels_t[l]
            
            
#             for l in range(1,Lmax+1):
#                 print('finding segments for {} domains'.format(l))
#                 segs=find_segments_from_dp(error_mat, seg_map, l, xcoords=xcoords_rotated)
#                 dp_labels=np.zeros(N)
#                 c=0
#                 for seg in segs:
#                     dp_labels[seg]=c
#                     c+=1
                
#                 losses[ind_t,l-1] = error_mat[-1,l-1] / N
#                 labels[(ind_t,l)]=dp_labels
    return losses,labels

def find_segments_from_dp(error_mat, segment_map, l, xcoords=None):
    """
    Backtrack through DP output to find the l segments
    :param error_mat, segment_map: outputs of dp_bucketized or dp_raw
    :param l: number of segments
    :param xcoords: (Optional) x-coordinates for each spot
    
    :return: array of l segments of DP
    """
    num_times=error_mat.shape[0]
    
    segs=[[] for i in range(l)]
    seg_val=l-1
    time_val=num_times-1
    
    if xcoords is None:
        xcoords=np.arange(num_times)
    
    sorted_xcoords=np.sort(xcoords)
    sorted_xcoords_inds=np.argsort(xcoords)

    while seg_val > 0:
        new_time_val,new_seg_val=segment_map[(time_val,seg_val)]
        segs[seg_val]=sorted_xcoords_inds[new_time_val+1:time_val+1]
        time_val=new_time_val
        seg_val=new_seg_val
    segs[0]=np.arange(0,time_val+1)
    return segs


# model: pytorch model
# A: numpy array of size N x F (rows are spots, cols are GLM-PCs)
# S: numpy array of size N x 2 (rows are spots, cols are coords)

# OUTPUT: labels
def get_isodepth_labels(model, A, S, num_domains, num_buckets=50, num_pcs_A=None):
    N=A.shape[0]
    S_torch=torch.Tensor(S)
    gaston_isodepth=model.spatial_embedding(torch.Tensor(S)).detach().numpy().flatten()
    
    kmax=num_domains
    
    if num_pcs_A is not None:
        A=A[:,:num_pcs_A]
    
    bin_endpoints=np.linspace(np.min(gaston_isodepth),np.max(gaston_isodepth)+0.01,num_buckets+1)
    error_mat,seg_map=dp_bucketized(A.T, bin_endpoints, kmax, xcoords=gaston_isodepth)
    bin_labels=np.digitize(gaston_isodepth,bin_endpoints)

    segs=find_segments_from_dp(error_mat, seg_map, num_domains)
    gaston_labels=np.zeros(N)
    c=0
    for seg in segs:
        for s in seg:
            gaston_labels[ np.where(bin_labels==s+1)[0] ] = c
        c+=1
    
    return gaston_isodepth,gaston_labels

def filter_rescale_boundary(counts_mat, coords_mat, gaston_isodepth, gaston_labels, isodepth_min=0, isodepth_max=1):
    locs=np.array( [i for i in range(len(gaston_isodepth)) if isodepth_min < gaston_isodepth[i] < isodepth_max] )
    
    counts_mat_subset=counts_mat[locs,:]
    print(f'restricting to {coords_mat.shape} spots')
    coords_mat_subset=coords_mat[locs,:]
    gaston_labels_subset=gaston_labels[locs]
    gaston_labels_subset -= np.min( np.unique( gaston_labels[locs] ) )
    gaston_isodepth_subset=gaston_isodepth[locs]
    
    ####
    
    gaston_isodepth_subset0=gaston_isodepth_subset[gaston_labels_subset==0]
    gaston_isodepth_subset1=gaston_isodepth_subset[gaston_labels_subset==1]

    isodepth_ratio=np.std( gaston_isodepth_subset1 ) / np.std( gaston_isodepth_subset0 )

    gaston_isodepth_subset0_new=((gaston_isodepth_subset0 - np.max(gaston_isodepth_subset0)) * isodepth_ratio) + np.max(gaston_isodepth_subset0)

    gaston_isodepth_subset[gaston_labels_subset==0] = gaston_isodepth_subset0_new
    
    return counts_mat_subset, coords_mat_subset, gaston_isodepth_subset, gaston_labels_subset
    
