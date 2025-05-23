#!/usr/bin/env python

from itertools import combinations, starmap
from functools import partial
import numpy as np
from numba import njit
import qcnico.qchemMAC as qcm
from qcnico.graph_tools import components
from scipy import sparse

@njit
def energy_distance(ei, ej, mu, T):
    kB = 8.617e-5
    return np.exp( (np.abs(ei-mu) + np.abs(ej-mu) + np.abs(ej-ei)) / (2*kB*T) )

@njit
def miller_abrahams_distance(ei, ej, ri, rj, mu, T, a0):
    kB = 8.617e-5 #eV/K
    return np.exp( ((np.abs(ei-mu) + np.abs(ej-mu) + np.abs(ej-ei)) / (2*kB*T)) + 2*np.linalg.norm(ri - rj)/a0 )

@njit
def log_miller_abrahams_distance(ei, ej, ri, rj, mu, T, a0):
    kB = 8.617e-5 #eV/K
    return ((np.abs(ei-mu) + np.abs(ej-mu) + np.abs(ej-ei)) / (2*kB*T)) + 2*np.linalg.norm(ri - rj)/a0

@njit
def diff_arrs(e, coms, a0, eF=0, E=np.array([0.0,0.0])):
    N = e.shape[0]
    ddarr = np.zeros(int(N*(N-1)/2))
    edarr = np.zeros(int(N*(N-1)/2))
    k = 0
    for i in range(N):
        ei = e[i] - coms[i].dot(E)
        for j in range(i):
            ej = e[j] - coms[j].dot(E)
            edarr[k] = (np.abs(ei-eF) + np.abs(ej-eF) + np.abs(ei - ej)) * 0.5
            ddarr[k] = 2*np.linalg.norm(coms[i]-coms[j])/a0
            k += 1
    return edarr, ddarr

@njit
def diff_arrs_w_inds(e, coms, a0, eF=0, E=np.array([0.0,0.0])):
    N = e.shape[0]
    ddarr = np.zeros(int(N*(N-1)/2))
    edarr = np.zeros(int(N*(N-1)/2))
    inds = np.zeros((int(N*(N-1)/2),2), dtype='int')
    k = 0
    for i in range(N):
        ei = e[i] - coms[i].dot(E)
        for j in range(i):
            ej = e[j] - coms[j].dot(E)
            edarr[k] = (np.abs(ei-eF) + np.abs(ej-eF) + np.abs(ei - ej)) * 0.5
            ddarr[k] = 2*np.linalg.norm(coms[i]-coms[j])/a0
            inds[k,0] = i
            inds[k,1] = j
            k += 1
    return edarr, ddarr, inds

@njit 
def diff_arrs_var_a(e,coms,radii,eF=0,E=np.array([0.0,0.0]),include_r_prefactor=True):
    N = e.shape[0]
    ddarr = np.zeros(int(N*(N-1)/2))
    edarr = np.zeros(int(N*(N-1)/2))
    inds = np.zeros((int(N*(N-1)/2),2), dtype='int')
    k = 0
    for i in range(N):
        ei = e[i] - coms[i].dot(E)
        for j in range(i):
            ej = e[j] - coms[j].dot(E)
            edarr[k] = (np.abs(ei-eF) + np.abs(ej-eF) + np.abs(ei - ej)) * 0.5
            sig_i = radii[i]
            sig_j = radii[j]
            radial_ij = (np.sum((coms[i,:] - coms[j,:])**2)) / (sig_i**2 + sig_j**2)
            if include_r_prefactor:
                radial_ij -= 2*np.log(sig_i*sig_j/(2*np.pi*(sig_i**2 + sig_j**2)))
            ddarr[k] = radial_ij
            inds[k,0] = i
            inds[k,1] = j
            k += 1
    return edarr, ddarr, inds



@njit
def dArray_energy(e, T, eF=None):
    N = e.size
    if eF is None:
        eF = 0.5 * (e[N//2 - 1] + e[N//2])
    darr = np.zeros(N*(N-1)//2)
    k = 0

    for i in range(1,N):
        for j in range(i):
            darr[k] = energy_distance(e[i], e[j], eF, T)
            k += 1
    
    return darr

@njit
def dArray_MA(e, coms, T, a0=1, eF=None):
    N = e.size
    if eF is None:
        eF = 0.5 * (e[N//2 - 1] + e[N//2])
    darr = np.zeros(N*(N-1)//2)
    k = 0

    for i in range(1,N):
        for j in range(i):
            darr[k] = miller_abrahams_distance(e[i], e[j], coms[i], coms[j], eF, T, a0)
            k += 1

    return darr

@njit
def dArray_logMA(e, coms, T, a0=1, eF=None):
    N = e.size
    if eF is None:
        eF = 0.5 * (e[N//2 - 1] + e[N//2])
    darr = np.zeros(N*(N-1)//2)
    k = 0

    for i in range(1,N):
        for j in range(i):
            darr[k] = log_miller_abrahams_distance(e[i], e[j], coms[i], coms[j], eF, T, a0)
            k += 1

    return darr

def distance_array_itertools(e, T):
    N = e.size
    print(N)
    eF = 0.5 * (e[N//2 - 1] + e[N//2])
    
    e_pairs = combinations(e,2)
    return np.array(list(starmap(partial(energy_distance, mu=eF, T=T), e_pairs)))

@njit
def pair_inds(n,N):
    zero_k_inds = np.array([k*(k-1)//2 for k in range(1,N)])
    i_inds = np.array([np.sum(nn >= zero_k_inds) for nn in n])
    return i_inds, (n - zero_k_inds[i_inds-1])

def k_ind(i,j): return int(i*(i-1)/2 + j)

# @njit
def percolate(darr, dpair_inds, L, R, return_adjmat=False, prev_d_ind=0): 
    """
    This function actually runs the percolation theory calculation.
    
    !!! *** VERY IMPORTANT *** !!!
    Only properties relative to the hopping sites (i.e. not the MOs,if we're gridifying the MOs to get sites out
    of them) should enter this function.

    Parameters
    ----------
    darr: `np.ndarray`, shape = (N*(N-1)/2,), dtype = `float`
        Flattened array of of pairwise 'distances' (usually log Miller-Abrahams rates) between sites
    dpair_inds: `np.ndarray`, shape = (N*(N-1)/2,2), dtype = `int`
        Indices of sites whose distances are stored in `darr`: `darr[k]` is the distance between sites `i` and `j`
        where  `(i,j) = dpair_inds[k,:]`
    L: `set` of `ints`
        Set of indices of sites strongly coupled to the left lead
    R: `set` of `int`s
        Set of indices of sites strongly coupled to the left lead
    return_adjmat: `bool` default=False
        If `True`, this function will return the adjacency matrix of the hopping sites once percolation is achieved 
        (useful for plotting the final clusters after the run).
    
    Returns
    -------
        spanning_clusters: `list` of `set`s of `int`s
            List of sets of indices of sites which form the spanning clusters. 1 set <--> 1 cluster
        d: `float`
            Percolation threshold distance
        adjmat: `np.ndarray`, shape = (N,N), dtype=`bool`
            Adjacency matrix of the sites at the percolation threshold. Returned only if `return_adjmat = True`.
    """
    
    percolated = False
    darr_sorted = np.sort(darr)
    ndists = darr.shape[0]
    n = int( (1+np.sqrt(1 + 8*ndists))/2 )
    adj_mat = np.zeros((n,n),dtype=bool)
    spanning_clusters = []
    d_ind = prev_d_ind
    while (not percolated) and (d_ind < ndists):
        d = darr_sorted[d_ind] #start with smallest distance and move up                                                                                                                                              
        print('d = ', d)       
        connected_inds = (darr < d).nonzero()[0] #darr is 1D array     
        print('Nb. of connected pairs = ', len(connected_inds))
        # ij = pair_inds(connected_inds,n)
        # print(ij)
        ij = dpair_inds[connected_inds,:]
        adj_mat[ij[:,0],ij[:,1]] = True
        adj_mat  += adj_mat.T
        print(adj_mat.astype(int))
        print('\n')
        relevant_MOs = set(np.unique(ij))
        coupledL = not L.isdisjoint(relevant_MOs)
        coupledR = not R.isdisjoint(relevant_MOs)
        if coupledL and coupledR:
            print('Getting clusters...')
            # LR_connected_inds = (L|R) & relevant_MOs 
            clusters = components(adj_mat) #seed only at edges of connected cluster
            print('Done! Now looping over clusters...')
            print(f'Nb. of clusters with more than 1 MO = {np.sum(np.array([len(c) for c in clusters])>1)}')
            for c in clusters:
                #print('Size of cluster: ', len(c))
                #print('Cluster: ', c)
                if (not c.isdisjoint(L)) and (not c.isdisjoint(R)):
                    spanning_clusters.append(c)
                    percolated = True
                    # if first_try:
                    #     print('First try!')
            print('Done with clusters loop!')
        
        d_ind += 1
        # d = darr_sorted[d_ind]
        # first_try = False
    
    if d_ind == ndists-1:
        print(f'No spanning clusters found!')
        return clusters, d, adj_mat, d_ind
    
    if return_adjmat:
        return spanning_clusters, d, adj_mat, d_ind -1
    else:
        return spanning_clusters, d, d_ind - 1


@njit
def avg_nb_neighbours(energies, dists, erange, urange):
    ebins = np.searchsorted(erange, energies) # figure out where the structure's energies fall in erange
    print('ebins     = ', ebins)
    n_e = erange.shape[0] + 1 # add 1 in case elements of energies are greater than max(erange)
    n_u = urange.shape[0]
    B = np.zeros((n_e,n_u))

    for k in range(n_u):
        print(k)
        u = urange[k]
        n_accessible = (dists <= u).sum(axis=1) 
        print(n_accessible)
        
        for n in range(n_e):
            imatch = (ebins == n).nonzero()[0]
            noccurences = imatch.shape[0]
            
            if noccurences == 0:
                # if no energies from the strcuture at hand fall into the nth energy bin, keep going
                continue
            # B[n,k] = np.array([n_accessible[i] for i in imatch]).sum() / noccurences
            B[n,k] = n_accessible[imatch].sum() / noccurences
    return B


# ------------------------ Jitted copies of percolation code here (above is obsolete) ------------------------

@njit
def jitted_components(M, seed_nodes=None):
    '''The same function as `components` in `qcnico.graph_tools`, but modified to work with Numba: now only works with dense matrix arrays.'''

    N = M.shape[0]
    seen = set()
    clusters = []
    if seed_nodes is None:
        seed_nodes = set(range(N)) #if no seed nodes are given, initiate cluster search over the entire graph
    for i in seed_nodes:
        if i not in seen:
            # do a breadth-first search starting from each unseen node
            c = set()
            nextlevel = {i}
            while nextlevel:
                thislevel = nextlevel
                nextlevel = set()
                for j in thislevel:
                    if j not in c:
                        c.add(j)
                        nextlevel.update(M[j,:].nonzero()[0])
            seen.update(c)
            clusters.append(c)
    return clusters


@njit
def jitted_percolate(darr, dpair_inds, L, R): 
    """
    This function actually runs the percolation theory calculation.
    It is a modified version of the above function `percolate`, massaged into working with Numba. Should
    also be more memory-efficient.
    
    !!! *** VERY IMPORTANT *** !!!
    Only properties relative to the hopping sites (i.e. not the MOs,if we're gridifying the MOs to get sites out
    of them) should enter this function.

    ~~~ Possible improvement ~~~
    Instead of scanning over each distance in ascending order, perform a binary search.

    Parameters
    ----------
    darr: `np.ndarray`, shape = (N*(N-1)/2,), dtype = `float`
        Flattened array of of pairwise 'distances' (usually log Miller-Abrahams rates) between sites
    L: `set` of `ints`
        Set of indices of sites strongly coupled to the left lead
    R: `set` of `int`s
        Set of indices of sites strongly coupled to the left lead
    return_adjmat: `bool` default=False
        If `True`, this function will return the adjacency matrix of the hopping sites once percolation is achieved 
        (useful for plotting the final clusters after the run).
    
    Returns
    -------
        spanning_clusters: `list` of `set`s of `int`s
            List of sets of indices of sites which form the spanning clusters. 1 set <--> 1 cluster
        d: `float`
            Percolation threshold distance
        adjmat: `np.ndarray`, shape = (N,N), dtype=`bool`
            Adjacency matrix of the sites at the percolation threshold. Returned only if `return_adjmat = True`.
    """

    percolated = False
    darr_sorted = np.sort(darr)
    ndists = darr.shape[0]
    n = int( (1+np.sqrt(1 + 8*ndists))/2 )
    adj_mat = np.zeros((n,n),dtype='bool')
    spanning_clusters = []
    d_ind = 0
    while (not percolated) and (d_ind < ndists):
        d = darr_sorted[d_ind] #start with smallest distance and move up                                                                                                                                              
        connected_inds = (darr < d).nonzero()[0] #darr is 1D array     
        ij = dpair_inds[connected_inds,:]
        for n in range(ij.shape[0]):
            i,j = ij[n,:]
            adj_mat[i,j] = True
            adj_mat[j,i] = True
        relevant_MOs = set(np.unique(ij))
        coupledL = not L.isdisjoint(relevant_MOs)
        coupledR = not R.isdisjoint(relevant_MOs)
        if coupledL and coupledR:
            LR_connected_inds = (L|R) & relevant_MOs 
            clusters = jitted_components(adj_mat, seed_nodes=LR_connected_inds) #seed only at edges of connected cluster
            for c in clusters:
                if (not c.isdisjoint(L)) and (not c.isdisjoint(R)):
                    spanning_clusters.append(c)
                    percolated = True
        
        d_ind += 1
    
    if d_ind == ndists-1:
        print(f'No spanning clusters found!')
        return clusters, d, adj_mat
    
    return spanning_clusters, d, adj_mat

if __name__ == '__main__':
    from os import path
    import sys
    import pickle
    from qcnico.qcffpi_io import read_energies, read_MO_file

    qcffpidir = path.expanduser('~/scratch/MO_airebo_dynamics/qcffpi')
    mo_file = 'MO_coefs.dat'
    energy_file = 'orb_energy.dat'
    
    n1 = int(sys.argv[1])
    n2 = int(sys.argv[2])
    step = int(sys.argv[3])

    frames = np.arange(n1,n2,step)
    nframes = frames.size

    cluster_list = [None] * nframes
    ds = np.zeros((nframes,2))
    for k, i in range(frames):
        framedir= f'300K_initplanar_norotate/frame-{i}'
        mo_path = path.join(qcffpidir,framedir,mo_file)
        energy_path = path.join(qcffpidir,framedir,energy_file)

        energies = read_energies(energy_path)
        pos, M = read_MO_file(mo_path)

        clusters, d = percolate(energies,pos,M)

        cluster_list[k] = clusters
        ds[k] = i, d
    
    np.save(ds, f'ds_frames_{n1}-{n2}-{step}.npy')
    with open('clusters_frames_{n1}-{n2}-{step}.pkl', 'wb') as fo:
        pickle.dump(cluster_list, fo)

