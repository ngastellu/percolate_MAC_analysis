#!/usr/bin/env python

from os import path
import pickle
import numpy as np

# @njit
def network_data(nsample, T, rmax, run_name,gated=False):
    '''Obtains the radii of the sites in a given structure's percolating cluster, at a given temperature'''

    sites_dir = f'sample-{nsample}/sites_data_0.00105_psi_pow2'
    pkl_dir = f'sample-{nsample}/{run_name}_pkls'
    pkl_prefix = f'out_percolate_{run_name}'

    try:
        pkl_file = path.join(pkl_dir ,f'{pkl_prefix}-{T}K.pkl')
        fo = open(pkl_file, 'rb')
    except FileNotFoundError as e:
        #some of the runs have the 'K' after the run missing
        pkl_file = path.join(pkl_dir , f'{pkl_prefix}-{T}.pkl')
        fo = open(pkl_file,'rb')

    pkl = pickle.load(fo)
    fo.close()

    clusters = pkl[0]
    
    # Only expect one percolating cluster, but if there's more than one, take the uniion
    if len(clusters) > 1:
        icluster = np.array(list(clusters[0].union(*clusters[1:])))
    else:
        icluster = np.array(list(clusters[0])) 
    
    print(icluster.dtype)
    
    radii = np.load(path.join(sites_dir, 'radii.npy'))
    energies = np.load(path.join(sites_dir, 'ee.npy'))

    if gated:
        energies -= np.min(energies)

    # Percolation code runs on pre-filtered sites; therefore we must also filter the sites
    # to ensure that the indices in the cluster correspond to the correct radii
    filter = radii < rmax
    radii = radii[filter]
    energies = energies[filter]
    
    radii = radii[icluster]
    energies = energies[icluster]
    d = pkl[1]
    A = pkl[2]
    print(f'[{nsample}] Adjmat is symmetric: ', np.all((A == A.T)))
    
    nb_neighbours = A[icluster].sum(axis=1)
    #print(radii)

    return radii, energies, d , nb_neighbours


# @njit
def gather_data(nn, T, rmax, run_name, gated=False):
    nMACs = len(nn)
    ntot = 250*nMACs # expect about 250 sites per cluster for each structure

    all_radii = np.zeros(ntot)
    all_energies = np.zeros(ntot)
    dcrits = np.zeros(nMACs)
    nb_neighbours = np.zeros(ntot,dtype='int')

    nsites = 0
    for k, n in enumerate(nn):
        radii, energies, d, degs = network_data(n,T,rmax,run_name, gated=gated)
        #print(radii)
        n_new = radii.shape[0]
        dcrits[k] = d

        if nsites+n_new < ntot:
            all_radii[nsites:nsites+n_new] = radii
            all_energies[nsites:nsites+n_new] = energies
            nb_neighbours[nsites:nsites+n_new] = degs
        
        else:
            ntot += 50*nMACs

            tmp = np.zeros(ntot)
            tmp[:nsites] = all_radii[:nsites]
            tmp[nsites:nsites+n_new] = radii
            all_radii = tmp
            
            tmp = np.zeros(ntot)
            tmp[:nsites] = all_energies[:nsites]
            tmp[nsites:nsites+n_new] = energies
            all_energies = tmp
            
            tmp = np.zeros(ntot,dtype='int')
            tmp[:nsites] = nb_neighbours[:nsites]
            tmp[nsites:nsites+n_new] = degs
            nb_neighbours = tmp
        
        nsites += n_new
    
    return all_radii[:nsites], all_energies[:nsites],dcrits, nb_neighbours[:nsites]

def gather_radii(nn, T, rmax, run_name, gated=False):
    return gather_data(nn, T, rmax, run_name, gated=gated)[0] 

def gather_energies(nn, T, rmax, run_name, gated=False):
    return gather_data(nn, T, rmax, run_name, gated=gated)[1] 



if __name__ == "__main__":
    import sys
    import os


    structype=sys.argv[1]
    motype=sys.argv[2]
    T = int(sys.argv[3])
    gated=True

    if structype == '40x40':
        rmax = 18.03
    elif structype == 'tempdot6':
        rmax=121.2
    elif structype == 'tempdot5':
        rmax = 198.69
    else:
        print(f'Structure type {structype} is invalid. Exiting angrily.')
        sys.exit()
    

    run_name = f'rmax_{rmax}_psipow2_sites_gammas'#_{motype}'
    
    with open(f'{run_name}_symlinks/good_runs_{run_name}.txt') as fo:
    #with open(f'to_local_{run_name}/good_runs_{run_name}.txt') as fo:
        nn = [int(l.strip()) for l in fo.readlines()]
    
    cluster_energies = gather_energies(nn, T, rmax, run_name,gated=gated)

    outdir = f'cluster_energies/'

    if not os.path.exists(outdir):
        os.mkdir(outdir)
    
    print(cluster_energies)
    print(cluster_energies.shape)
    np.save(os.path.join(outdir, f'clust_energies_{run_name}-{T}K.npy'), cluster_energies)

    
