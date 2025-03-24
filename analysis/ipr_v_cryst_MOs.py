#!/usr/bin/env python

from time import perf_counter
import numpy as np
import os
from qcnico.qchemMAC import inverse_participation_ratios
from MO_crystallinity import MOs_crystallinity

ensemble = os.path.basename(os.path.getcwd())
renormalise = False # renormalise crystallinity by fraction of crystalline atoms in each structure

if ensemble == '40x40':
    lbls = np.arange(1,300)
elif ensemble == 'tempdot6':
    lbls = np.arange(218)
elif ensemble == 'tempdot5':
    lbls = np.arange(217)
else:
    print(f'{ensemble} is an invalid ensemble name.')


for n in lbls:
    start = perf_counter()
    try:
        cryst_mask_dir = os.path.expanduser(f'~/scratch/structural_characteristics_MAC/labelled_ring_centers/{ensemble}/sample-{n}/')
        cryst_mask = np.load(cryst_mask_dir +  f'crystalline_atoms_mask-{n}.npy')
    except FileNotFoundError as e:
        print(e)
        continue
    for motype in ['lo', 'virtual_w_HOMO', 'hi']:
        try:
            if motype == 'virtual_w_HOMO':
                M = np.load(os.path.expanduser(f'~/scratch/ArpackMAC/{ensemble}/MOs/{motype}/MOs_ARPACK_bigMAC-{n}.npy'))
            else:
                M = np.load(os.path.expanduser(f'~/scratch/ArpackMAC/{ensemble}/MOs/{motype}/MOs_ARPACK_{motype}_{ensemble}-{n}.npy'))
            N = M.shape[0]
        except FileNotFoundError as e:
            print(e)
            continue

    

    electronic_crystallinity = MOs_crystallinity(M,cryst_mask,renormalise)
    iprs = inverse_participation_ratios(M)

    np.save(f'MO_ipr_v_MO_cryst/{motype}/{n}.npy', np.hstack((iprs, electronic_crystallinity)))
    end = perf_counter()
    print(f'{n} [{end-start} seconds]', flush=True)