#!/usr/bin/env python


import numpy as np
import os
import sys
from qcnico.coords_io import read_xyz
from qcnico.data_utils import save_npy
from qcnico.qchemMAC import MO_rgyr



# ---------- Initialise stuff ----------

structype = os.path.basename(os.getcwd())
runtype = 'sites_pure'

motype = sys.argv[1]


if structype == '40x40':
    rmax = 18.03
    lbls = np.arange(1,301)
elif structype == 'tempdot6':
    rmax=121.2
    lbls = np.arange(132)
elif structype == 'tempdot5':
    rmax = 198.69
    lbls = np.arange(117)
else:
    print(f'Structure type {structype} is invalid. Exiting angrily.')
    sys.exit()

for nn in lbls:
    print(f'Doing {nn}...', end=' ', flush=True)
    if runtype not in ['MOs_pure', 'sites_pure', 'mixed']:
        print(f'Invalid `runtype` {runtype}. Valid entries are:\n* "mixed": sites created by k-clustering, MO gammas;\n*"MOs_pure": MOs directly as sites;\n*"sites_pure": sites created by k-clustering, gammas obtained from site kets.\nExiting angrily.')
        sys.exit()

    gammas_dir = f'sample-{nn}/gammas/'
    gamma_npy_suffix = f'{runtype}_{motype}'

    # ---------- Get radii ----------

    if runtype == 'sites_pure' or runtype == 'mixed':

        gamma_npy_suffix = f'rmax_{rmax}_' + gamma_npy_suffix

        if motype == 'virtual':
            sitesdir = f'sample-{nn}/sites_data_0.00105_psi_pow2/'
        else:
            sitesdir = f'sample-{nn}/sites_data_0.00105_psi_pow2_{motype}/'

        try:
            radii = np.load(os.path.join(sitesdir, f'radii.npy'))   
            radii = radii[radii < rmax]
        except FileNotFoundError as e:
            print('Radii NPY missing!')
            continue

    else:
        try:
            if motype == 'virtual':
                M = np.load(os.path.expanduser(f"~/scratch/ArpackMAC/{structype}/MOs/{motype}/MOs_ARPACK_bigMAC-{nn}.npy"))
            else:
                M = np.load(os.path.expanduser(f"~/scratch/ArpackMAC/{structype}/MOs/{motype}/MOs_ARPACK_{motype}_{structype}-{nn}.npy"))
        except FileNotFoundError as e:
            print('MO file missing!')
            continue

        strucdir = os.path.expanduser(f"~/scratch/clean_bigMAC/{structype}/relaxed_structures_no_dangle/")

        if structype == '40x40':
            rmax = 18.03
            strucfile = os.path.join(strucdir, f'bigMAC-{nn}_relaxed_no-dangle.xyz')
        elif structype == 'tempdot6':
            rmax=121.2
            strucfile = os.path.join(strucdir, f'{structype}n{nn}_relaxed_no-dangle.xyz')
        elif structype == 'tempdot5':
            rmax = 198.69
            strucfile = os.path.join(strucdir, f'{structype}n{nn}_relaxed_no-dangle.xyz')
        else:
            print(f'Structure type {structype} is invalid. Exiting angrily.')
            sys.exit()

        pos = read_xyz(strucfile)

        radii = MO_rgyr(pos, M) 

    # ---------- Load gammas ----------
    try:
        gamL = np.load(os.path.join(gammas_dir,f'gamL_{gamma_npy_suffix}.npy'))
        gamR = np.load(os.path.join(gammas_dir,f'gamR_{gamma_npy_suffix}.npy'))
    except FileNotFoundError as e:
        print(f'Gammas file missing!')
        continue

    outdir = f'gamma_v_radii_{runtype}/'
    save_npy(np.vstack((radii, gamL)).T, f'gamL_{gamma_npy_suffix}-{nn}', npydir=outdir)
    save_npy(np.vstack((radii, gamR)).T, f'gamR_{gamma_npy_suffix}-{nn}', npydir=outdir)
    print('Done!')
