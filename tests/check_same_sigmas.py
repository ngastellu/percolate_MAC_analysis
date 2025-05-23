
#!/usr/bin/env python

import numpy as np


kB = 8.617e-5
e2C = 1.602177e-19 # elementary charge to Coulomb
w0 = 1e15
T = 300

# This factor combines the hop frequency with the unit conversion (to yield conductivity in siemens)
# w0 is chosen such that the final result matches the results from the AMC paper.
conv_factor = e2C*w0

sigmasdir = '/Users/nico/Desktop/simulation_outputs/percolation/sigmas_v_T/'
# motypes = ['kBTlo_dcut500','virtual','kBThi']
motypes = ['lo','virtual_w_HOMO','hi']

# temps = np.arange(40,440,10)[14:]

r_maxs = ['18.03', '136.47', '199.33']
ensemble_lbls = ['sAMC-500', 'sAMC-q400', 'sAMC-300']
structypes=['40x40', 'tempdot6', 'tempdot5']
for st, st2, rmax in zip(structypes,ensemble_lbls,r_maxs):
    print(f'\n---------- {st} ----------')
    for mt in motypes:
        try:
            # out1 = np.load(f'/Users/nico/Desktop/simulation_outputs/percolation/sigmas_v_T/sigma_v_T_w_err_{st}_rmax_{rmax}_{mt}.npy').T
            out1 = np.load(f'/Users/nico/Desktop/simulation_outputs/percolation/sigmas_v_T/sandbox_run/{st2}/sigma_v_T_w_err_production_{mt}.npy')[0,1:]
        except FileNotFoundError as e:
            out1 = np.load(f'/Users/nico/Desktop/simulation_outputs/percolation/sigmas_v_T/sigma_v_T_w_err_{st}_rmax_{rmax}_sites_gammas_{mt}.npy')
            out1 = out1[out1[:,0]==T][0] # in this case, out1 contains data for multiple temperatures
            # print(out1)
        if st == '40x40':
            out2 = np.load(f'/Users/nico/Desktop/simulation_outputs/percolation/sigmas_v_T/sigma_v_T_w_err_rerun2_rmax_{rmax}_{mt}.npy')[0,1:] * conv_factor/(kB*T)
        else:
            out2 = np.load(f'/Users/nico/Desktop/simulation_outputs/percolation/sigmas_v_T/sigma_v_T_w_err_{st}_rerun_rmax_{rmax}_{mt}.npy')[0,1:] * conv_factor/(kB*T)        
        
        match = np.all(out1 == out2)
        print(f'\nOutputs for {mt} match: ', match)
        if not match:
            print('out1: ', out1)
            print('out2: ', out2)
            print('diff: ', np.abs(out1-out2))

