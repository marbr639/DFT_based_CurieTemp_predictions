import fileinput
import os
from shutil import copyfile
import subprocess
import numpy as np  # Import numpy as np for consistency
import random
from distutils.dir_util import copy_tree
from distutils.dir_util import remove_tree
import shutil
import ase.io
import ase.io.vasp
from ase.neighborlist import neighbor_list
import time
import re
import sys
import math
from get_SRO import Structure, calc_neighbors, readInitialStructure, set_magmoms,spin_spin_correlation



if len(sys.argv) >= 3:
    input_lst = [float(i) for i in sys.argv[1:]]
    length=len(input_lst)
    if (length % 2) == 0:
        nAtoms = input_lst[:int(length/2)]
        averageSizes = input_lst[int(length/2):]




mag_configs = []
NUNCONSTR = 0

for conf in range(0,50):
    c = {'config MAGMOM': [], 'config M_CONSTR': [], 'ssCorr': None}
    for n, mag in zip(nAtoms, averageSizes):
        for i in range(0,int(n)):
            point=np.array((0.0,0.0,0.0));


            x=np.random.rand(4)*2-1;
            x2=np.square(x)

            while (np.sum(x2)>=1.0):
                x=np.random.rand(4)*2-1;
                x2=np.square(x);

            point[0]=2*(x[1]*x[3]+x[0]*x[2])/np.sum(x2);
            point[1]=2*(x[2]*x[3]-x[0]*x[1])/np.sum(x2);
            point[2]=(x2[0]+x2[3]-x2[1]-x2[2])/np.sum(x2);
            if mag < 0.75:
                c['config M_CONSTR'].append([0.0, 0.0, 0.0])
                NUNCONSTR += 1
                c['config MAGMOM'].append([point[0]*mag, point[1]*mag, point[2]*mag])
            else:
                c['config M_CONSTR'].append([point[0]*mag, point[1]*mag, point[2]*mag])
                c['config MAGMOM'].append([point[0]*mag, point[1]*mag, point[2]*mag])
    mag_configs.append(c)


for ind, c in enumerate(mag_configs):

    structure = Structure()
    readInitialStructure(structure, 'POSCAR', '')

    # Taking M_CONSTR since here the small moments are set to zero and will not
    # contribute to the determination of smalles SRO

    set_magmoms(structure, '', 'array', c['config M_CONSTR'])

    NNList = calc_neighbors(structure, 'POSCAR', tol=0.1)

    if NNList == "NOMAGMOM":
        ssCorr = 1
        mag_configs = []
    else:
        ssCorr = spin_spin_correlation(structure, structure.initial_Structure, NNList)
        mag_configs[ind]['ssCorr'] = ssCorr



best_config_MAGMOM = []
best_config_M_CONSTR = []
MAGMOM_str = ''
M_CONSTR_str = ''
smallest_ssCorr = None

for i in mag_configs:


    if smallest_ssCorr == None:
        smallest_ssCorr = abs(i['ssCorr'])
        best_config_MAGMOM = i['config MAGMOM']
        best_config_M_CONSTR = i['config M_CONSTR']

    else:
        if smallest_ssCorr > abs(i['ssCorr']):
            smallest_ssCorr = abs(i['ssCorr'])
            best_config_MAGMOM = i['config MAGMOM']
            best_config_M_CONSTR = i['config M_CONSTR']



for m, c in zip(best_config_MAGMOM, best_config_M_CONSTR):

    MAGMOM_str = MAGMOM_str + '{} {} {} '.format(m[0], m[1], m[2])
    M_CONSTR_str = M_CONSTR_str + '{} {} {} '.format(c[0], c[1], c[2])

NUNCONSTR_str = str(NUNCONSTR/50)

print(f'MAGMOM = {MAGMOM_str}')
print(f'M_CONSTR = {M_CONSTR_str}')
