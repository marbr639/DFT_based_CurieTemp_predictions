import fileinput
import os
from shutil import copyfile
import subprocess
import numpy as np  # Import numpy as np for consistency
import random
from distutils.dir_util import copy_tree
from distutils.dir_util import remove_tree
import shutil
import pandas
import ase.io
import ase.io.vasp
from ase.neighborlist import neighbor_list
import time
import re
import sys
import math
import bz2
from get_SRO import Structure,  calc_neighbors, readInitialStructure, set_magmoms

if len(sys.argv) >= 3:
    POSCAR = str(sys.argv[1])
    INCAR = str(sys.argv[2])
else: 
    print("Write the names of the following files in this order: \n - DLM POSCAR and DLM INCAR. \nAll files must be in the same place as this script. \nMake sure the script gather_mag_final.sh is in the same place as this script.")
    sys.exit(1)

def parse_poscar(filename):
    _, ext = os.path.splitext(filename)
    if ext == '.bz2':
        with bz2.open(filename, 'rt') as f:
            lines = f.readlines()
    else:
        with open(filename, 'r') as f:
            lines = f.readlines()
    element = lines[5].split()
    distr = lines[6].split()
    formula = ''.join([e + d for e, d in zip(element, distr)])
    return formula


def main():

    NN_file = open('Number_of_neighbors', 'w+')
    formula = parse_poscar(POSCAR)
    print('Working on ', formula)
    structure = Structure()
    readInitialStructure(structure, 'POSCAR', '')

    set_magmoms(structure, '', 'INCAR', INCAR)
    NNList = calc_neighbors(structure, POSCAR, tol=0.1)
    natoms = len(NNList)
    print('Total number of atoms: ', natoms)
    nneighbors= []
    for ind, atom in enumerate(NNList):
        print('Number of nearest magnetic neighbors for atom',ind,' is ', len(atom))
        nneighbors.append(len(atom))

    mean_NN = np.mean(nneighbors)
    NN_file.write('{} {}\n'.format(formula, mean_NN))



if __name__ == "__main__":
    main()
