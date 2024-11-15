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


def parse_poscar(filename):
    with open(filename, 'r') as f:
        lines = f.readlines()
    element = lines[5].split()
    distr = lines[6].split()
    formula = ''.join([e + d for e, d in zip(element, distr)])
    return formula

def read_poscar(filename):
    _, ext = os.path.splitext(filename)
    if ext == '.bz2':
        with bz2.open(filename, 'rt') as f:
            lines = f.readlines()
    else:
        with open(filename, 'r') as f:
            lines = f.readlines()

    # Extract atom types and atom counts based on POSCAR format
    atom_types = lines[5].split()
    atom_counts = list(map(int, lines[6].split()))

    return atom_types, atom_counts

def get_magmoms(outcar,script, num_atoms):
    # Define your Bash script as a list of command and arguments
    bash_script = ["./"+script, str(num_atoms), outcar]

    # Run the Bash script and capture its output
    process = subprocess.Popen(bash_script, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout, stderr = process.communicate()

    magmoms = []
    if process.returncode == 0:
        # Read the contents of magmom_size.txt into a list
        with open('magmom_size', 'r') as magmom_file:
            output_lines = magmom_file.readlines()
    
        # Process output_lines as needed
        for line in output_lines:
            magmoms.append(line.strip())  # Strip any leading/trailing whitespace
    else:
        # If there was an error, you can print the error message
        print("Error:", result.stderr)

    return magmoms

# Read data from files
def read_curie(filename):
    data = {}
    with open(filename, 'r') as file:
        for line in file:
            columns = line.strip().split()
            data[columns[0]] = float(columns[1])
    return data

Curie_temp_data = read_curie("Curie_temp")



kB = 8.617333262*10**(-5)


def main():
    folders = [f for f in os.listdir() if os.path.isdir(f)]

    NN_file = open('Number_of_neighbors', 'w+')

    for folder in folders:
        poscar_filename = os.path.join(folder, 'POSCAR')
        if not os.path.exists(poscar_filename):
            print('No ', poscar_filename, ' found!')
            continue
        formula = parse_poscar(poscar_filename)
        print('Working on ', formula)
        structure = Structure()
        readInitialStructure(structure, 'POSCAR', folder)

        run_folder = [f for f in os.listdir(folder) if f.startswith('ht.run')]
        if not run_folder:
            continue
        run_folder = run_folder[0]

        incar_filename = os.path.join(folder, run_folder, 'INCAR')
        if not os.path.exists(incar_filename):
            print('No ', incar_filename, ' found!')
            continue
        set_magmoms(structure, '', 'INCAR', incar_filename)
        NNList = calc_neighbors(structure, poscar_filename, tol=0.1)
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
