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
from ase.neighborlist import neighbor_list, build_neighbor_list, first_neighbors
import time
import re
import sys
import math
from collections import defaultdict

# Define a Structure class for handling atomic structures and magnetic moments
class Structure():

    def __init__(self, atomsList=None):
        if atomsList is None:
            atomsList = []
        self.atomsList = atomsList
        self.initial_Structure=None

    def __getitem__(self, i_atoms):
        return self.atomsList[i_atoms]

    def __len__(self):
        return len(self.atomsList)

    def append_atoms(self, atoms):
        if isinstance(atoms, list):
            for i_atoms in range(0,len(atoms)):
                self.atomsList.append(atoms)
        else:
            self.atomsList.append(atoms)

    def set_initial_Structure(self, Initial_Structure):
        self.initial_Structure = Initial_Structure  # Fix typo: newInital_Structure -> Initial_Structure
        return

    def get_initial_Structure(self):
        return self.initial_Structure

    def get_magnetic_Structure(self, *args):
        # Identify magnetic atoms based on their magnetic moments
        if self.initial_Structure.get_initial_magnetic_moments().shape[1] == 3:
            mag = np.linalg.norm(self.initial_Structure.get_initial_magnetic_moments(), axis=1)
        else:
            mag = self.initial_Structure.get_initial_magnetic_moments()

        idx = np.where(mag > 0.75)


        MagneticStructure = self.initial_Structure[idx[0]]

        # Assign unique tags to magnetic atoms
        for i_atoms, atom in enumerate(MagneticStructure):
            MagneticStructure[i_atoms].tag = i_atoms + 1



        return MagneticStructure, idx

    # Other methods...


def calculate_distances(indices, poscar_path):
    poscar = open(poscar_path)
    all_lines_poscar = poscar.readlines()
    poscar.close()
    lp = float(all_lines_poscar[1])
    if lp < 0:
        lp = 1
    a = list(map(float, all_lines_poscar[2].split()))
    b = list(map(float, all_lines_poscar[3].split()))
    c = list(map(float, all_lines_poscar[4].split()))
    dist_lst_total = []
    pos_direct = []
    Natoms = sum(list(map(int, all_lines_poscar[6].split())))
    for i in range(0, Natoms):
        pos_direct.append(list(map(float, all_lines_poscar[i + 8].split())))


    magnetic_positions = [pos_direct[i] for i in indices]
    distances_dict = defaultdict(list)
    for i in range(len(magnetic_positions)):
        for j in range(i + 1, len(magnetic_positions)):
            pos_i = np.array(magnetic_positions[i])
            pos_j = np.array(magnetic_positions[j])

            min_dist = float('inf')
            min_dist_vector = None
            for k in [1, 0, -1]:
                tmpx = pos_i[0] - (pos_j[0] + k)
                for l in [1, 0, -1]:
                    tmpy = pos_i[1] - (pos_j[1] + l)
                    for m in [1, 0, -1]:
                        tmpz = pos_i[2] - (pos_j[2] + m)

                        distx = (a[0] ** 2 + a[1] ** 2 + a[2] ** 2) * tmpx ** 2 * lp ** 2
                        disty = (b[0] ** 2 + b[1] ** 2 + b[2] ** 2) * tmpy ** 2 * lp ** 2
                        distz = (c[0] ** 2 + c[1] ** 2 + c[2] ** 2) * tmpz ** 2 * lp ** 2
                        distxy = 2 * tmpx * tmpy * (a[0] * b[0] + a[1] * b[1] + a[2] * b[2]) * lp ** 2
                        distxz = 2 * tmpx * tmpz * (a[0] * c[0] + a[1] * c[1] + a[2] * c[2]) * lp ** 2
                        distyz = 2 * tmpy * tmpz * (b[0] * c[0] + b[1] * c[1] + b[2] * c[2]) * lp ** 2
                        dist = np.sqrt(distx + disty + distz + distxy + distxz + distyz)

                        if dist < min_dist:
                            min_dist = dist
                            min_dist_vector = (tmpx, tmpy, tmpz)

            #print(i, indices[i])
            #print(j, indices[j])
            distances_dict[i].append((indices[j], round(min_dist,6), min_dist_vector))
            distances_dict[j].append((indices[i], round(min_dist,6), min_dist_vector))

    # Sort distances for each atom index
    for atom_index, neighbors in distances_dict.items():
        distances_dict[atom_index] = sorted(neighbors, key=lambda x: x[1])

    return distances_dict





# Define a function to calculate the nearest neighbors of magnetic atoms
def calc_neighbors(structure, poscar_path, tol):
    # Get magnetic sites of the Structure
    MagneticAtoms = structure.get_magnetic_Structure()[0]
    indices_magnetic_atoms = structure.get_magnetic_Structure()[1][0]

    if len(indices_magnetic_atoms)==0:
        # If all magnetic moments are smaller than 0.75 (update!)
        return "NOMAGMOM"
    else:
        magmoms = [atom.magmom for atom in MagneticAtoms]
        positions = [atom.position for atom in MagneticAtoms]
        test_dist_lst = calculate_distances(list(indices_magnetic_atoms), poscar_path)
        dist_lst = [i[1] for i in test_dist_lst[0]]






        # Calculate the maximum distance
        maxDist = min(dist_lst)+tol


        initial_NNList = []

        for i in range(len(test_dist_lst)):
            tmp_dist = 0
            ind = 0
            initial_NNList.append([])

            while tmp_dist < maxDist and ind < len(test_dist_lst[i]) - 1:
                vector = test_dist_lst[i][ind][2]
                initial_NNList[i].append([test_dist_lst[i][ind][0], vector[0], vector[1], vector[2]])
                tmp_dist = test_dist_lst[i][ind + 1][1]
                ind += 1

        return initial_NNList

# Define a function to read the initial atomic structure from a file
def readInitialStructure(structure, filename, base_path):
    file_h = os.path.join(base_path, filename)
    tmp_f = open(file_h)
    lines = tmp_f.readlines()
    tmp_f.close()
    lat_param = abs(float(lines[1]))

    structure.initial_Structure = ase.io.vasp.read_vasp(file_h)

    cell=structure.initial_Structure.cell.lengths()
    structure.initial_Structure.set_cell(cell/lat_param, scale_atoms=True)

    return

def set_magmoms(structure, base_path, input_type, input_config):

    if input_type == 'INCAR':

        file_h = os.path.join(base_path, input_config)

        if not os.path.isfile(os.path.join(base_path, input_config)):
            print("No INCAR_INIT found! The file should have been created in the first run!")
            exit()

        for line in fileinput.input(file_h, inplace=False):
            if line.strip().startswith('LNONCOLLINEAR'):

                if '.TRUE.' in line:
                    N = 3
                    print('NONCOLLINEAR')
                else:
                    N = 1
                    print('COLLINEAR')
            if line.strip().startswith('MAGMOM'):
                allMoments = np.reshape(re.findall(r'[-+]?\d*\.\d+(?:[eE][-+]?\d+)?', line), (-1, N))

    elif input_type == 'magmoms':

        file_h = os.path.join(base_path, input_config)


        if not os.path.isfile(os.path.join(base_path, input_config)):
            print("No magmom-file found!")
            exit()

        allMoments = []


        tmp = open(file_h, 'r')
        lines = tmp.readlines()
        tmp.close()
        for line in lines:
            allMoments.append([line.split()[0], line.split()[1], line.split()[2]])


    elif input_type == 'array':
        allMoments = input_config

    for i_atom in range(0, len(structure.initial_Structure)):
        
        # Ensure that allMoments has the correct shape (3 elements)
        if len(allMoments[i_atom]) != 3:
            print(f"Error: Incorrect number of elements in MAGMOM for atom {i_atom + 1}")
            exit()

        # Assign the magnetic moments to the atom
        structure.initial_Structure[i_atom].magmom = allMoments[i_atom]

    return


# Define a function to calculate spin-spin correlations
def spin_spin_correlation(structure, initial_Structure, initial_NNList):
    # Initialize an array to store spin-spin correlations for each pair
    SRO=0
    natoms=len(initial_NNList)
    MagneticAtoms = structure.get_magnetic_Structure()[0]
    indices_magnetic_atoms = structure.get_magnetic_Structure()[1][0]

    # Iterate over each magnetic atom
    for i, atoms in enumerate(initial_NNList):

        # Calculate spin-spin correlation with neighboring magnetic atoms
        for j, neighbor_info in enumerate(atoms):

            neighbor_index, dx, dy, dz = neighbor_info

            # Access magnetic moments for the current atom and its neighbor
            mag_moment_i = initial_Structure[indices_magnetic_atoms[i]].magmom
            e_i = mag_moment_i/np.linalg.norm(mag_moment_i)
            mag_moment_j = initial_Structure[int(neighbor_index)].magmom
            e_j = mag_moment_j/np.linalg.norm(mag_moment_j)


            # Calculate the dot product of magnetic moments
            spin_corr_ij = (1/len(atoms))*np.dot(e_i, e_j)
            SRO +=spin_corr_ij


    SRO= SRO/natoms

    return SRO
