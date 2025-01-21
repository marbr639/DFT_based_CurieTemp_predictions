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

class Structure():
    """
    Class Structure
    A representation of an atomic structure, including information about the constituent atoms, the initial 
    structure, and methods to analyze magnetic properties.

    Attributes:
    - atomsList: A list of atoms in the structure. Defaults to an empty list if not provided.
    - initial_Structure: An ASE Atoms object representing the initial atomic structure, which can include 
      information such as atomic positions, lattice vectors, and magnetic moments.

    Methods:
    - __init__(atomsList=None): Initializes the structure with an optional list of atoms.
    - __getitem__(i_atoms): Retrieves the atom at the specified index in atomsList.
    - __len__(): Returns the total number of atoms in the structure.
    - append_atoms(atoms): Adds one or more atoms to the atomsList. Supports appending a single atom or a list of atoms.
    - set_initial_Structure(Initial_Structure): Sets the initial_Structure attribute with an ASE Atoms object.
    - get_initial_Structure(): Returns the initial_Structure attribute.
    - get_magnetic_Structure(*args): Identifies magnetic atoms based on their magnetic moments. 
      - If magnetic moments are vectors (3D), their magnitudes are calculated.
      - Atoms with magnetic moments larger than 0.75 are considered magnetic.
      - Returns:
        - MagneticStructure: An ASE Atoms object containing only the magnetic atoms.
        - idx: Indices of the magnetic atoms in the initial structure.

    Notes:
    - Magnetic atoms are assigned unique tags starting from 1 for identification.
    - This class provides basic functionality for storing, modifying, and analyzing atomic structures,
      with a specific focus on identifying magnetic atoms based on their magnetic moments.
    """

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

def calculate_distances(indices, poscar_path):
    """
    Calculates the pairwise distances between specified atoms in a structure, considering periodic boundary conditions.

    This function computes the shortest distances between pairs of atoms specified by their indices. 
    Distances are calculated using the lattice vectors from a POSCAR file, accounting for periodic 
    boundary conditions by considering translations in all directions. The function also identifies 
    the relative vectors corresponding to the minimum distances.

    Parameters:
    - indices: A list of atom indices for which distances need to be calculated.
    - poscar_path: Path to the POSCAR file containing lattice vectors and atomic positions.

    Returns:
    - A dictionary where each key is an atom index, and the value is a sorted list of tuples. 
      Each tuple contains:
      - The index of a neighboring atom.
      - The shortest distance to that neighbor (rounded to six decimal places).
      - The relative vector corresponding to the minimum distance.
    """
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

            distances_dict[i].append((indices[j], round(min_dist,6), min_dist_vector))
            distances_dict[j].append((indices[i], round(min_dist,6), min_dist_vector))

    # Sort distances for each atom index
    for atom_index, neighbors in distances_dict.items():
        distances_dict[atom_index] = sorted(neighbors, key=lambda x: x[1])

    return distances_dict

def calc_neighbors(structure, poscar_path, tol):
    """
    Calculates the nearest neighbors for magnetic atoms in a given structure.

    This function identifies magnetic atoms in the structure based on their magnetic moments 
    and computes their nearest neighbors within a specified tolerance. The positions of the 
    magnetic atoms are used to calculate pairwise distances, and neighbors are determined 
    based on a distance threshold, which is the minimum interatomic distance plus the 
    tolerance.

    Parameters:
    - structure: An object representing the structure, including magnetic properties.
    - poscar_path: Path to the POSCAR file containing atomic positions and lattice information.
    - tol: Tolerance value to determine the neighbor cutoff distance.

    Returns:
    - A list of nearest neighbors for each magnetic atom, where each neighbor is represented 
      by its index and relative vector.
    - Returns "NOMAGMOM" if no significant magnetic moments are detected (default threshold: 0.75).
    """
    # Get magnetic sites of the Structure
    MagneticAtoms = structure.get_magnetic_Structure()[0]
    indices_magnetic_atoms = structure.get_magnetic_Structure()[1][0]

    if len(indices_magnetic_atoms)==0:
        # If all magnetic moments are smaller than 0.75 (update!)
        return "NOMAGMOM"
    else:
        magmoms = [atom.magmom for atom in MagneticAtoms]
        positions = [atom.position for atom in MagneticAtoms]
        complete_dictionary = calculate_distances(list(indices_magnetic_atoms), poscar_path)
        # it is enough to look at the first atom to pick out the maximum distance for our neighbor sphere
        dist_lst = [i[1] for i in complete_dictionary[0]]

        # Calculate the maximum distance
        maxDist = min(dist_lst)+tol


        initial_NNList = []

        for i in range(len(complete_dictionary)):
            tmp_dist = 0
            ind = 0
            initial_NNList.append([])

            while tmp_dist < maxDist and ind < len(complete_dictionary[i]) - 1:
                vector = complete_dictionary[i][ind][2]
                initial_NNList[i].append([complete_dictionary[i][ind][0], vector[0], vector[1], vector[2]])
                tmp_dist = complete_dictionary[i][ind + 1][1]
                ind += 1

        return initial_NNList


def readInitialStructure(structure, filename, base_path):
    """
    Reads and initializes the atomic structure from a VASP POSCAR/CONTCAR file.

    This function reads the structure file specified by the filename and base path, 
    normalizes the lattice parameters by dividing by the absolute value of the scaling 
    factor (read from the second line of the file), and updates the structure's cell 
    accordingly. It ensures that the atomic positions are scaled consistently with the 
    updated cell dimensions.

    Parameters:
    - structure: An object to hold the initial structure data.
    - filename: Name of the VASP POSCAR/CONTCAR file to read.
    - base_path: Base directory where the structure file is located.

    Returns:
    None
    """
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
    """
    Sets the magnetic configuration for a given structure based on the specified input type. 
    The magnetic moments can be initialized using one of the following methods:

    1. INCAR file ('INCAR'): Extracts magnetic moments based on the M_CONSTR line in the INCAR file (i.e. unconstrained magnetic moments have value zero!). 
       Handles both collinear and non-collinear setups.

    2. External file ('magmoms'): Reads an external file containing an array of magnetic moments (x-, y-, and z-components of each magnetic moment must be included).

    3. Direct array ('array'): Uses a provided array of magnetic moments directly.

     For each atom in the structure, the magnetic moments are assigned after ensuring the input format is valid.

    Parameters:
    - structure: The structure object whose magnetic moments will be modified.
    - base_path: The base directory where input files are located.
    - input_type: Specifies the source of the magnetic moment data ('INCAR', 'magmoms', or 'array').
    - input_config: Path to the input file or the magnetic moment array.

    Returns:
    None
    """
    if input_type == 'INCAR':

        file_h = os.path.join(base_path, input_config)

        if not os.path.isfile(os.path.join(base_path, input_config)):
            print("No INCAR found!")
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


def spin_spin_correlation(structure, initial_Structure, initial_NNList):
    """

    Calculates the average spin-spin correlation (short-range order, SRO) for a magnetic structure based 
    on its nearest-neighbor magnetic interactions.

    Parameters:
    - structure: A Structure object containing the full atomic configuration, including magnetic atoms.
    - initial_Structure: The initial atomic structure (ASE Atoms object) with magnetic moment information.
    - initial_NNList: A list of nearest-neighbor information for each magnetic atom. 
      - Each element corresponds to a magnetic atom and contains information about its neighboring atoms,
        including their indices and relative displacement vectors (dx, dy, dz).

    Returns:
    - SRO: The average spin-spin correlation, calculated as the normalized dot product of unit vectors 
      representing the magnetic moments of neighboring magnetic atoms.

    Steps:
    1. Extract the magnetic atoms and their indices from the structure using `get_magnetic_Structure`.
    2. Loop over each magnetic atom and its list of neighbors from `initial_NNList`.
    3. For each pair (magnetic atom and its neighbor):
       - Normalize their magnetic moments to obtain unit vectors.
       - Compute the dot product of the unit vectors to determine the spin correlation.
       - Normalize the spin correlation by the number of neighbors and add to the total SRO.
    4. Average the SRO by dividing by the total number of magnetic atoms.
    5. Return the computed SRO.

    Notes:
    - Magnetic moments are normalized to ensure that the calculation uses unit vectors, which avoids 
      bias due to the magnitude of magnetic moments.
    - The contribution of each pair to the SRO is weighted equally among the neighbors.
    - This function provides insight into the local magnetic order in the system, which can be used 
      to assess the degree of alignment or disorder in the magnetic structure.
    """

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
