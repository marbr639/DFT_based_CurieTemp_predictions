#!/usr/bin/env python3

import os
import numpy as np
import bz2
import subprocess
import sys

if len(sys.argv) >= 7:
    OSZICAR_GS = str(sys.argv[1])
    OSZICAR_DLM = str(sys.argv[2])
    OUTCAR_GS = str(sys.argv[3])
    OUTCAR_DLM = str(sys.argv[4])
    POSCAR_GS = str(sys.argv[5])
    POSCAR_DLM = str(sys.argv[6])
else: 
    print("Write the names of the following files: \n - GS OSZICAR, DLM OSZICAR, GS OUTCAR, DLM OUTCAR, GS POSCAR, and DLM POSCAR. \nAll files must be in the same place as this script. \nMake sure the script gather_mag_final.sh is in the same place as this script.")
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
    bash_script = ["./"+script, str(num_atoms), outcar]

    process = subprocess.Popen(bash_script, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout, stderr = process.communicate()

    magmoms = []
    if process.returncode == 0:
        with open('magmom_size', 'r') as magmom_file:
            output_lines = magmom_file.readlines()

        for line in output_lines:
            magmoms.append(line.strip()) 
    else:
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



kB = 8.617333262*10**(-5)



def main():

    energies = open('energy_entropy', 'w+')

    formula = parse_poscar(POSCAR_DLM)

    oszicar_gs_filename = OSZICAR_GS
    oszicar_lambda_filename = OSZICAR_DLM

    _, ext = os.path.splitext(oszicar_gs_filename)
    if ext == '.bz2':
        with bz2.open(oszicar_gs_filename, 'rt') as f:
            oszicar_gs = f.readlines()

    else:
         with open(oszicar_gs_filename, 'r') as f:
            oszicar_gs = f.readlines() 

    _, ext = os.path.splitext(oszicar_lambda_filename)
    if ext == '.bz2':       
        with bz2.open(oszicar_lambda_filename, 'rt') as f:
            oszicar_lambda = f.readlines()

    else: 
        with open(oszicar_lambda_filename, 'r') as f:
            oszicar_lambda = f.readlines()
   
            
 

    atom_types_gs, atom_counts_gs = read_poscar(POSCAR_GS)
    atom_types_lambda, atom_counts_lambda = read_poscar(POSCAR_DLM)


    # Create a list of atom types repeated according to their counts
    complete_atom_types_gs = [atom_type for atom_type, count in zip(atom_types_gs, atom_counts_gs) for _ in range(count)]
    complete_atom_types_lambda = [atom_type for atom_type, count in zip(atom_types_lambda, atom_counts_lambda) for _ in range(count)]
        
    tot_atoms_gs = sum(atom_counts_gs)
    tot_atoms_lambda = sum(atom_counts_lambda)

    tot_magmom_sizes_gs = np.array(get_magmoms(OUTCAR_GS, 'gather_mag_final.sh', tot_atoms_gs), dtype=float)

    tot_magmom_sizes_lambda = np.array(get_magmoms(OUTCAR_DLM, 'gather_mag_final.sh', tot_atoms_lambda), dtype=float)
        


    idx_constr_moments = []

    for idx, line in enumerate(oszicar_lambda):
        if 'lambda*MW_perp' in line:
            for i in range(idx + 1, len(oszicar_lambda)):
                idx_constr = oszicar_lambda[i].split()
                idx_constr_moments.append(int(idx_constr[0])-1)


    type_constr_moments = [complete_atom_types_lambda[idx] for idx in idx_constr_moments]
    idx_constr_moments_gs = [i for i, atom_type in enumerate(complete_atom_types_gs) if atom_type in type_constr_moments]

    constr_moments_gs = tot_magmom_sizes_gs[idx_constr_moments_gs]
    constr_moments_lambda = [tot_magmom_sizes_lambda[idx] for idx in idx_constr_moments]

    num_atoms_gs = len(constr_moments_gs)
    num_atoms_lambda = len(constr_moments_lambda)


    print('{} has  {} constrained magnetic moments. Number of > 0.75 $\mu_B$ moments in GS is {}'.format(formula,num_atoms_lambda, num_atoms_gs ))
    print('The constrained moments of {} in the DLM state are: {}'.format(formula, constr_moments_lambda))

    energy_gs=0
    for idx,line in enumerate(oszicar_gs):
        if 'F=' in line:
            energy_gs = float(line.split()[4])

    energy_lambda=0
    for idx,line in enumerate(oszicar_lambda):
        if 'F=' in line:
            energy_lambda = float(line.split()[4])

        
    mean_ln = sum([np.log(float(m)+1) for m in constr_moments_lambda])/len(constr_moments_lambda)
    S = kB*mean_ln      

    diff_energy = (energy_lambda / num_atoms_lambda) - (energy_gs / num_atoms_gs)
    diff_energy_S = ((energy_lambda/num_atoms_lambda) - (energy_gs/num_atoms_gs))/S
    energies.write('{}  {} {} {}\n'.format(formula, diff_energy, S, diff_energy_S))

    energies.close()

if __name__ == "__main__":
    main()
