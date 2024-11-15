#!/usr/bin/env python3

import os
import numpy as np
import bz2
import subprocess

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

    energies = open('energies', 'w+')

    desired_atom_types = ['Fe', 'Mn', 'Gd', 'Co', 'Cr']

    for folder in folders:
        poscar_filename = os.path.join(folder, 'POSCAR')
        if not os.path.exists(poscar_filename):
            continue

        formula = parse_poscar(poscar_filename)

        run_folder = [f for f in os.listdir(folder) if f.startswith('ht.run')]
        if not run_folder:
            continue
        run_folder = run_folder[0]

        oszicar_gs_filename = os.path.join(folder, run_folder, 'GS_energies.bz2')
        oszicar_lambda_filename = os.path.join(folder, run_folder, 'OSZICAR.lambda-final.bz2')

        if not (os.path.exists(oszicar_gs_filename) and os.path.exists(oszicar_lambda_filename)):
            continue

        with bz2.open(oszicar_gs_filename, 'rt') as f:
            oszicar_gs = f.readlines()
        with bz2.open(oszicar_lambda_filename, 'rt') as f:
            oszicar_lambda = f.readlines()

        atoms_lambda_poscar_filename = os.path.join(folder, run_folder, 'POSCAR.bz2')
        atoms_GS_poscar_filename = os.path.join(folder, 'Primitive.vasp')

        atom_types_gs, atom_counts_gs = read_poscar(atoms_GS_poscar_filename)
        atom_types_lambda, atom_counts_lambda = read_poscar(atoms_lambda_poscar_filename)


        # Create a list of atom types repeated according to their counts
        complete_atom_types_gs = [atom_type for atom_type, count in zip(atom_types_gs, atom_counts_gs) for _ in range(count)]
        complete_atom_types_lambda = [atom_type for atom_type, count in zip(atom_types_lambda, atom_counts_lambda) for _ in range(count)]


        
        tot_atoms_gs = sum(atom_counts_gs)
        tot_atoms_lambda = sum(atom_counts_lambda)

        outcar_path_gs = os.path.join(folder, run_folder, 'OUTCAR.cleaned.GS-final.bz2')
        if not os.path.exists(outcar_path_gs):
            print('No OUTCAR.cleaned.GS-final.bz2 file found!')

        tot_magmom_sizes_gs = np.array(get_magmoms(outcar_path_gs, 'gather_mag_final.sh', tot_atoms_gs), dtype=float)

        outcar_path_lambda = os.path.join(folder, run_folder, 'OUTCAR.cleaned.lambda-final.bz2')
        if not os.path.exists(outcar_path_lambda):
             outcar_path_lambda = os.path.join(folder, run_folder, 'OUTCAR.bz2')
             
             
        

        tot_magmom_sizes_lambda = np.array(get_magmoms(outcar_path_lambda, 'gather_mag_final.sh', tot_atoms_lambda), dtype=float)
        


        #idx_big_moments_gs = np.where(tot_magmom_sizes_gs > 0.75)[0] 


       # big_moments_gs = tot_magmom_sizes_gs[idx_big_moments_gs]


        idx_constr_moments = []
        tmp = bz2.open(folder + '/' + run_folder + '/OSZICAR.lambda-final.bz2', 'r')
        OSZICAR = tmp.readlines()
        tmp.close()

        for idx, line in enumerate(OSZICAR):
            if b'lambda*MW_perp' in line:
                for i in range(idx + 1, len(OSZICAR)):
                    idx_constr = OSZICAR[i].split()
                    idx_constr_moments.append(int(idx_constr[0])-1)


        #print(complete_atom_types_lambda)

        type_constr_moments = [complete_atom_types_lambda[idx] for idx in idx_constr_moments]
        

        idx_big_moments_gs = [i for i, atom_type in enumerate(complete_atom_types_gs) if atom_type in type_constr_moments]


        big_moments_gs = tot_magmom_sizes_gs[idx_big_moments_gs]
        #type_big_moments = [complete_atom_types_gs[idx] for idx in idx_big_moments_gs]
        
        #type_constr_moments = [complete_atom_types_lambda[idx] for idx in idx_constr_moments]
       
        #idx_big_moments_lambda = np.where(complete_atom_types_lambda == type_constr_moments)[0]

        # Get indices in complete_atom_types_lambda where the entry matches any entry in type_big_moments
        #idx_big_moments_lambda = [index for index, atom_type in enumerate(complete_atom_types_lambda) if atom_type in type_big_moments]
        #idx_big_moments_gs = np.where(tot_magmom_sizes_lambda > 0.75)[0] 

        #big_moments_lambda = [tot_magmom_sizes_lambda[idx] for idx in idx_big_moments_lambda]
        
        big_moments_lambda = [tot_magmom_sizes_lambda[idx] for idx in idx_constr_moments]

        
        

        num_atoms_gs = len(big_moments_gs)
        num_atoms_lambda = len(big_moments_lambda)


        print('{} has  {} constrained magnetic moments. Number of big moments in GS is {}'.format(formula,num_atoms_lambda, num_atoms_gs ))
        print('The big moments of {} in the DLM state are: {}'.format(formula, big_moments_lambda))
        #print(formula, num_atoms_gs, num_atoms_lambda)
        #print(big_moments_gs, big_moments_lambda)

       
        #for atom_type_gs, atom_count_gs in zip(atom_types_gs, atom_counts_gs):
        #    if atom_type_gs in desired_atom_types:
        #        num_atoms_gs += atom_count_gs


        #for atom_type_lambda, atom_count_lambda in zip(atom_types_lambda, atom_counts_lambda):
        #    if atom_type_lambda in desired_atom_types:
        #        num_atoms_lambda += atom_count_lambda
        #    tot_atoms_lambda += atom_count_lambda

        # Extracting energies and atom counts

        energy_gs_lst = [float(item.strip()) for item in oszicar_gs if item.strip()]
        #print(formula)
        #print(energy_gs_lst)
        energy_gs = min(energy_gs_lst)



        energy_lambda=0
        for idx,line in enumerate(oszicar_lambda):
            if 'F=' in line:
                energy_lambda = float(line.split()[2])



        #mean_mag = sum([float(m)/len(big_moments_lambda) for m in big_moments_lambda])
        #lnW=np.log(mean_mag + 1)
        #S= kB*lnW
        
        mean_ln = sum([np.log(float(m)+1) for m in big_moments_lambda])/len(big_moments_lambda)
        S = kB*mean_ln
        

        # Calculating and writing energies
        # Add your energy calculation and writing logic here

        diff_energy = (energy_lambda / num_atoms_lambda) - (energy_gs / num_atoms_gs)
        diff_energy_S = ((energy_lambda/num_atoms_lambda) - (energy_gs/num_atoms_gs))/S
        energies.write('{}  {} {} \n'.format(formula, diff_energy, diff_energy_S))

    energies.close()

if __name__ == "__main__":
    main()
