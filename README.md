# DFT Based Predictions for Curie Temperatures

This repository contains all the scripts necessary for setting up and extracting results from **VASP** calculations to predict the Curie temperature of magnetic materials using the DFT-based method developed in [Marian Arale Brännvall, Gabriel Persson, Luis Casillas-Trujillo, Rickard Armiento, and Björn Alling, *Predicting the Curie temperature of magnetic materials with automated calculations across chemistries and structures*, Phys. Rev. Materials 8, 114417, 2024.](https://journals.aps.org/prmaterials/abstract/10.1103/PhysRevMaterials.8.114417)

## Scripts Overview

### **Setup Scripts**
These scripts are used to prepare the input configurations needed for the VASP calculations:

1. **`get_random_magnetic_configuration.py`**  
   - **Description**: Generates a random noncollinear magnetic configuration for a crystal.  
   - **Input**:  
     - Number of atoms in the system.
     - Example: python get_random_magnetic_configuration.py 24 ; will create a MAGMOM line for 24 atoms.
   - **Output**:  
     - A noncollinear magnetic moment configuration where all moments have a size of 3 Bohr magnetons.  
   - **Use Case**: This configuration can be used as input in the `INCAR` file for a ground state search, set the MAGMOM-flag to this configuration.

2. **`get_best_DLM_config.py`**  
   - **Description**: Generates a disordered local moment (DLM) configuration with minimized short-range order (SRO).  
   - **Input**:  
     - Number of each species in the crystal.  
     - Initial magnetic moment sizes for each species.
     - Example: python get_best_DLM_config.py 48 16 0.5 2.2 ;This will produce MAGMOM and M_CONSTR lines where the first 48 atoms have unconstrained moments.
   - **Functionality**:  
     - Determines whether the moments for each species are constrained based on their input sizes (less than 0.75 Bohr magnetons = unconstrained).
     - There is a line, NNList = calc_neighbors(structure, 'POSCAR', tol=0.1), which determines between which neighboring atoms to calculate the SRO, the tol determines the radius within which the neighbors are found. This might need to be adjusted, check with get_NN.py where a similr line is found.   
     - Requires access to:  
       - `get_SRO.py`: Calculates the short-range order.  
       - `POSCAR`: Crystal structure information.  
   - **Output**:  
     - A noncollinear configuration with the lowest SRO among 50 random configurations.  
     - Configuration for the `M_CONSTR` flag in `INCAR`, which specifies which magnetic moments are constrained.  

### **Extraction Scripts**
These scripts are used to analyze the output from VASP calculations and gather the necessary data for Curie temperature predictions:

1. **`gather_energy_diff.py`**  
   - **Description**: Extracts the energy of the ground state and the DLM state respectively from VASP output files and calculates the energy differences. It also extracts the magnetic moment sizes from the DLM configuration and calculates the magnetic entropy. This is then used to calculate the energy difference divided by the magnetic entropy. The energy difference, the magentic entropy, and the energy difference divided by the magnetic entropy are printed out into a files called energies.
   - **Input**: python gather_energy_diff.py OSZICAR.GS OSZICAR.DLM OUTCAR.GS OUTCAR.DLM POSCAR.GS POSCAR.DLM ; Set the filenames appropriately so that the GS files are the names of the files from the ground-state search and the DLM-files from the final DLM-run. The files must be in the same place as this script.
   - **Dependencies**:
     - Requires `gather_mag_final.sh`. Make sure the script is executable!

2. **`get_NN.py`**  
   - **Description**: Identifies the number of nearest neighbors (NN) in the crystal structure. Prints it into a file called Number_of_nearest_neighbors.
   - **Input**: python get_NN.py POSCAR INCAR ; Set the file names to he POSCAR and INCAR corresponding to the DLM run. The files must be in the same place as this script. 
   - **Dependencies**:  
     - Requires `get_SRO.py` for short-range order calculations.
     - Requires `gather_mag_final.sh`. Make sure the script is executable!
       

---


## How to Use

1. **Setup Calculations**:
   - Use `get_random_magnetic_configuration.py` to create a starting configuration for ground state calculations.
   - Use `get_best_DLM_config.py` to prepare a DLM configuration with minimized SRO.

2. **Run VASP Calculations**:
   - Use the generated configurations in the `INCAR` and the relevant `POSCAR` files for VASP calculations.

3. **Extract Data**:
   - Use `gather_energy_diff.py` and `get_NN.py` to extract the required outputs from VASP calculations.

4. **Predict Curie Temperature**:
   - Use the extracted data in conjunction with the DFT-based method to predict the Curie temperature.

---

## Dependencies

- **Python**: Requires a version **greater than 3.8 but less than 3.12**  
- **Required Packages**: Listed in [`requirements.txt`](requirements.txt)  

---

## Notes

- Ensure all scripts are in the same directory for dependency.
- Input files (e.g., `POSCAR`, `INCAR`) must follow the VASP standard format.
