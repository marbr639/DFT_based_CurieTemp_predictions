# DFT Based Predicions for Curie Temperatures

This repository contains all the scripts necessary for setting up and extracting results from **VASP** calculations to predict the Curie temperature of magnetic materials using the **BCPAA-method** developed [Marian Arale Brännvall, Gabriel Persson, Luis Casillas-Trujillo, Rickard Armiento, and Björn Alling, *Predicting the Curie temperature of magnetic materials with automated calculations across chemistries and structures*, Phys. Rev. Materials 8, 114417, 2024.](https://journals.aps.org/prmaterials/abstract/10.1103/PhysRevMaterials.8.114417) 

## Scripts Overview

### **Setup Scripts**
These scripts are used to prepare the input configurations needed for the VASP calculations:

1. **`get_random_magnetic_configuration.py`**  
   - **Description**: Generates a random noncollinear magnetic configuration for a crystal.  
   - **Input**:  
     - Number of atoms in the system.  
   - **Output**:  
     - A noncollinear magnetic moment configuration where all moments have a size of 3 Bohr magnetons.  
   - **Use Case**: This configuration can be used as input in the `INCAR` file for a ground state search, set the MAGMOM-flag to this configuration.

2. **`get_best_DLM_config.py`**  
   - **Description**: Generates a disordered local moment (DLM) configuration with minimized short-range order (SRO).  
   - **Input**:  
     - Number of each species in the crystal.  
     - Initial magnetic moment sizes for each species.  
   - **Functionality**:  
     - Determines whether the moments for each species are constrained based on their input sizes (less than 0.75 Bohr magnetons = unconstrained).  
     - Requires access to:  
       - `get_SRO.py`: Calculates the short-range order.  
       - `POSCAR`: Crystal structure information.  
   - **Output**:  
     - A noncollinear configuration with the lowest SRO among 50 random configurations.  
     - Configuration for the `M_CONSTR` flag in `INCAR`, which specifies which magnetic moments are constrained.  

### **Extraction Scripts**
These scripts are used to analyze the output from VASP calculations and gather the necessary data for Curie temperature predictions:

1. **`gather_energy_diff.py`**  
   - **Description**: Extracts the energy of the groudn state and the DLM state respectivelyfrom VASP output files and calculates the energy differences. It also extracts the magnetic moment sizes from the DLM configuration and calculates the magnetic entropy. This is then used to calculate the energy difference divided by the magnetic entropy. The energy difference and the energy difference divided by the magnetic entropy are printed out into a files called energies.

2. **`get_NN.py`**  
   - **Description**: Identifies the number of nearest neighbors (NN) in the crystal structure. Prints it into a file called Number_of_nearest_neighbors.
   - **Dependencies**:  
     - Requires `get_SRO.py` for short-range order calculations.

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
   - Use the extracted data in conjunction with the BCPAA method to predict the Curie temperature.

---

## Dependencies

- Python 3.x
- VASP (Vienna Ab initio Simulation Package)

---

## Notes

- Ensure all scripts are in the same directory for dependency resolution (e.g., `get_SRO.py`).
- Input files (e.g., `POSCAR`, `INCAR`) must follow the VASP standard format.
