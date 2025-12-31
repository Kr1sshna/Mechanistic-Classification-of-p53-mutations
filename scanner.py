import ctypes
import numpy as np
import os
import forcefield
import pandas as pd

# Load C Engine
lib_path = os.path.join(os.getcwd(), 'project/libengine.so')
lib = ctypes.CDLL(lib_path)

#define pointer to double datatype
ptr_double = np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS') 

#define the arguments datatypes of the c engine function
lib.calculate_interface_energy.argtypes = [
    ctypes.c_int, ptr_double, ptr_double, ptr_double, ptr_double, # Protein A arrays
    ctypes.c_int, ptr_double, ptr_double, ptr_double, ptr_double, # Protein B arrays
    ptr_double # Results array
]

#parses the pdb file and extract charge, sigma and epsilon value of every atom
def parse_pdb_to_arrays(filename, target_chain):
    coords = [] #initiates empty arrays for our parameters
    charges = []
    sigmas = []
    epsilons = []
    
    with open(filename, 'r') as f:
        for line in f:
            if line.startswith("ATOM") or line.startswith("HETATM"): #reads only the line which contains our required value
                chain_id = line[21]
                if chain_id not in target_chain:
                    continue
                
                x = float(line[30:38]) #extracts the x, y and z coordinate
                y = float(line[38:46])
                z = float(line[46:54])
                
                res_name = line[17:20] #extracts the res_name, atom_name and element
                atom_name = line[12:16]
                element = line[76:78]
                
                q, s, e = forcefield.get_atom_parameters(res_name, atom_name, element) #calls the function get_atom_parameters declared in forcefield.py

                coords.extend([x, y, z])  #appends all the parameters in their respective arrays
                charges.append(q)
                sigmas.append(s)
                epsilons.append(e)
                
    return (np.array(coords, dtype=np.float64), #returns the array for each parameters
            np.array(charges, dtype=np.float64),
            np.array(sigmas, dtype=np.float64),
            np.array(epsilons, dtype=np.float64))

def main():
    
    pdb_file = "project/1TSR_Repair.pdb" #path of the pdb file
    print(f"Loading {pdb_file}...")
    
    chain_p53 = ["A","B","C"] #chain id of the trimer of p53 protein in the 1TSR.pdb
    chain_dna = ["E","F"] #chain id of the DNA sequence
    
    print(f"Parsing Chain {chain_p53} (p53)...\n")
    cA, qA, sA, eA = parse_pdb_to_arrays(pdb_file, chain_p53)
    nA = len(qA)
    df_p53 = pd.DataFrame({
        "Charges_p53": qA,
        "Sigmas_p53": sA,
        "Epsilons_p53": eA
    })

    print()
    print(f"Parsing Chain {chain_dna} (DNA)..].\n")
    cB, qB, sB, eB = parse_pdb_to_arrays(pdb_file, chain_dna) #extracts p53 and DNA atoms's parameter
    nB = len(qB)
    df_dna = pd.DataFrame({
        "Charges_DNA": qB,
        "Sigmas_DNA": sB,
        "Epsilons_DNA": eB
    })
    
    df_p53.to_csv("p53_data.csv") 
    df_dna.to_csv("dna_data.csv")

    print(f"Setup Complete: {nA} atoms vs {nB} atoms.")
    
    results = np.zeros(3, dtype=np.float64)
    
    # Calling the C engine function
    print("Running Physics Engine...")
    lib.calculate_interface_energy(nA, cA, qA, sA, eA, nB, cB, qB, sB, eB, results)
    
    print("\n--- RESULTS ---")
    print(f"Van der Waals (Shape):  {results[0]:.4f} kcal/mol")
    print(f"Electrostatics (Charge):{results[1]:.4f} kcal/mol")
    print(f"TOTAL BINDING ENERGY:   {results[2]:.4f} kcal/mol")

if __name__ == "__main__":
    main()
