import pandas as pd
import numpy as np
import scanner 

#Charges of amino acids
aa_chages = {
    "R": 1.0,  "K": 1.0, "D": -1.0, "E": -1.0, 
    "H": 0.5, "A": 0.0, "Q": 0.0, "N": 0.0,
    "S": 0.0, "T": 0.0, "G": 0.0, "V": 0.0, 
    "L": 0.0, "I": 0.0, "F": 0.0, "Y": 0.0, 
    "W": 0.0, "M": 0.0, "C": 0.0, "P": 0.0
}

#Dictionary of mutations used
mutations = [
    ("WT", "WT", "WT", "WT"),
    (248, "ARG", "Q", "Pathogenic"), 
    (273, "ARG", "H", "Pathogenic"), 
    (280, "ARG", "K", "Pathogenic"), 
    (120, "LYS", "A", "Pathogenic"), 
    (241, "SER", "A", "Pathogenic"), 
    (175, "ARG", "H", "Pathogenic"), 
    (245, "GLY", "S", "Pathogenic"), 
    (176, "CYS", "F", "Pathogenic"), 
    (179, "HIS", "R", "Pathogenic"), 
    (220, "TYR", "C", "Pathogenic"), 
    (166, "SER", "A", "Benign"),
    (213, "ARG", "Q", "Pathogenic"),
    (282, "ARG", "W", "Pathogenic"),
    (249, "ARG", "S", "Pathogenic"),
    (128, "PRO", "T", "Benign"),
    (110, "ARG", "C", "Pathogenic"),
    (215, "SER", "R", "Pathogenic"),
    (191, "PRO", "A", "Benign"),
    (229, "CYS", "F", "Benign"),
    (261, "SER", "G", "Benign")
]

def run_experiment():
    
    pdb_file = "project/1TSR.pdb"
    p53_chains = ["A", "B", "C"]
    dna_chains = ["E", "F"]
    
    cA, qA, sA, eA = scanner.parse_pdb_to_arrays(pdb_file, p53_chains) #extracts the parameter
    cB, qB, sB, eB = scanner.parse_pdb_to_arrays(pdb_file, dna_chains)
    nA, nB = len(qA), len(qB)
    
    results_log = []

    for res_idx, original_aa, target_aa, type in mutations:
        if res_idx == "WT":
            current_qA = qA  # Copies the original charge array to be used
        else:
            print(f"Mutating {original_aa}{res_idx} -> {target_aa}...")
            
            current_qA = np.copy(qA)
            new_charge_val = aa_chages.get(target_aa, 0.0) #gets the charge for the mutated residue
            
            indices = []
            cursor = 0 #cursor to trace the pdb file and find the res index
            found_check = 0 # To check if the residue index is not in the pdb file
            with open(pdb_file, 'r') as f:
                for line in f:
                    if line.startswith("ATOM") and line[21] in p53_chains:
                        if int(line[22:26]) == res_idx:
                            indices.append(cursor) #appends the indices for all the atoms having the residue index as res_idx
                            found_check+=1 #increments the found_check
                        cursor += 1 #increments the cursor
            
            if found_check == 0: 
                print(f"{res_idx} not found in {pdb_file}")

            if indices:
                for i in indices: #zeroes the charge of all the atoms of res_idx
                    current_qA[i] = 0.0
                
                last_atom_idx = indices[-1] #assigns the overall charge of the mutated residue to the last atom of the res to simulate the change in residue charge
                current_qA[last_atom_idx] = new_charge_val

        results = np.zeros(3, dtype=np.float64)
        scanner.lib.calculate_interface_energy( #calculates the interface energy of the mutated protein
            nA, cA, current_qA, sA, eA, 
            nB, cB, qB, sB, eB, 
            results
        )
        
        results_log.append({
            "Mutation": f"{original_aa}{res_idx}{target_aa}" if res_idx != "WT" else "Wild Type",
            "Target_AA": target_aa,
            "Type": type,
            "Total_Energy": results[2],
            "Elec_Energy": results[1]
        })

    df = pd.DataFrame(results_log)
    
    wt_energy = df[df['Mutation'] == 'Wild Type']['Total_Energy'].values[0]
    df['Delta_Delta_E'] = df['Total_Energy'] - wt_energy #gives us $\Delta\Delta E$ for each mutation
    print(df)
    df.to_csv("scan_results.csv", index=False) #saves it in the file

if __name__ == "__main__":
    run_experiment()
