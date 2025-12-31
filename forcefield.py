# Charge parameters for atoms in the protein and DNA structure

ATOM_CHARGES = {
    # Positives
    "ARG": {"NH1": 1.0, "NH2": 1.0},
    "LYS": {"NZ": 1.0},
    "HIS": {"ND1": 0.8, "NE2": 0.2},
    
    # Negatives
    "ASP": {"OD1": -0.8, "OD2": -0.8, "CG": 0.6},
    "GLU": {"OE1": -0.8, "OE2": -0.8, "CD": 0.6},

    # DNA Residues
    "DA": {
        "O1P": -0.7, "O2P": -0.7, "OP1": -0.7, "OP2": -0.7, "P": 1.1,
        "O5'": 0.0, "C5'": 0.0, "C4'": 0.0, "O4'": 0.0, "C3'": 0.0, "O3'": 0.0, "C2'": 0.0, "C1'": 0.0,
        "N9": 0.0, "C8": 0.0, "N7": 0.0, "C5": 0.0, "C6": 0.0, "N6": 0.0, "N1": 0.0, "C2": 0.0, "N3": 0.0, "C4": 0.0
    },
    "DC": {
        "O1P": -0.7, "O2P": -0.7, "OP1": -0.7, "OP2": -0.7, "P": 1.1,
        "O5'": 0.0, "C5'": 0.0, "C4'": 0.0, "O4'": 0.0, "C3'": 0.0, "O3'": 0.0, "C2'": 0.0, "C1'": 0.0,
        "N1": 0.0, "C2": 0.0, "O2": 0.0, "N3": 0.0, "C4": 0.0, "N4": 0.0, "C5": 0.0, "C6": 0.0
    },
    "DG": {
        "O1P": -0.7, "O2P": -0.7, "OP1": -0.7, "OP2": -0.7, "P": 1.1,
        "O5'": 0.0, "C5'": 0.0, "C4'": 0.0, "O4'":0.0, "C3'":0.0, "O3'":0.0, "C2'":0.0, "C1'":0.0,
        "N9": 0.0, "C8": 0.0, "N7": 0.0, "C5": 0.0, "C6": 0.0, "O6": 0.0, "N1": 0.0, "C2": 0.0, "N2": 0.0, "N3": 0.0, "C4": 0.0
    },
    "DT": {
        "O1P": -0.7, "O2P": -0.7, "OP1": -0.7, "OP2": -0.7, "P": 1.1,
        "O5'": 0.0, "C5'": 0.0, "C4'": 0.0, "O4'": 0.0, "C3'": 0.0, "O3'": 0.0, "C2'": 0.0, "C1'": 0.0,
        "N1": 0.0, "C2": 0.0, "O2": 0.0, "N3": 0.0, "C4": 0.0, "O4": 0.0, "C5": 0.0, "C7": 0.0, "C6": 0.0
    },
}

# LJ sigma and epsilon value for each atom
ATOM_VDW = {
    "C": (3.4, 0.11),
    "N": (3.2, 0.17),
    "O": (2.9, 0.21),
    "S": (3.5, 0.25),
    "P": (4.2, 0.20),
    "H": (1.0, 0.02),
    "X": (3.0, 0.10)
}

#Function to extract partial charge, sigma and epsilon value of the atom
def get_atom_parameters(res_name, atom_name, element):
    q = 0.0             #initiates charge value to 0
    res = res_name.strip()  #strips the res_name from any space character
    atom = atom_name.strip() #strps the atom_name from any space character

    if res in ATOM_CHARGES:  #checks if the atom(of residue res) is in the dict ATOM_CHARGES, if satisfied, assigns charge to the value of the key or else 0
        if atom in ATOM_CHARGES[res].keys():
            q = ATOM_CHARGES[res][atom]
            
    elem = element.strip().upper()  #strips the element from any space character and converts it into uppercase
    sig = 0.0 #initiates sigma and epsilon value to 0
    eps = 0.0
    if elem in ATOM_VDW: #if elem is in the dict ATOM_VDW, assigns sigma and epsilon value to the elem key of the dict
        sig, eps = ATOM_VDW[elem]
    else: #else assign it to an hypothised value
        sig, eps = ATOM_VDW["X"]
            
    return q, sig, eps #return the paraeters