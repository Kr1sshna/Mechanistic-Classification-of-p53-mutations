#include <math.h>
#include <stdio.h>

/**
 * PHYSICAL CONSTANTS
 * 
 * K_COULOMB: Coulomb's constant in kcal/mol·Å·e²
 *             Used to convert electrostatic interactions to energy units
 * 
 * DIELECTRIC: Relative dielectric constant of the medium (water ~80, protein ~8)
 *             Higher values reduce electrostatic interactions
 * 
 * CUTOFF_SQ: Squared distance cutoff for interaction calculations (in Å²)
 *             Pairs beyond this distance are ignored for computational efficiency
 *             CUTOFF = sqrt(200) ≈ 14.14 Å
 */
#define K_COULOMB 332.0 
#define DIELECTRIC 8.0 
#define CUTOFF_SQ 200.0


//lennard_jones - Calculates Van der Waals interaction energy

double lennard_jones(double r2, double sigma, double epsilon) {
    // Calculates r^6 and r^12 from r²
    double r6 = r2 * r2 * r2;
    double r12 = r6 * r6;
    
    // Calculates sigma^6, sigma^12
    double s2 = sigma * sigma;
    double s6 = s2 * s2 * s2;
    double s12 = s6 * s6;

    // Apply LJ potential formula
    return 4.0 * epsilon * ((s12 / r12) - (s6 / r6));
}


//coulomb - Calculates electrostatic interaction energy
 
double coulomb(double r, double q1, double q2) {
    // Energy = K_COULOMB * q1 * q2 / (DIELECTRIC * r)
    return (K_COULOMB * q1 * q2) / (DIELECTRIC * r);
}

/**
 * calculate_interface_energy - Computes total interaction energy between two molecules
 * 
 * This function performs a pairwise calculation between all atoms in molecule A
 * and all atoms in molecule B. For each atom pair
 
 * Arguments:
 *   n_A & n_B             - Number of atoms in molecule A and molecule B
 *   coords_A & coords_B   - 1D array of x,y,z coordinates for molecule A and molecule B
 *   charges_A & charges_B - Partial charges for each atom in molecule A and molecule B
 *   sig_A & sig_B         - Lennard-Jones sigma parameters for molecule A and molecule B
 *   eps_A & eps_B         - Lennard-Jones epsilon parameters for molecule A and molecule B
 *   results               - Output array to store computed energies:
 *                           results[0] = total van der Waals energy
 *                           results[1] = total electrostatic energy
 *                           results[2] = total interface energy (VDW + electrostatic)
 */
void calculate_interface_energy(
    int n_A, double *coords_A, double *charges_A, double *sig_A, double *eps_A,
    int n_B, double *coords_B, double *charges_B, double *sig_B, double *eps_B,
    double *results
) {
    // Accumulators for total energies (initialized to zero)
    double total_vdw = 0.0;    // Van der Waals (Lennard-Jones) energy
    double total_elec = 0.0;   // Electrostatic (Coulomb) energy
    
    // Iterate through all atoms in molecule A
    for (int i = 0; i < n_A; i++) {
        // x, y, z coordinates for atom i of molecule A
        double ax = coords_A[3*i + 0];  // x-coordinate of atom i
        double ay = coords_A[3*i + 1];  // y-coordinate of atom i
        double az = coords_A[3*i + 2];  // z-coordinate of atom i
        
        // Get force field parameters for atom i
        double qa = charges_A[i];  // Partial charge of atom i
        double sa = sig_A[i];      // sigma of atom i
        double ea = eps_A[i];      // epsilon of atom i

        // Iterate through all atoms in molecule B
        for (int j = 0; j < n_B; j++) {
            // x, y, z coordinates for atom j of molecule B
            double bx = coords_B[3*j + 0];  // x-coordinate of atom j
            double by = coords_B[3*j + 1];  // y-coordinate of atom j
            double bz = coords_B[3*j + 2];  // z-coordinate of atom j
            
            // Calculate distance vector components between atoms i and j
            double dx = ax - bx;  // Difference in x-coordinates
            double dy = ay - by;  // Difference in y-coordinates
            double dz = az - bz;  // Difference in z-coordinates
            
            // Calculate squared distance (avoids expensive sqrt for cutoff check)
            double r2 = dx*dx + dy*dy + dz*dz;
            
            // Only compute energies if atoms are within cutoff and not overlapping
            // r2 < CUTOFF_SQ: Ignore distant atom pairs for efficiency
            // r2 > 0.1: Prevent numerical instability from overlapping atoms and steric clashes
            if (r2 < CUTOFF_SQ && r2 > 0.1) { 
                // Calculate the distance for energy calculations
                double r = sqrt(r2);
                
                // Get force field parameters for atom j
                double qb = charges_B[j];  // Partial charge of atom j
                double sb = sig_B[j];      // sigma of atom j
                double eb = eps_B[j];      // epsilon of atom j
                
                // Lorentz-Berthelot mixing rules:
                // Combine parameters from atoms i and j to get pairwise parameters
                double sig_mix = (sa + sb) * 0.5;
                double eps_mix = sqrt(ea * eb);

                // Calculate interaction energies and sum them up every iteration
                total_vdw += lennard_jones(r2, sig_mix, eps_mix);  // VDW contribution
                total_elec += coulomb(r, qa, qb);                   // Electrostatic contribution
            }
        }
    }
    
    // Store results in output array
    results[0] = total_vdw;              // Total Van der Waals energy
    results[1] = total_elec;             // Total Electrostatic energy
    results[2] = total_vdw + total_elec; // Final interface energy
}
