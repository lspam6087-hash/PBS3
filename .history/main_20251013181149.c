/******************************************************************************/ 
/*                                                                            */
/*  A Molecular Dynamics simulation of Lennard-Jones particles                */
/*                                                                            */
/*	This code is part of the course "Particle-based Simulations"              */
/*  taught at Eindhoven University of Technology.                             */
/*  No part of this code may be reproduced without permission of the author:  */
/*  Dr. Ir. E.A.J.F. Peters                                                   */
/*                                                                            */
/*  Dr. Ir. J.T. Padding:    version 1.1, 30/1/2013                           */
/*  Jeroen Hofman:           version 1.2, 28/7/2015                           */
/*  Dr. Ir. E.A.J.F. Peters: version 6.1, 17/9/2025                           */
/******************************************************************************/ 

/** 
 * For the 2025 PBS assignment, the code needs to be extended:
 * 
 * - Initialize vectors.type so particles get the proper type 
 * - Implement bonds in initialise_bonds in file initialise.c 
 * - Implement the bonded and non-bonded force in forces.c (Make forces type-dependent) 
 * - Change the particle position initialization to account for bond lengths, angles and dihedrals 
 * - Implement a Berendsen thermostat in dynamics.c 
 * - Implement the needed on-the-fly data analysis
 */ 

#include <stdio.h> 
#include <stdlib.h> 
#include <math.h> 
#include "constants.h" 
#include "structs.h" 
#include "setparameters.h" 
#include "initialise.h" 
#include "nbrlist.h" 
#include "forces.h" 
#include "dynamics.h" 
#include "memory.h" 
#include "fileoutput.h"
#include "vector_functions.h"
#include "grf.h" 
#include "histogram.h"
#include "density.h"
#include "parametercalc.h"

#define INCLUDE_FORCES // If B2 needs to be checked. This needs to be commented
#define HISTOGRAM
#define NUMPART_CALC

/** 
 * @brief Main MD simulation code. After initialization, 
 * a velocity-Verlet scheme is executed for a specified number of time steps. 
 * 
 * This simulation runs for a binary mixture of methane and ethane. 
 * 
 * Tasks implemented in this assignment:
 * 
 * 1. Add support for multiple particle types (CH3 and CH2) and handle non-bonded forces accordingly.
 * 2. Modify the code to include bonded interactions for n-butane molecules.
 * 3. Change the particle position initialization to account for bond lengths, angles and dihedrals.
 * 4. Implement the Berendsen thermostat to maintain the system temperature (NVT ensemble).
 * 5. Perform on-the-fly analysis for velocity distribution, torsion angle distribution and mean-square displacement.
 * 
 * @return int 0 on success, non-zero on failure. 
 */ 
int main(void) 
{ 
    struct Vectors vectors; 
    struct Parameters parameters; 
    struct Nbrlist nbrlist; 
    struct VelHist p_vhist;
    size_t step; 
    double Ekin, Epot, time, grcount; 

    // Step 1: Set the simulation parameters from input files
    set_parameters(&parameters, &p_vhist); 

    #ifdef NUMPART_CALC
        parameters.num_dt_steps = 100000;
        parameters.reset_chi_file = 0;
        parameters.delta_a = 12.0;
        parameters.L = (struct Vec3D){8.0, 8.0, 20.0}; // Set box dimensions for number of particles calculation
        num_part_calc(&parameters);
    #endif

    // Step 2: Allocate memory for particles, forces, and neighbor lists
    alloc_memory(&parameters, &vectors, &nbrlist, &p_vhist); 

    // Check if a restart is required
    if (parameters.load_restart == 1) 
    { 
        load_restart(&parameters, &vectors, &p_vhist); 
        initialise_structure(&parameters, &vectors, &nbrlist); 
        step = 0; 
        time = 0.0; 
    } 
    else 
    {   
    /// \todo Initialize particle types (CH3 and CH2) in vectors.type array
    /// \todo Implement the bonds between the UA of n-butane in initialise_bonds (initialise.c)

        initialise(&parameters, &vectors, &nbrlist, &step, &time); 
    }

    // Step 3: Build the neighbor list for non-bonded interactions
    build_nbrlist(&parameters, &vectors, &nbrlist); 

    // Step 4: Calculate initial forces (non-bonded and bonded if implemented)
    Epot = calculate_forces(&parameters, &nbrlist, &vectors); 

    // Output initial particle positions in PDB format
    record_trajectories_pdb(1, &parameters, &vectors, time); 

    initialise_grf(&parameters, &vectors);

    #ifdef HISTOGRAM
        initialize_hist(&parameters, &vectors, step, &p_vhist);
    #endif

    // Main MD loop using velocity-Verlet integration
    while (step < parameters.num_dt_steps) 
    { 
        step++; 
        time += parameters.dt; 

        // Update velocities (half-step)
        /// \todo Implement the use of type-dependent masses
        Ekin = update_velocities_half_dt(&parameters, &nbrlist, &vectors); 

        /// \todo Implement and apply the Berendsen thermostat to maintain temperature (dynamics.c)
        thermostat(&parameters, &vectors, Ekin); 

        // Update positions
        update_positions(&parameters, &nbrlist, &vectors); 

        // Apply boundary conditions
        boundary_conditions(&parameters, &vectors); 

        // Rebuild neighbor list if needed
        update_nbrlist(&parameters, &vectors, &nbrlist); 

        // Calculate forces for the current configuration (bonded forces if implemented)
        Epot = calculate_forces(&parameters, &nbrlist, &vectors); 

        // Final velocity update (half-step)
        Ekin = update_velocities_half_dt(&parameters, &nbrlist, &vectors); 

        // Update the GRF
        // grcount = update_grf(&parameters, &vectors);
        grcount = update_grf2(&parameters, &nbrlist, &vectors);

        // Update the Histogram
        update_hist(&parameters, &vectors, step, &p_vhist);

        // Output system state every 'num_dt_pdb' steps
        if (step % parameters.num_dt_pdb == 0) 
            record_trajectories_pdb(0, &parameters, &vectors, time); 

        // Save restart file every 'num_dt_restart' steps
        if (step % parameters.num_dt_restart == 0) 
            save_restart(&parameters, &vectors); 

        /// \todo Implement on-the-fly analysis of velocity distribution, torsion angle distribution and mean-square displacement
        #ifdef HISTOGRAM
            // Update the Histogram
            update_hist(&parameters, &vectors, step, &p_vhist);
            update_phi_hist(&parameters, &vectors, &p_vhist);
        #endif

        // Print to the screen to monitor the progress of the simulation
        /// \todo Write the output (also) to file, and extend the output with temperature
        printf("Step %lu, Time %.4f, Epot %f, Ekin %f, Etot %f\n", (long unsigned)step, time, Epot, Ekin, Epot + Ekin);

        // Save the relevant parameters for later data analysis
        record_diagnostics_csv((step == 1) ? 1 : 0, &parameters, time ,Ekin, Epot); 
    } 

    // 
    finalise_grf(&parameters, &vectors, grcount);

    record_histogram_csv(&parameters, &p_vhist, step);

    // Density profile
    initialize_density_histograms(&parameters, &vectors);
    accumulate_density_histogram(&parameters, &vectors);
    write_density_histograms(&parameters, &vectors);

    // Chi parameter calculation and recording
    initialize_phi_hist(&parameters, &vectors, &p_vhist);
    record_phi_histogram_csv(&parameters, &p_vhist);
    print_chi_csv(&parameters, &p_vhist);

    // Save final state
    save_restart(&parameters, &vectors); 

    // Step 5: Free memory and clean up
    free_memory(&vectors, &nbrlist, &p_vhist); 

    return 0; 
}
