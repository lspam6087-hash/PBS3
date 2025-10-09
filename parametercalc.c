#include "parametercalc.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>
#include "structs.h"

// This function calculates the volume fraction of type A particles
// Assumes that the volume of all of the beads (monomers) is the same
double volfrac_calc(struct Parameters *p_parameters, struct Vectors *p_vectors){
    double phi_A, N_A = 0.0;
    double num_part = p_parameters->num_part;

    for (size_t i = 0; i < num_part; i++){
        if (p_vectors->type[i] == 0){ // Assumes type A is represented by 0
            N_A++;
        }
    }

    phi_A = N_A / num_part;
    return phi_A, N_A;
}

double chi_calculation(struct Parameters *p_parameters, struct Vectors *p_vectors){
    double chi = 0.0;

    double phi_A, N_A = volfrac_calc(p_parameters, p_vectors);

    chi = (1/N_A) * (log((1-phi_A)/phi_A))/(1-2*phi_A);

    return chi;
}    
