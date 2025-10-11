#include "parametercalc.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>
#include "structs.h"

void num_part_calc(struct Parameters *p_parameters){
    size_t n_part = 0;

    double rho = p_parameters->rho;
    double N_A = 6.022e23;
    struct Vec3D L = p_parameters->L;

    double V = L.x * L.y * L.z;

    //Calculation of number of particles
    n_part = (size_t) (rho * V);

    p_parameters->num_part = n_part;
}

// This function calculates the volume fraction of type A particles
// Assumes that the volume of all of the beads (monomers) is the same
void volfrac_calc(struct Parameters *p_parameters, struct Vectors *p_vectors){
    double phi_A, N_A = 0.0;
    double num_part = p_parameters->num_part;

    for (size_t i = 0; i < num_part; i++){
        if (p_vectors->type[i] == 0){ // Assumes type A is represented by 0
            N_A++;
        }
    }

    phi_A = N_A / num_part;

    p_parameters->phi_A = phi_A;
    p_parameters->N_A = N_A;
}

void chi_calculation(struct Parameters *p_parameters, struct Vectors *p_vectors){
    // Calculate the volume fraction first
    volfrac_calc(p_parameters, p_vectors);

    double chi = 0.0;
    double phi_A = p_parameters->phi_A;
    double N_A = p_parameters->N_A;

    chi = (1/N_A) * (log((1-phi_A)/phi_A))/(1-2*phi_A);

    p_parameters->chi = chi;
}    

void record_chi_csv(struct Parameters *p_parameters){
    FILE *fp_chi;
    char filename[1024];
    snprintf(filename, 1024, "%s", p_parameters->filename_chi_data);
    size_t reset = p_parameters->reset_chi_file;

    double chi = p_parameters->chi;
    double rho = p_parameters->rho;
    double delta_a = p_parameters->delta_a;
    size_t N = p_parameters->amount_mon;

    if (reset == 1){
        // If reset is 1, we overwrite the file
        fp_chi = fopen(filename, "w");
        fprintf(fp_chi, "Chi, Rho, delta_a, N\n");
    } else {
        // If reset is 0, we append to the file
        fp_chi = fopen(filename, "a");
    }
    
    if (!fp_chi) {
        printf("Error opening histogram file %s!\n", filename);
        return;
    }

    // Write the chi value, density and delta_a to the file
    fprintf(fp_chi, "%f, %f, %f, %d\n", chi, rho, delta_a, N);

    fclose(fp_chi);
}
