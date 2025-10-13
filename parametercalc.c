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

// // This function calculates the volume fraction of type A particles
// // Assumes that the volume of all of the beads (monomers) is the same
// void volfrac_calc(struct Parameters *p_parameters, struct Vectors *p_vectors){
//     double phi_A, N_A = 0.0;
//     double num_part = p_parameters->num_part;
//     for (size_t i = 0; i < num_part; i++){
//         if (p_vectors->type[i] == 0){ // Assumes type A is represented by 0
//             N_A++;
//         }
//     }
//     phi_A = N_A / num_part;
//     p_parameters->phi_A = phi_A;
//     p_parameters->N_A = N_A;
// }

// void chi_calculation(struct Parameters *p_parameters, struct Vectors *p_vectors){
//     // Calculate the volume fraction first
//     volfrac_calc(p_parameters, p_vectors);
//     double chi = 0.0;
//     double phi_A = p_parameters->phi_A;
//     double N_A = p_parameters->N_A;
//     chi = (1/N_A) * (log((1-phi_A)/phi_A))/(1-2*phi_A);
//     p_parameters->chi = chi;
// }    

double average_phi_A_first_half(struct VelHist *p_vhist) {
    size_t half = p_vhist->nbins / 2;
    double sum = 0.0;
    size_t valid_bins = 0;
    for (size_t i = 0; i < half; i++) {
        if (p_vhist->counts[i] > 0) {
            sum += (double)p_vhist->typeA_counts[i] / p_vhist->counts[i];
            valid_bins++;
        }
    }
    return valid_bins > 0 ? sum / valid_bins : 0.0;
}

double chi_calculation(struct Parameters *p_parameters, struct VelHist *p_vhist){
    // Calculate the average volume fraction from the first half of bins
    double phi_A = average_phi_A_first_half(p_vhist);
    double N_A = p_parameters->amount_mon; // or sum typeA_counts in first half

    double chi = 0.0;
    chi = (1/N_A) * (log((1-phi_A)/phi_A))/(1-2*phi_A);

    p_parameters->chi = chi;
    p_parameters->phi_A = phi_A;

    return chi;
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

void initialize_phi_hist(struct Parameters *p_parameters, struct Vectors *p_vectors, struct VelHist *p_vhist){
    size_t nbins = p_vhist->nbins;
    double xmax = p_parameters->L.x;
    double xmin = 0.0;
    double dx = (xmax - xmin)/nbins;
    p_vhist->dv = dx; // rename to dx for clarity

    for (size_t i=0; i<nbins; i++)
    {
        p_vhist->bin_centers[i] = xmin + (i + 0.5) * dx;
        p_vhist->typeA_counts[i] = 0;
        p_vhist->counts[i] = 0;
    }
}

void update_phi_hist(struct Parameters *p_parameters, struct Vectors *p_vectors, struct VelHist *p_vhist){
    size_t nbins = p_vhist->nbins;
    size_t n = p_parameters->num_part;
    double dx = p_vhist->dv;

    for (size_t i = 0; i < n; i++){
        double x = p_vectors->r[i].x;
        int bin = (int)floor(x / dx);
        if (bin >= 0 && bin < nbins) {
            p_vhist->counts[bin]++;
            if (p_vectors->type[i] == 0) // type A
                p_vhist->typeA_counts[bin]++;
        }
    }
    p_vhist->total_counts++; // total samples
}

void record_phi_histogram_csv(struct Parameters *p_parameters,struct VelHist *p_vhist)
{
    FILE *fp_hist;
    char filename[1024];
    snprintf(filename, 1024, "%s", p_parameters->filename_hist);

    fp_hist = fopen(filename, "w");
    if (!fp_hist) {
        printf("Error opening histogram file %s!\n", filename);
        return;
    }

    //print header
    fprintf(fp_hist, "# samples=%zu bins=%d xmax=%g\n", p_vhist->total_counts, p_vhist->nbins, p_parameters->L.x);
    fprintf(fp_hist, "x_center, volfrac_A, count_A, count_total\n");

    for (int i = 0; i < p_vhist->nbins; i++) 
    {
        size_t countA = p_vhist->typeA_counts[i];
        size_t countTotal = p_vhist->counts[i];
        double center = p_vhist->bin_centers[i];
        double volfrac = countTotal > 0 ? ((double)countA / countTotal) : 0.0;
        fprintf(fp_hist, "%.6e,%.6e,%zu,%zu\n", center, volfrac, countA, countTotal);
    }
    fclose(fp_hist);
}


void print_chi_csv(struct Parameters *p_parameters, struct VelHist *p_vhist) {
    // Calculate average phi_A from first half of bins
    double phi_A = average_phi_A_first_half(p_vhist);
    double N_A = p_parameters->N_A;
    double chi = chi_calculation(p_parameters, p_vhist);

    double rho = p_parameters->rho;
    double delta_a = p_parameters->delta_a;
    size_t N = p_parameters->amount_mon;

    FILE *fp_chi;
    char filename[1024];
    snprintf(filename, 1024, "%s", p_parameters->filename_chi_data);

    if (p_parameters->reset_chi_file == 1){
        fp_chi = fopen(filename, "w");
        fprintf(fp_chi,"Chi, Rho, delta_a, N\n");
    } else {
        fp_chi = fopen(filename, "a");
    }

    // Print in CSV format
    fprintf(fp_chi,"%.6f, %.6f, %.6f, %zu\n", chi, rho, delta_a, N);
}