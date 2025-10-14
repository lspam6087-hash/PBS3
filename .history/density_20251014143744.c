#include "density.h"
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "structs.h"

void initialize_density_histograms(struct Parameters *p_parameters, struct Vectors *p_vectors){
    for (int i = 0; i < p_parameters->nbins_dens; i++) {
        p_vectors->density_A[i] = 0.0;
        p_vectors->density_B[i] = 0.0;
        p_vectors->density_total[i] = 0.0;
    }
    p_vectors->density_samples = 0;
}


void accumulate_density_histogram(struct Parameters *p_parameters, struct Vectors *p_vectors){
                     
    size_t nbins = p_parameters->nbins_dens;
    size_t num_particles = p_parameters->num_part;
    struct Vec3D L = p_parameters->L;
    double bin_width = L.x / nbins;
    double bin_volume = bin_width * L.y * L.z;

    double *density_A = p_vectors->density_A;
    double *density_B = p_vectors->density_B;
    double *density_total = p_vectors->density_total;

    struct Vec3D *r = p_vectors->r;
    int *type = p_vectors->type;

    // Looping over all particles to accumalate densities
    for (int i = 0; i < num_particles; i++) {
        double x = r[i].x;
        int bin = (int)floor(r[i].x / bin_width);

        // Assigning densities to respective bins
        if (bin >= 0 && bin < nbins) {
            if (type[i] == 0){
                density_A[bin] += 1.0;
            }
            else if (type[i] == 1){
                density_B[bin] += 1.0;
            }
            density_total[bin] += 1.0;
        }
    }

    // Normalize densities by bin volume
    // for (int i = 0; i < nbins; i++) {
    //     density_A[i] /= bin_volume;
    //     density_B[i] /= bin_volume;
    //     density_total[i] /= bin_volume;
    // }

    // double total_N = 0.0;
    // for (int i = 0; i < nbins; i++) {
    //     total_N += density_total[i] * bin_volume;
    // }
    // // printing the total_N and num_part to check if they are equal
    // printf("Total N from density histogram: %.f\nNumber of particles: %lu\n", total_N, (long unsigned)num_particles);
}

void write_density_histograms(struct Parameters *p_parameters, struct Vectors *p_vectors) {
    size_t nbins = p_parameters->nbins_dens;

    double *density_A = p_vectors->density_A;
    double *density_B = p_vectors->density_B;
    double *density_total = p_vectors->density_total;

    struct Vec3D L = p_parameters->L;
    double bin_width = L.x / nbins;

    FILE *fp = fopen(p_parameters->filename_hist_dens, "w");
    if (!fp) {
        perror("Error opening density histogram file");
        return;
    }

    fprintf(fp, "x_center, density_A, density_B, density_total\n");

    for (int i = 0; i < nbins; i++) {
        double x_center = (i + 0.5) * bin_width;
        fprintf(fp, "%f, %f, %f, %f\n", x_center, density_A[i], density_B[i], density_total[i]);
    }

    fclose(fp);
}