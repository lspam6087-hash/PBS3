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

// Initialize
void initialise_grf(struct Parameters *p_parameter, struct Vectors *p_vectors){
    // struct Vec3D L = p_parameter->L;
    // size_t grcount = p_parameter->grcount;
    size_t nbin = p_parameter->nbin;
    // double dbin = (0.5 * L.x)/nbin;
    double *grbin = p_vectors->grbin;

    for (size_t i = 0; i <nbin; i++)
    {
        grbin[i] = 0.0;
    }
}

// Update
void update_grf(struct Parameters *p_parameter, struct Vectors *p_vectors){
    size_t grcount = p_parameter->grcount;
    struct Bond *bonds = p_vectors->bonds;
    size_t num_bonds = p_vectors->num_bonds;
    struct Vec3D *r = p_vectors->r;
    struct Vec3D L = p_parameter->L;
    struct Vec3D rij;
    size_t nbin = p_parameter->nbin;
    double dbin = p_parameter->dbin;

    double *grbin = p_vectors->grbin;

    // Loop through each bond and calculate the forces
    for (size_t q = 0; q < num_bonds; ++q){

        size_t i = bonds[q].i;
        size_t j = bonds[q].j;

        // Apply the minimum image convention for calculating distances
        rij.x = r[i].x - r[j].x;
        rij.x = rij.x - L.x * floor(rij.x / L.x + 0.5);
        rij.y = r[i].y - r[j].y;
        rij.y = rij.y - L.y * floor(rij.y / L.y + 0.5);
        rij.z = r[i].z - r[j].z;
        rij.z = rij.z - L.z * floor(rij.z / L.z + 0.5);

        double rij_norm = norm(rij);

        int ibin = (int)(rij_norm / dbin);
        if (ibin < nbin){
            grbin[ibin] += 2;
        }

        grcount += 1;
    }
}

// Finalize
void finalise_grf(struct Parameters *p_parameter, struct Vectors *p_vectors){
    size_t N = p_parameter->num_part;
    struct Vec3D L = p_parameter->L;
    size_t nbin = p_parameter->nbin;
    double dbin = p_parameter->dbin;
    double *grbin = p_vectors->grbin;
    size_t grcount = p_parameter->grcount;
    double volume = 0.0;
    double rho = (N-1)/(L.x*L.y*L.z);

    for(size_t i=0; i<(nbin-1); i++){
        volume = 4/3 * PI * ((i+1)*(i+1)*(i+1) - (i)*(i)*(i)) * dbin*dbin*dbin;
        grbin[i] = grbin[i] / (rho * volume * grcount * N);
        printf("%d, %f\n", (i+0.5)*dbin, grbin[i]);
    }
}