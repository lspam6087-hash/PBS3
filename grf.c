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
    p_parameter->grcount =0.0;
    size_t nbin = p_parameter->nbin;
    // double dbin = (0.5 * L.x)/nbin;
    double *grbin = p_vectors->grbin;

    for (size_t i = 0; i <nbin; i++)
    {
        grbin[i] = 0.0;
    }
}

// Update
double update_grf(struct Parameters *p_parameter, struct Vectors *p_vectors){
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
    return grcount;
}

// Finalize
void finalise_grf(struct Parameters *p_parameter, struct Vectors *p_vectors, double grcount){
    size_t N = p_parameter->num_part;
    struct Vec3D L = p_parameter->L;
    size_t nbin = p_parameter->nbin;
    double dbin = p_parameter->dbin;
    double *grbin = p_vectors->grbin;
    double volume = 0.0;
    double rho = (N-1)/(L.x*L.y*L.z);


    for(size_t i=0; i<(nbin-1); i++){
        volume = 4.0/3.0 * PI * ((i+1)*(i+1)*(i+1) - (i)*(i)*(i)) * dbin*dbin*dbin;
        grbin[i] = grbin[i] / (rho * volume * grcount * N);

        if(i == 0){
            FILE* fp = fopen("data/grf_data.dat", "w");
            fprintf(fp, "%.8f %.12f\n", (i+0.5)*dbin, grbin[i]);
        } else {
            FILE* fp = fopen("data/grf_data.dat", "a");
            fprintf(fp, "%.8f %.12f\n", (i+0.5)*dbin, grbin[i]);
        }
    }
}


// This function calculates non-bonded forces between particles using the neighbor list.
// The potential energy and forces are calculated using the Lennard-Jones potential.
double update_grf2(struct Parameters *p_parameters, struct Nbrlist *p_nbrlist, struct Vectors *p_vectors)
{
    double r_cutsq;
    struct DeltaR rij;
    struct Pair *nbr = p_nbrlist->nbr;
    const size_t num_nbrs = p_nbrlist->num_nbrs;

    struct Vec3D *r = p_vectors->r;
    struct Vec3D *v = p_vectors->v;
    struct Vec3D L = p_parameters->L;

    r_cutsq = p_parameters->r_cut * p_parameters->r_cut;

    double Epot = 0.0;

    size_t nbin = p_parameters->nbin;
    double dbin = p_parameters->dbin;

    double *grbin = p_vectors->grbin;
    

    // Loop through the neighbor list and calculate the forces for each particle pair
    for (size_t k = 0; k < num_nbrs; k++)
    {
        rij = nbr[k].rij;
        size_t i = nbr[k].i;
        size_t j = nbr[k].j;

        double type_i = p_vectors->type[i];
        double type_j = p_vectors->type[j];

        // Compute forces if the distance is smaller than the cutoff distance
        if (sqrt(rij.sq) < sqrt(r_cutsq))
        {
            double rij_norm = sqrt(rij.sq);
            int ibin = (int)(rij_norm / dbin);
            if (ibin < nbin){
                grbin[ibin] += 2;
            }
            
        }          
    }
    p_parameters->grcount += 1;
    return p_parameters->grcount;  // Return the potential energy due to non-bonded interactions
}