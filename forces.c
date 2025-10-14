#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "forces.h"
#include "constants.h"
#include "structs.h"
#include "nbrlist.h"
#include "vector_functions.h"
#include "random.h"

#define INCLUDE_CONSERVATIVE
#define INCLUDE_DISSIPATIVE
#define INCLUDE_RANDOM

// This function calculates all forces acting on the particles (bonded and non-bonded).
// It initializes the forces array, then calculates bond-stretch, angle-bend, dihedral-torsion,
// and non-bonded forces. The total potential energy is returned.
double calculate_forces(struct Parameters *p_parameters, struct Nbrlist *p_nbrlist, struct Vectors *p_vectors)
{
    struct Vec3D *f = p_vectors->f;
    size_t num_part = p_parameters->num_part;

    // Initialize the forces to zero for all particles
    for (size_t i = 0; i < num_part; i++)
        f[i] = (struct Vec3D){0.0, 0.0, 0.0};

    // Calculate the forces and accumulate the potential energy from each type of interaction
    double Epot = 0.0;
    Epot = calculate_forces_bond(p_parameters, p_vectors);
    Epot += calculate_forces_nb(p_parameters, p_nbrlist, p_vectors);
    
    return Epot;
}

// This function calculates bond-stretch forces based on the current positions of the bonded particles.
// It applies the minimum image convention to calculate the distance between bonded pairs and then
// computes the force and potential energy due to the bond interaction.
double calculate_forces_bond(struct Parameters *p_parameters, struct Vectors *p_vectors)
{
    double Epot = 0.0;
    struct Bond *bonds = p_vectors->bonds;
    size_t num_bonds = p_vectors->num_bonds;
    struct Vec3D *f = p_vectors->f;
    struct Vec3D *r = p_vectors->r;
    struct Vec3D L = p_parameters->L;
    struct Vec3D rij;
    struct Vec3D fi = {0.0, 0.0, 0.0};

    // Loop through each bond and calculate the forces
    for (size_t q = 0; q < num_bonds; ++q)
    {
        size_t i = bonds[q].i;
        size_t j = bonds[q].j;

        // Apply the minimum image convention for calculating distances
        rij.x = r[i].x - r[j].x;
        rij.x = rij.x - L.x * floor(rij.x / L.x + 0.5);
        rij.y = r[i].y - r[j].y;
        rij.y = rij.y - L.y * floor(rij.y / L.y + 0.5);
        rij.z = r[i].z - r[j].z;
        rij.z = rij.z - L.z * floor(rij.z / L.z + 0.5);

        /// \todo Provide the bond force calculation and assign forces to particles i and j

        // Determining the bead-spring force
        double r_cut = p_parameters->r_cut;
        double KbT = 1;                     // \todo Change this parameter if necessary
        double C = 2 * KbT * pow(r_cut, -2);

        fi.x = -C * rij.x;
        fi.y = -C * rij.y;
        fi.z = -C * rij.z;

        f[i].x += fi.x;
        f[i].y += fi.y;
        f[i].z += fi.z;
        f[j].x -= fi.x;
        f[j].y -= fi.y;
        f[j].z -= fi.z;
    }

    return Epot;  // Return the potential energy due to bond-stretch interactions
}

// This function calculates non-bonded forces between particles using the neighbor list.
// The potential energy and forces are calculated using the Lennard-Jones potential.
double calculate_forces_nb(struct Parameters *p_parameters, struct Nbrlist *p_nbrlist, struct Vectors *p_vectors)
{
    struct Vec3D df_c, df_d, df_r;
    double r_cutsq, fr;
    struct DeltaR rij;
    struct Pair *nbr = p_nbrlist->nbr;
    const size_t num_nbrs = p_nbrlist->num_nbrs;
    struct Vec3D *f = p_vectors->f;
    size_t num_part = p_parameters->num_part;

    struct Vec3D *r = p_vectors->r;
    struct Vec3D *v = p_vectors->v;
    struct Vec3D L = p_parameters->L;

    double a_ij = p_parameters->a_ij;
    double T = p_parameters->T;
    double delta_a = p_parameters->delta_a;

    r_cutsq = p_parameters->r_cut * p_parameters->r_cut;

    double Epot = 0.0;

    double dt = p_parameters->dt;
    double dt_sq = pow(dt, -0.5);
    

    // Loop through the neighbor list and calculate the forces for each particle pair
    for (size_t k = 0; k < num_nbrs; k++)
    {
        rij = nbr[k].rij;
        size_t i = nbr[k].i;
        size_t j = nbr[k].j;

        double type_i = p_vectors->type[i];
        double type_j = p_vectors->type[j];
        
        double rij_norm = sqrt(rij.x*rij.x + rij.y*rij.y + rij.z*rij.z);
        struct Vec3D rij_hat = {0.0, 0.0, 0.0};
        rij_hat.x = (1.0/rij_norm) * rij.x;
        rij_hat.y = (1.0/rij_norm) * rij.y;
        rij_hat.z = (1.0/rij_norm) * rij.z;
        
        if (type_i == type_j){
            a_ij = 75.0/3.0; 
        }
        else {
            a_ij = 75.0/3.0 + delta_a;
        }

        // Compute forces if the distance is smaller than the cutoff distance
        if (rij.sq < r_cutsq)
        {
            #ifdef INCLUDE_CONSERVATIVE
            /// \todo Make the LJ parameters type-dependent (CH3 and CH2)

            df_c.x = a_ij * (1 - rij_norm) * rij_hat.x;
            df_c.y = a_ij * (1 - rij_norm) * rij_hat.y;
            df_c.z = a_ij * (1 - rij_norm) * rij_hat.z;

            // Compute the force and apply it to both particles
            f[i].x += df_c.x;
            f[i].y += df_c.y;
            f[i].z += df_c.z;
            f[j].x -= df_c.x;
            f[j].y -= df_c.y;
            f[j].z -= df_c.z; 
            #endif

            // Determine the potential energy of the conservative forces
            Epot += a_ij / 2 * (1-rij_norm)*(1-rij_norm);

            #ifdef INCLUDE_DISSIPATIVE
            // Dissipative forces
            struct Vec3D vij = {0.0, 0.0, 0.0};
            double gamma = 4.5; // From the assignment
            double w_D = pow(1-rij_norm, 2); //Weight function

            vij.x = v[i].x - v[j].x;
            vij.y = v[i].y - v[j].y;
            vij.z = v[i].z - v[j].z;

            df_d.x = -gamma * w_D * dot(rij_hat, vij) * rij_hat.x;
            df_d.y = -gamma * w_D * dot(rij_hat, vij) * rij_hat.y;
            df_d.z = -gamma * w_D * dot(rij_hat, vij) * rij_hat.z;

            // Compute the force and apply it to both particles

            f[i].x += df_d.x;
            f[i].y += df_d.y;
            f[i].z += df_d.z;
            f[j].x -= df_d.x;
            f[j].y -= df_d.y;
            f[j].z -= df_d.z;
            #endif

            #ifdef INCLUDE_RANDOM
            // Random forces
            double sigma = 3; // From the assignment
            double w_R = pow(w_D, 0.5); //Weight function
            double theta_ij = gauss();
            // double theta_ji = dist_rand_uniform();
            
            df_r.x = sigma * w_R * theta_ij * dt_sq * rij_hat.x;
            df_r.y = sigma * w_R * theta_ij * dt_sq * rij_hat.y;
            df_r.z = sigma * w_R * theta_ij * dt_sq * rij_hat.z;

            f[i].x += df_r.x;
            f[i].y += df_r.y;
            f[i].z += df_r.z;
            f[j].x -= df_r.x;
            f[j].y -= df_r.y;
            f[j].z -= df_r.z;
            #endif
        }          
    }
    return Epot;  // Return the potential energy due to non-bonded interactions
}
