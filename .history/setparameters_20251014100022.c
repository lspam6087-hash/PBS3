#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "constants.h"
#include "structs.h"

// double box_size(struct Parameters *p_parameters)
// {
//   size_t n_part = p_parameters->num_part;
//   double rho = p_parameters->rho;
//   double N_A = 6.022e23;
//   double total_mass = (2*p_parameters->mass[0] + 2*p_parameters->mass[1])*1.0e-3;
//   double r_cut = p_parameters->r_cut;
//   double L, L3;

//   //Calculation of box size
//   L3 = total_mass * n_part/(rho * N_A);
//   L = pow(L3, 1.0/3.0)/(1.0e-10); //Conversion to [A]

//   if (L < 2*r_cut)
//     L = 2*r_cut;

//   return L;
// }

// Set the parameters of this simulation
void set_parameters(struct Parameters *p_parameters, struct VelHist *p_vhist)
{
// The parameters first 5 parameters are only used for demonstration puprposes
  p_parameters->kT = 1.0;                                   //thermal energy
  p_parameters->mass = 1.0;                                 //mass of a particle
  p_parameters->epsilon = 1.0;                              //LJ interaction strength
  p_parameters->sigma = 1.0;                                //LJ particle diameter
  p_parameters->a_ij = 0.0;                                 //Repulsion parameter
  p_parameters->T = 298.0;                                  //[K] Temperature
  p_parameters->delta_a = 12.0;                              //Excess repulsion
  p_parameters->rho = 3;                                    //Number density

// The parameters below control core functionalities of the code, but many values will need to be changed
  p_parameters->num_part = 4;                            //number of particles
  p_parameters->num_dt_steps = 25000;                        //number of time steps
  p_parameters->sample_interval = 10;                       // Sample interval
  p_parameters->exclude_12_nb = 0;                          // 1-2 connected atoms exluded from non-bonded interactions 
  p_parameters->exclude_13_nb = 0;                          // 1-3 connected atoms exluded from non-bonded interactions    
  p_parameters->dt = 0.01;                                  //integration time step
  p_parameters->r_cut = 1.0;                                //cut-off distance used for neigbor list
  // double Lbox = box_size(p_parameters);
  p_parameters->L = (struct Vec3D){2.1, 2.1, 2.1};          //Box size
  p_parameters->r_shell = 0.4;                              //shell thickness for neighbor list
  p_parameters->num_dt_pdb = 50;                           //number of time steps in between pdb outputs
  strcpy(p_parameters->filename_pdb, "trajectories");       //filename (without extension) for pdb file
  strcpy(p_parameters->filename_hist, "data/vel_histogramB3.csv");//filename histogram velocity
  strcpy(p_parameters->filename_hist_dens, "data/dens_histogramC2.csv");//filename histogram density
  p_parameters->rescale_output = 1;                         //factor used to rescale output lengthscale (Most visualisation programs identify bonds based on distances of order 1)
  p_parameters->load_restart = 0;                           //if equal 1 restart file is loaded
  strcpy(p_parameters->restart_in_filename, "restart.dat"); //filename for loaded restart file
  p_parameters->num_dt_restart = 1000;                      // number of time steps between saves
  strcpy(p_parameters->restart_out_filename, "restart.dat");//filename for saved restart file

  p_parameters->nbin = 100;//(int)sqrt(p_parameters->num_part);                                //Number of bins
  p_parameters->grcount = 0;                                //Counter of calls to update_gr 
  p_parameters->dbin = p_parameters->r_cut/p_parameters->nbin; //Bin width

  // These parameters are used for density profile calculation and plotting
  p_parameters->nbins_dens = 250;                            // Number of bins for density histogram
  p_parameters->amount_mon = 2;                             // Amount of monomers in a polymer
  p_parameters->binary_mix = 1;                             // 1 for binary mixture, 0 for a single component system

  // Parameters for chi calculation and recording
  p_parameters->chi = 0.0;
  p_parameters->phi_A = 0.0;
  p_parameters->N_A = 0.0;
  p_parameters->reset_chi_file = 0;
  strcpy(p_parameters->filename_chi_data, "data/chi_dataC4.csv");//filename histogram density




  // Parameters for the velocity histogram
  p_vhist->nbins = 1000;                                      //Number of bins
  p_vhist->vmin = 0.0;                                      //Min velocity
  p_vhist->vmax = 30.0;                                     //Max velocity
  p_vhist->bin_width = (p_vhist->vmax-p_vhist->vmin)/p_vhist->nbins;
  p_vhist->total_counts = 0.0;

  if (p_parameters->r_cut > p_parameters->L.x / 2.0)
    fprintf(stderr, "Warning! r_cut > Lx/2");
  if (p_parameters->r_cut > p_parameters->L.y / 2.0)
    fprintf(stderr, "Warning! r_cut > Ly/2");
  if (p_parameters->r_cut > p_parameters->L.z / 2.0)
    fprintf(stderr, "Warning! r_cut > Lz/2");
}
