#ifndef PARAMETERCALC_H_
#define PARAMETERCALC_H_
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>
#include "structs.h"

void num_part_calc(struct Parameters *p_parameters);

/**
 * @brief Calculates the volume fraction of the left size of the x-direction
 * 
 * @return volume fraction of type A particles and number of type A particles
 */
double average_phi_A_first_half(struct VelHist *p_vhist);

/**
 * @brief Calculates the Chi parameter wih the volume fraction of type A particles 
 * 
 * @return Chi parameter
 */
double chi_calculation(struct Parameters *p_parameters, struct VelHist *p_vhist);

void record_chi_csv(struct Parameters *p_parameters);

void initialize_phi_hist(struct Parameters *p_parameters, struct Vectors *p_vectors, struct VelHist *p_vhist);

void update_phi_hist(struct Parameters *p_parameters, struct Vectors *p_vectors, struct VelHist *p_vhist);

void record_phi_histogram_csv(struct Parameters *p_parameters,struct VelHist *p_vhist);

void print_chi_csv(struct Parameters *p_parameters, struct VelHist *p_vhist);
#endif 
