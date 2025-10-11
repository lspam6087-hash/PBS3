#ifndef PARAMETERCALC_H_
#define PARAMETERCALC_H_
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>
#include "structs.h"

void num_part_calc(struct Parameters *p_parameters);

/**
 * @brief Calculates the volume fraction
 * 
 * @return volume fraction of type A particles and number of type A particles
 */
void volfrac_calc(struct Parameters *p_parameters, struct Vectors *p_vectors);

/**
 * @brief Calculates the Chi parameter wih the volume fraction of type A particles 
 * 
 * @return Chi parameter
 */
void chi_calculation(struct Parameters *p_parameters, struct Vectors *p_vectors);

void record_chi_csv(struct Parameters *p_parameters);
#endif 
