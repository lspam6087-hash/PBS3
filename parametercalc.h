#ifndef PARAMETERCALC_H_
#define PARAMETERCALC_H_
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>
#include "structs.h"
/**
 * @brief Calculates the volume 
 * 
 * @return volume fraction of type A particles and number of type A particles
 */
double volfrac_calc(struct Parameters *p_parameters, struct Vectors *p_vectors);

/**
 * @brief Calculates the Chi parameter wih the volume fraction of type A particles 
 * 
 * @return Chi parameter
 */
double chi_calculation(struct Parameters *p_parameters, struct Vectors *p_vectors);

#endif 
