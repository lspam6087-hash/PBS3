/* density.h */
#ifndef DENSITY_H
#define DENSITY_H

#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "structs.h"

void initialize_density_histograms(struct Parameters *p_parameters, struct Vectors *p_vectors);

void accumulate_density_histogram(struct Parameters *p_parameters, struct Vectors *p_vectors);

void write_density_histograms(struct Parameters *p_parameters, struct Vectors *p_vectors);

#endif