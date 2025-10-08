/* velocity_histogram.h */
#ifndef VELOCITY_HISTOGRAM_H
#define VELOCITY_HISTOGRAM_H

#include <stddef.h>
#include "structs.h"


void initialize_hist(struct Parameters *p_parameters, struct Vectors *p_vectors, size_t step, struct VelHist *p_vhist);

void update_hist(struct Parameters *p_parameters, struct Vectors *p_vectors, size_t step, struct VelHist *p_vhist);

void record_histogram_csv(struct Parameters *p_parameters,struct VelHist *p_vhist, size_t step);

#endif