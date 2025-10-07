/* velocity_histogram.h */
#ifndef VELOCITY_HISTOGRAM_H
#define VELOCITY_HISTOGRAM_H

#include <stddef.h>
#include "structs.h"

/**
 * @brief Set the histogram struct used in the simulation.
 * 
 * @param nbins number of bins
 * @param vmin minimum speed
 * @param vmax maximum speed
 */
struct VelHist *vh_create(size_t nbins, double vmin, double vmax);
void vh_destroy(struct VelHist *h);
void vh_reset(struct VelHist *h);

void vh_write_ascii(struct VelHist *h, const char *filename); /* writes bin_center count */


void vh_add(struct VelHist *h, double speed);

/**
 * @brief Include the velocity of the particles in the histogram.
 * 
 * @param p_parameters used members: nbins, hist_vmax, num_part, L
 * @param p_vectors used members: v
 * @param step current step
 */
void write_hist(struct Parameters *p_parameters, struct Vectors *p_vectors, size_t step, struct VelHist *p_vhist);
#endif