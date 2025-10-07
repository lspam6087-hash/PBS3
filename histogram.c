/* velocity_histogram.c */
#include "histogram.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>

// \todo make this an on-the-fly histogram
// Setting the parameters
struct VelHist *vh_create(size_t nbins, double vmin, double vmax) {
    if (nbins == 0 || vmax <= vmin) return NULL;
    struct VelHist *h = (struct VelHist*) malloc(sizeof(struct VelHist));
    if (!h) return NULL;
    h->nbins = nbins;
    h->vmin = vmin;
    h->vmax = vmax;
    h->bin_width = (vmax - vmin) / (double)nbins;
    h->counts = (double*) calloc(nbins, sizeof(double));
    h->total_counts = 0;
    if (!h->counts) { free(h); return NULL; }
    return h;
}

void vh_add(struct VelHist *h, double value) {
    if (value < h->vmin || value >= h->vmax) return; // ignore out of range
    size_t bin = (size_t)((value - h->vmin) / h->bin_width);
    if (bin >= h->nbins) bin = h->nbins - 1; // clamp just in case
    h->counts[bin]++;
    h->total_counts++;
}

void write_hist(struct Parameters *p_parameters, struct Vectors *p_vectors, size_t step, struct VelHist *p_vhist){

    size_t nbins = p_parameters->nbins;
    double vmin = 0.0;
    double vmax = p_parameters->hist_vmax; 
    int N = p_parameters->num_part;
    struct Vec3D *v = p_vectors->v;
    size_t sample_interval = p_parameters->sample_interval; 

    // if (!p_vhist) {
    //     p_vhist = vh_create(nbins, vmin, vmax);
    // }

    if (step % sample_interval != 0) return;

    // Accumulate into histogram instead of writing speeds directly
    for (size_t i = 0; i < N; ++i) {
        double speed = sqrt(v[i].x*v[i].x + v[i].y*v[i].y + v[i].z*v[i].z);
        vh_add(p_vhist, speed);
    }

    // Example: every 100 intervals, dump the histogram
    if (step % sample_interval == 0) {
        if (step == 0){
            FILE *fp = fopen("data/hist.dat", "w");
            for (size_t b = 0; b < p_vhist->nbins; ++b) {
                double bin_center = p_vhist->vmin + (b + 0.5) * p_vhist->bin_width;
                fprintf(fp, "%f %f\n", bin_center, p_vhist->counts[b]);
                // fclose(fp);
            }
        } else {
            FILE *fp = fopen("data/hist.dat", "a");
            for (size_t b = 0; b < p_vhist->nbins; ++b) {
                double bin_center = p_vhist->vmin + (b + 0.5) * p_vhist->bin_width;
                fprintf(fp, "%.6f %lu\n", bin_center, p_vhist->counts[b]);
                // fclose(fp);
            }
        }
    }
}
