/* velocity_histogram.c */
#include "histogram.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>
#include "vector_functions.h"


void initialize_hist(struct Parameters *p_parameters, struct Vectors *p_vectors, size_t step, struct VelHist *p_vhist){
    size_t nbins = p_vhist->nbins;
    double vmax = p_vhist->vmax;
    double vmin = p_vhist->vmin;

    double dv = (vmax - vmin)/nbins;
    p_vhist->dv = dv;

    for (size_t i=0; i<nbins; i++)
    {
        p_vhist->bin_centers[i] = (i * 0.5) + dv;
        p_vhist->counts[i] = 0;
    }
}

void update_hist(struct Parameters *p_parameters, struct Vectors *p_vectors, size_t step, struct VelHist *p_vhist){
    size_t nbins = p_vhist->nbins;
    size_t n = p_parameters->num_part;
    double dv = p_vhist->dv;

    for (size_t i; i < n; i++){
        double vel = norm(p_vectors->v[i]);

        //Calculate the bin index and check if it is within bounds
        int bin = (int)floor(vel / dv);
        if (bin >= 0 && bin < nbins) 
            p_vhist->counts[bin]++;
    }
    p_vhist->total_counts++;
}

void record_histogram_csv(struct Parameters *p_parameters,struct VelHist *p_vhist, size_t step)
{
    FILE *fp_hist;
    char filename[1024];
    snprintf(filename, 1024, "%s", p_parameters->filename_hist);

    fp_hist = fopen(filename, "w");
    if (!fp_hist) {
        printf("Error opening histogram file %s!\n", filename);
        return;
    }

    //print header
    fprintf(fp_hist, "# step=%zu samples=%zu bins=%d vmax=%g\n", step, p_vhist->total_counts, p_vhist->nbins, p_vhist->vmax);
    fprintf(fp_hist, "v_center, density, count\n");

    //Calculate normalization
    double norm = 1.0 / (p_vhist->total_counts * p_parameters->num_part * p_vhist->dv);

    for (int i = 0; i < p_vhist->nbins; i++) 
    {
        size_t counts = p_vhist->counts[i];
        double centers = p_vhist->bin_centers[i];
        double density = counts * norm;
        fprintf(fp_hist, "%.6e,%.6e,%zu\n", centers, density, counts);
    }
    fclose(fp_hist);
}
