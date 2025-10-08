#ifndef GRF_H_
#define GRF_H_

#include "structs.h"

void initialise_grf(struct Parameters *p_parameter, struct Vectors *p_vectors);

double update_grf(struct Parameters *p_parameter, struct Vectors *p_vectors);

double update_grf2(struct Parameters *p_parameters, struct Nbrlist *p_nbrlist, struct Vectors *p_vectors);

void finalise_grf(struct Parameters *p_parameter, struct Vectors *p_vectors, double grcount);

#endif /* GRF_H_ */
