#ifndef GRF_H_
#define GRF_H_

#include "structs.h"

void initialise_grf(struct Parameters *p_parameter, struct Vectors *p_vectors);

void update_grf(struct Parameters *p_parameter, struct Vectors *p_vectors);

void finalise_grf(struct Parameters *p_parameter, struct Vectors *p_vectors);

#endif /* GRF_H_ */
