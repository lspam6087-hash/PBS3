#ifndef VECTOR_FUNCTIONS_H
#define VECTOR_FUNCTIONS_H

#include "structs.h"
/* Header file of the functions.
   These functions contain elementary vector operations that everyone should know by now :D*/

struct Vec3D v3(double x, double y, double z);
struct Vec3D add(struct Vec3D a, struct Vec3D b);
struct Vec3D sub(struct Vec3D a, struct Vec3D b);
struct Vec3D scl(double s, struct Vec3D a);
double dot(struct Vec3D a, struct Vec3D b);
struct Vec3D cross(struct Vec3D a, struct Vec3D b);
double norm(struct Vec3D a);

#endif