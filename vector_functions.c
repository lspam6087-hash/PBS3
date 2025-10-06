#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "forces.h"
#include "constants.h"
#include "structs.h"
#include "nbrlist.h"
#include "vector_functions.h"

/* This file contains the definition of the vector operations*/

 inline struct Vec3D v3(double x, double y, double z) {
    /*Creating a 3D vector*/
    return (struct Vec3D){x, y, z};
    }

 inline struct Vec3D add(struct Vec3D a, struct Vec3D b) {
    /* Addition in 3D vectors*/
    return v3(a.x + b.x, a.y + b.y, a.z + b.z);
}

 inline struct Vec3D sub(struct Vec3D a, struct Vec3D b) {
    /* Subtraction in 3D vectors*/
    return v3(a.x - b.x, a.y - b.y, a.z - b.z);   
}

 inline struct Vec3D scl(double s, struct Vec3D a) {
    /* Scaling a 3D vector*/
    return v3(s*a.x, s*a.y, s*a.z);
}

 inline double dot(struct Vec3D a, struct Vec3D b) {
    /* Dot product between two 3D vectors*/
    return(a.x*b.x + a.y*b.y + a.z*b.z);
}

 inline struct Vec3D cross(struct Vec3D a, struct Vec3D b) {
    /* Cross product between two 3D vectors*/
    return v3(a.y*b.z - a.z*b.y, a.x*b.z - a.z*b.x, a.x*b.y - a.y*b.x);   
}

 inline double norm(struct Vec3D a) {
    /* Magnitude of a 3D vector*/
    return sqrt(a.x*a.x + a.y*a.y + a.z*a.z);
}