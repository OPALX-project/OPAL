#ifndef OPAL_VEKTOR_HH
#define OPAL_VEKTOR_HH

#include "AppTypes/Vektor.h"

typedef Vektor<double, 3> Vector_t;

/// comment: this should go to AppTypes/Vektor.h

inline
double euclidean_norm(Vector_t a) {
    return sqrt(dot(a,a));
}

#endif