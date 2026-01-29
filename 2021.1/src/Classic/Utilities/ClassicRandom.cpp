// ------------------------------------------------------------------------
// $RCSfile: Random.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: Random:
//   A standard random generator.
// ------------------------------------------------------------------------
// Class category: Utilities
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:32:38 $
// $Author: fci $
//
// ------------------------------------------------------------------------

#include "Utilities/ClassicRandom.h"
#include <cmath>


// using a default seed value, initialize the pseudo-random number
// generator and load the array irn with a set of random numbers in
// [0, maxran)
// input: a seed value in the range [0..maxran)
Random::Random() {
    init55(123456789);
}


// using the given seed value, initialize the pseudo-random number
// generator and load the array irn with a set of random numbers in
// [0, maxran)
// input: a seed value in the range [0..maxran)
Random::Random(int seed) {
    init55(seed);
}


// destructor for random generator.
Random::~Random()
{}


// set new seed.
void Random::reseed(int seed) {
    init55(seed);
}


// return a real pseudo-random number in the range [0,1).
double Random::uniform() {
    const double scale = 1.0 / maxran;

    if(next >= nr) irngen();

    return scale * (double) irn[next++];
}


// return two double Gaussioan pseudo-random numbers
void Random::gauss(double &gr1, double &gr2) {
    double xi1, xi2, zzr;

    // construct two random numbers in range [-1.0,+1.0)
    // reject, if zzr > 1
    do {
        xi1 = uniform() * 2.0 - 1.0;
        xi2 = uniform() * 2.0 - 1.0;
        zzr = xi1 * xi1 + xi2 * xi2;
    } while(zzr > 1.0);

    // transform accepted point to Gaussian distribution
    zzr = sqrt(-2.0 * log(zzr) / zzr);
    gr1 = xi1 * zzr;
    gr2 = xi2 * zzr;
}


// return an integer pseudo-random number in the range [0,maxran]
int Random::integer() {
    if(next > nr) irngen();

    return irn[next++];
}


// initializes the random generator after setting a seed.
void Random::init55(int seed) {
    static const int nd = 21;
    next = 0;

    if(seed < 0)
        seed = - seed;
    else if(seed == 0)
        seed = 1234556789;

    int j = seed % maxran;
    irn[nr-1] = j;
    int k = 1;

    for(int i = 1; i < nr; i++) {
        int ii = (nd * i) % nr;
        irn[ii-1] = k;
        k = j - k;

        if(k < 0) k += maxran;

        j = irn[ii-1];
    }

    // call irngen a few times to "warm it up"
    irngen();
    irngen();
    irngen();
}


// generate the next "nr" elements in the pseudo-random sequence.
void Random::irngen() {
    static const int nj = 24;
    int i, j;

    for(i = 0; i < nj; i++) {
        j = irn[i] - irn[i+nr-nj];

        if(j < 0) j += maxran;

        irn[i] = j;
    }

    for(i = nj; i < nr; i++) {
        j = irn[i] - irn[i-nj];

        if(j < 0) j += maxran;

        irn[i] = j;
    }

    next = 0;
}
