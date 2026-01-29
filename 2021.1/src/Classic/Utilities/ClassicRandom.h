#ifndef MAD_Random_HH
#define MAD_Random_HH

// ------------------------------------------------------------------------
// $RCSfile: ClassicRandom.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: Random
//
// ------------------------------------------------------------------------
// Class category: Utilities
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:32:38 $
// $Author: fci $
//
// ------------------------------------------------------------------------

// the random numbers are generated in batches of 55
const int nr = 55;


// the random integers are generated in the range [0, MAXRAN)
static const int maxran = 1000000000;

/// The CLASSIC random generator.
// This generator is based on:
// [center]
// D. Knuth: The Art of Computer Programming, Vol. 2, 2nd edition.
// [/center]

class Random {

public:

    /// Constructor with default seed.
    Random();

    /// Constructor with given seed.
    //  The argument is the new seed, an integer in the range [0, 1000000000].
    Random(int seed);

    ~Random();

    /// Set a new seed.
    //  The argument is the new seed, an integer in the range [0, 1000000000].
    void reseed(int seed = 123456789);

    /// Uniform distribution.
    //  Return a real pseudo-random number in the range [0,1).
    double uniform();

    /// Gaussian distribution.
    //  Return two real Gaussian pseudo-random numbers with unit standard
    //  deviation.  The values are placed in the two arguments.
    void gauss(double &gr1, double &gr2);

    /// Uniform distribution.
    //  Return an integer pseudo-random number in the range [0,1000000000).
    int integer();

    /// Initialise random number generator.
    void init55(int seed);

private:

    // generate next batch of random numbers
    void irngen();

    // a set of "NR" random integers
    int irn[nr];

    // the next item to be used in the random sequence
    int next;
};

#endif // MAD_Random_HH
