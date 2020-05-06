//
// Class PartBinsCyc
//   Defines a structure to hold energy bins and their
//   associated data for multi-bunch tracking in cyclotrons
//
// Copyright (c) 2010, Jianjun Yang, Paul Scherrer Institut, Villigen PSI, Switzerland
// Copyright (c) 2017-2019, Paul Scherrer Institut, Villigen PSI, Switzerland
// All rights reserved
//
// Implemented as part of the PhD thesis
// "Beam dynamics in high intensity cyclotrons including neighboring bunch effects"
// and the paper
// "Beam dynamics in high intensity cyclotrons including neighboring bunch effects:
// Model, implementation, and application"
// (https://journals.aps.org/prab/pdf/10.1103/PhysRevSTAB.13.064201)
//
// This file is part of OPAL.
//
// OPAL is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// You should have received a copy of the GNU General Public License
// along with OPAL. If not, see <https://www.gnu.org/licenses/>.
//
#ifndef OPAL_BinsCyc_HH
#define OPAL_BinsCyc_HH

#ifndef PartBinTest
#else
#include "ranlib.h"
#define Inform ostream
#endif

#include "Algorithms/PartBins.h"
#include <gsl/gsl_rng.h>
#include <gsl/gsl_histogram.h>
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_randist.h>

#include <cassert>

class PartBinsCyc: public PartBins {

public:


    /** constructer function for cyclotron*/
    PartBinsCyc(int bunches, int bins, size_t  partInBin[]);
    PartBinsCyc(int specifiedNumBins, int bins);

    /** get the number of used bin */
    int getNBins() {return bins_m; }

    /** \brief How many particles are on one bin */
    inline size_t getGlobalBinCount(int bin) {
      size_t a = nBin_m[bin];
      reduce(a, a, OpAddAssign());
      return a;
    }

    /** \brief How many particles are on one energy bin */
    inline size_t getLocalBinCount(int bin) {
        assert(bin < bins_m);
        return nBin_m[bin];
    }

    bool weHaveBins() {
      return  ( nemittedBins_m > 0 );
    }
};

#endif // OPAL_BinsCyc_HH
