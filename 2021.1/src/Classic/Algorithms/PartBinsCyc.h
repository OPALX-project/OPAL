//
// Class PartBinsCyc
//   Defines a structure to hold energy bins and their
//   associated data for multi-bunch tracking in cyclotrons
//
// Copyright (c) 2010, Jianjun Yang, Paul Scherrer Institut, Villigen PSI, Switzerland
// Copyright (c) 2017-2020, Paul Scherrer Institut, Villigen PSI, Switzerland
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

#include "Algorithms/PartBins.h"

class PartBinsCyc: public PartBins {

public:


    /** constructor function for cyclotron*/
    PartBinsCyc(int bunches, int bins, size_t  partInBin[]);
    PartBinsCyc(int specifiedNumBins, int bins);

    /** get the number of used bin */
    virtual int getNBins() {return bins_m; }

    virtual bool weHaveBins() {
        return ( nemittedBins_m > 0 );
    }
};

#endif // OPAL_BinsCyc_HH
