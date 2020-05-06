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
#include "Algorithms/PartBinsCyc.h"
#include "Physics/Physics.h"
#include <cfloat>
#include <vector>
extern Inform *gmsg;

// constructer function for cyclotron
PartBinsCyc::PartBinsCyc(int specifiedNumBins, int bins, size_t  partInBin[])
    : PartBins(specifiedNumBins, 0) {

    bins_m = specifiedNumBins;        // max bin number
    nemittedBins_m = bins;            // the bin number with particles
    
    for(int i = 0; i < nemittedBins_m; i++) {
        nBin_m[i] = partInBin[i];

        *gmsg << "Read in: Bin=" << i << " Particles Num=" << getTotalNumPerBin(i) << endl;
        binsEmitted_m[i] = true;
    }
}

// constructer function for cyclotron for restart run.
PartBinsCyc::PartBinsCyc(int specifiedNumBins, int bins)
    : PartBins(specifiedNumBins, 0) {

    bins_m = specifiedNumBins;        // max bin number
    nemittedBins_m = bins;            // the bin number with particles
    
    for(int i = 0; i < nemittedBins_m; i++) {
      binsEmitted_m[i] = true;
    }
}
