//
// Class PartBins
//   Defines a structure to hold energy bins and their associated data
//
// Copyright (c) 2007-2020, Paul Scherrer Institut, Villigen PSI, Switzerland
// All rights reserved
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

#include "Algorithms/PartBins.h"
#include "Algorithms/PBunchDefs.h"
#include "Physics/Physics.h"
#include "Physics/Units.h"

#include "Utility/Inform.h"

#include <cfloat>
#include <limits>
#include <vector>
#include <cmath>

PartBins::PartBins(int bins, int sbins) :
    gamma_m(1.0),
    bins_m(bins),
    sBins_m(sbins),
    xmin_m(0.0),
    xmax_m(0.0),
    hBin_m(0.0),
    nemittedBins_m(0) {

    // number of particles in the bins on the local node
    nBin_m        = std::unique_ptr<size_t[]>(new size_t[bins_m]);
    xbinmin_m     = std::unique_ptr<double[]>(new double[bins_m]);
    xbinmax_m     = std::unique_ptr<double[]>(new double[bins_m]);

    // flag whether the bin contain particles or not
    binsEmitted_m = std::unique_ptr<bool[]>(new bool[bins_m]);

    nDelBin_m     = std::unique_ptr<size_t[]>(new size_t[bins_m]);

    for(int i = 0; i < bins_m; i++) {
        nDelBin_m[i] = nBin_m[i] = 0;
        xbinmin_m[i] = std::numeric_limits<double>::max();
        xbinmax_m[i] = -xbinmin_m[i];
        binsEmitted_m[i] = false;
    }
}


size_t PartBins::getTotalNum() {
    size_t s = 0;
    size_t sd = 0;
    size_t gs = 0;

    for(int i = 0; i < getLastemittedBin(); i++) {
        s  += nBin_m[i];
        sd += nDelBin_m[i];
    }
    gs = s - sd;
    reduce(gs, gs, OpAddAssign());
    return gs;
}

size_t PartBins::getTotalNumPerBin(int b) {
    size_t s = 0;
    s  = nBin_m[b];
    reduce(s, s, OpAddAssign());
    return s;
}


void PartBins::updatePartInBin_cyc(size_t countLost[]) {

  for(int ii = 0; ii < nemittedBins_m; ii++) {
    if(countLost[ii] > 0)
      nBin_m[ii] -= countLost[ii];
  }
}

void PartBins::resetPartInBin_cyc(size_t newPartNum[], int maxbinIndex) {
    reduce(maxbinIndex, maxbinIndex, OpMaxAssign());
    nemittedBins_m =  maxbinIndex + 1;

    for(int ii = 0; ii < nemittedBins_m; ii++) {
        nBin_m[ii] = newPartNum[ii]; // only count particles on the local node
        setBinEmitted(ii);  // set true for this bin
    }
}


PartBins::~PartBins() {
    tmppart_m.clear();
    isEmitted_m.clear();
}


bool PartBins::getPart(size_t n, int bin, std::vector<double> &p) {

    if(tmppart_m[n][6] == bin) {
        p = tmppart_m[n];
        return true;
    } else
        return false;
}

/** /brief There is only a local sort, no global yet */
void PartBins::sortArray() {
    /** sort the vector of particles such that position of the particles decrease with increasing index.
        Then push the particles back by 1e-13 s * beta * c (approximately one step).
        In order that the method getBin(double x) works xmin_m has to be lowered a bit more.
    */

    double sshift = std::sqrt(1. - (1. / (gamma_m * gamma_m))) * Physics::c * 0.1 * Units::ps2s;
    std::sort(tmppart_m.begin(), tmppart_m.end(), DescendingLocationSort(2));
    xmax_m = tmppart_m[0][2];
    xmin_m = tmppart_m.back()[2];

    for(unsigned int n = 0; n < tmppart_m.size(); n++)
        tmppart_m[n][2] -= xmax_m + sshift; /* push particles back */

    xmin_m -= xmax_m + 0.0001 * (xmax_m - xmin_m) + sshift; /* lower the limits */
    xmax_m = -sshift;

    reduce(xmin_m, xmin_m, OpMinAssign());
    reduce(xmax_m, xmax_m, OpMaxAssign());

    hBin_m = (std::abs(xmax_m - xmin_m)) / (bins_m);
    calcHBins();
    for(int n = 0; n < bins_m; n++)
        if(nBin_m[n] == 0) setBinEmitted(n);
}


void PartBins::calcHBins() {

    for(unsigned int n = 0; n < tmppart_m.size(); n++)
        tmppart_m[n][6] = getBin(tmppart_m[n][2]);
    calcExtrema();
}

void PartBins::calcExtrema() {
    for(unsigned int n = 0; n < tmppart_m.size(); n++) {
        if(xbinmin_m[(int)tmppart_m[n][6]] >= tmppart_m[n][2])
            xbinmin_m[(int)tmppart_m[n][6]] = tmppart_m[n][2];

        if(xbinmax_m[(int)tmppart_m[n][6]] <= tmppart_m[n][2])
            xbinmax_m[(int)tmppart_m[n][6]] = tmppart_m[n][2];
    }
}

Inform &PartBins::print(Inform &os) {

    os << "-----------------------------------------" << endl;
    os << "     CREATE BINNED GAUSS DISTRIBUTION DONE        " << endl;

    os << "Bins= " << bins_m << " hBin= " << hBin_m << " Particle vector length " << tmppart_m.size() << endl;

    //for(int i = 0; i < gsl_histogram_bins(h_m); i++)
    //os << "Bin # " << i << " val " << gsl_histogram_get(h_m, i) << endl;
    for(int i = 0; i < bins_m; i++) {
        size_t msum = 0;
        for(int j = 0; j < sBins_m; j++)
            msum += gsl_histogram_get(h_m.get(), i * sBins_m + j);
        os << "Bin # " << i << " val " << msum << endl;
    }

    if(getLastemittedBin() >= 0)
        os << "Last emitted bin is " << getLastemittedBin() << endl;
    else
        os << "No bin is emitted !" << endl;
    return os;
}

int PartBins::getBin(double x) {
    /**
       returns the index of the bin to which the particle with z = 'x' belongs.
       If getBin returns b < 0 || b >= bins_m, then is x out of range!
    */
    int b = (int) std::floor(std::abs(xmax_m - x) / hBin_m);
    nBin_m[b]++;
    return b;
}