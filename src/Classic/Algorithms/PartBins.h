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

#ifndef OPAL_Bins_HH
#define OPAL_Bins_HH

#include <gsl/gsl_histogram.h>

#include <functional>
#include <memory>
#include <vector>

class Inform;

class PartBins {

    /**

    Definition of a bin:  low <= x < hi
    Low water mark is included, high water
    is excluded.

    */


public:


    PartBins(int bins, int sbins);

    virtual ~PartBins();


    /** \brief Add a particle to the temporary container */
    void fill(std::vector<double> &p) {
        tmppart_m.push_back(p);
        isEmitted_m.push_back(false);
    }

    /** \brief  get the number of particles in the temporary particle structure used for binning */
    size_t getNp() {
        //    size_t a = tmppart_m.size();
        // reduce(a, a, OpAddAssign());
        // return a;
        return tmppart_m.size();
    }

    /** set particles number in given bin */
    void setPartNum(int bin, long long num) {nBin_m[bin] = num;}

    /** \brief Is true if we still have particles to emit */
    bool doEmission() {return getNp() > 0;}

    bool isEmitted(int n, int /*bin*/) {
        return isEmitted_m[n]; //(isEmitted_m[n][0]==1) && (isEmitted_m[n][1] == bin);
    }

    /** assigns the proper position of particle n if it belongs to bin 'bin' */
    bool getPart(size_t n, int bin, std::vector<double> &p);

    /** sort the vector of particles such that positions of the particles decrease with increasing index.
        Then push the particles back by xmax_m + jifactor * bunch_length. In order that the method getBin(double x) works xmin_m has to be lowered a bit more.
    */
    void sortArray();

 private:
    /** assigns the particles to the bins */
    void calcHBins();
    void calcExtrema();
    /** assume we emit in monotonic increasing order */
    void setBinEmitted(int bin) {binsEmitted_m[bin] = true;}

 public:
    /** update local particles number in bin after reset Bin ID of PartBunch  */
    void resetPartInBin_cyc(size_t newPartNum[], int binID);
    /** update local particles number in bin after particle deletion */
    void updatePartInBin_cyc(size_t countLost[]);

    void setGamma(double gamma) { gamma_m = gamma;}
    double getGamma() {return gamma_m;}

protected:

    double gamma_m;
    /**
       returns the index of the bin to which the particle with z = 'x' belongs.
       If getBin returns b < 0 || b >= bins_m, then x is out of range!
    */
    int getBin(double x);

    int bins_m;
    int sBins_m;

    /** extremal particle positions */
    double xmin_m;
    double xmax_m;

    /** extremal particle position within the bins */
    std::unique_ptr<double[]> xbinmin_m;
    std::unique_ptr<double[]> xbinmax_m;

    /** bin size */
    double hBin_m;

    /** holds the particles not yet in the bunch */
    std::vector< std::vector<double> > tmppart_m;
    std::vector< bool > isEmitted_m;
    /** holds information whether all particles of a bin are emitted */
    std::unique_ptr<bool[]> binsEmitted_m;

public:

    Inform &print(Inform &os);

    int getSBins() { return sBins_m; };

    /** get the number of used bin */
    virtual int getNBins() {return gsl_histogram_bins(h_m.get()) / sBins_m; }

    /** the last emitted bin is always smaller or equal getNbins */
    int getLastemittedBin() {return nemittedBins_m; }

    /** \brief If the bunch object rebins we need to call resetBins() */
    void resetBins() {
        h_m.reset(nullptr);
    }

    virtual bool weHaveBins() {
        return h_m != nullptr;
    }

    /** \brief How many particles are in all the bins */
    size_t getTotalNum();

    /** \brief How many particles are in the bin b */
    size_t getTotalNumPerBin(int b);


protected:

    /** number of emitted bins */
    int nemittedBins_m;

    /** number of particles in the bins, the sum of all the nodes */
    std::unique_ptr<size_t[]> nBin_m;

    /** number of deleted particles in the bins */
    std::unique_ptr<size_t[]> nDelBin_m;

    std::unique_ptr<gsl_histogram> h_m;

};

inline Inform &operator<<(Inform &os, PartBins &p) {
    return p.print(os);
}


class AscendingLocationSort {
public:
    AscendingLocationSort(int direction = 0): direction_m(direction)
    {;}

    bool operator()(const std::vector<double>& first_part, const std::vector<double>& second_part) {
        return first_part[direction_m] < second_part[direction_m];
    }
private:
    int direction_m;
};

class DescendingLocationSort {
public:
    DescendingLocationSort(int direction = 0): direction_m(direction)
    {;}

    bool operator()(const std::vector<double>& first_part, const std::vector<double>& second_part) {
        return first_part[direction_m] > second_part[direction_m];
    }
private:
    int direction_m;
};

#endif // OPAL_Bins_HH
