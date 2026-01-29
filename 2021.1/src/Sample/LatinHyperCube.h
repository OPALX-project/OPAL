//
// Class LatinHyperCube
//   This class does Latin hypercube sampling.
//
// Copyright (c) 2018, Matthias Frey, Paul Scherrer Institut, Villigen PSI, Switzerland
// All rights reserved
//
// Implemented as part of the PhD thesis
// "Precise Simulations of Multibunches in High Intensity Cyclotrons"
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
#ifndef OPAL_LATIN_HYPERCUBE_H
#define OPAL_LATIN_HYPERCUBE_H

#include "Sample/SamplingMethod.h"
#include "Sample/RNGStream.h"

#include <algorithm>
#include <deque>

class LatinHyperCube : public SamplingMethod
{

public:
    typedef typename std::uniform_real_distribution<double> dist_t;

    LatinHyperCube(double lower, double upper)
        : binsize_m(0.0)
        , upper_m(upper)
        , lower_m(lower)
        , dist_m(0.0, 1.0)
        , RNGInstance_m(RNGStream::getInstance())
        , seed_m(RNGStream::getGlobalSeed())
    {}
    
    LatinHyperCube(double lower, double upper, int seed)
        : binsize_m(0.0)
        , upper_m(upper)
        , lower_m(lower)
        , dist_m(0.0, 1.0)
        , RNGInstance_m(nullptr)
        , seed_m(seed)
    {}

    ~LatinHyperCube() {
        if ( RNGInstance_m )
            RNGStream::deleteInstance(RNGInstance_m);
    }

    void create(boost::shared_ptr<SampleIndividual>& ind, std::size_t i) {
        /* values are created within [0, 1], thus, they need to be mapped
         * the domain [lower, upper]
         */
        ind->genes[i] = map2domain_m(RNGInstance_m->getNext(dist_m));
    }
    
    void allocate(const CmdArguments_t& args, const Comm::Bundle_t& comm) {
        int id = comm.island_id;
        
        if ( !RNGInstance_m )
            RNGInstance_m = RNGStream::getInstance(seed_m + id);
        
        int nSamples = args->getArg<int>("nsamples", true);
        int nMasters = args->getArg<int>("num-masters", true);
        
        int nLocSamples = nSamples / nMasters;
        int rest = nSamples - nMasters * nLocSamples;
        
        if ( id < rest )
            nLocSamples++;
        
        int startBin = 0;
        
        if ( rest == 0 )
            startBin = nLocSamples * id;
        else {
            if ( id < rest ) {
                startBin = nLocSamples * id;
            } else {
                startBin = (nLocSamples + 1) * rest + (id - rest) * nLocSamples;
            }
        }
        
        binsize_m = ( upper_m - lower_m ) / double(nSamples);
        
        this->fillBins_m(nSamples, nLocSamples, startBin, seed_m);
    }
    
private:
    double map2domain_m(double val) {
        /* y = mx + q
         * 
         * [0, 1] --> [a, b]
         * 
         * y = (b - a) * x + a
         * 
         * where a and b are the lower, respectively, upper
         * bound of the current bin.
         */
        
        std::size_t bin = bin_m.back();
        bin_m.pop_back();
        
        return  binsize_m * (val + bin) + lower_m;
    }
    
    void fillBins_m(std::size_t nTotal, std::size_t nLocal, int startBin,
                    std::size_t seed)
    {
        std::deque<std::size_t> tmp;
        tmp.resize(nTotal);
        std::iota(tmp.begin(), tmp.end(), 0);

        // all masters need to shuffle the same way
        std::mt19937_64 eng(seed);
        std::shuffle(tmp.begin(), tmp.end(), eng);

        // each master takes its bins
        std::copy(tmp.begin()+startBin,
                  tmp.begin()+startBin+nLocal,
                  std::back_inserter(bin_m));
    }

private:
    std::deque<std::size_t> bin_m;
    double binsize_m;
    
    double upper_m;
    double lower_m;
    
    dist_t dist_m;
    
    RNGStream *RNGInstance_m;
    
    std::size_t seed_m;
};

#endif
