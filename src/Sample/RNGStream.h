//
// Class RNGStream
//   This class takes care of RNG generator instances.
//
// Copyright (c) 2018, Christof Metzger-Kraus, Open Sourcerer
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
#ifndef RNGSTREAM_H
#define RNGSTREAM_H

#include <random>

class RNGStream
{
public:
    static RNGStream* getInstance();
    static RNGStream* getInstance(unsigned int seed);
    static void deleteInstance(RNGStream* & generator);

    static void setGlobalSeed(unsigned int seed);

    static unsigned int getGlobalSeed();
    
    std::mt19937_64 & getGenerator();

    template <class DISTR>
    typename DISTR::result_type getNext(DISTR & RNGDist) {
        return RNGDist(RNGenerator_m);
    }
    
private:
    RNGStream():
        RNGenerator_m(globalSeed_sm),
        isGlobal_m(true)
    { }

    RNGStream(unsigned int seed):
        RNGenerator_m(seed),
        isGlobal_m(false)
    { }

    ~RNGStream()
    { }

    static RNGStream *globalInstance_sm;
    static unsigned int globalSeed_sm;
    static unsigned int numGlobalInstances_sm;
    std::mt19937_64 RNGenerator_m;
    bool isGlobal_m;

};
#endif