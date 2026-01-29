//
// Class Normal
//   This class provides normally distributed samples.
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
#ifndef OPAL_NORMAL_RANDOM_SAMPLING_H
#define OPAL_NORMAL_RANDOM_SAMPLING_H

#include "Sample/SamplingMethod.h"
#include "Sample/RNGStream.h"

#include <type_traits>

class Normal : public SamplingMethod
{

public:
    typedef std::normal_distribution<double> dist_t;


    Normal(double lower, double upper)
        : dist_m(0.5 * (lower + upper), (upper - lower) / 10)
        , RNGInstance_m(RNGStream::getInstance())
        , seed_m(RNGStream::getGlobalSeed())
    {}

    Normal(double lower, double upper, std::size_t seed)
        : dist_m(0.5 * (lower + upper), (upper - lower) / 10)
        , RNGInstance_m(nullptr)
        , seed_m(seed)
    {}

    ~Normal() {
        if ( RNGInstance_m)
            RNGStream::deleteInstance(RNGInstance_m);
    }

    void create(boost::shared_ptr<SampleIndividual>& ind, size_t i) {
        ind->genes[i] = RNGInstance_m->getNext(dist_m);
    }
    
    void allocate(const CmdArguments_t& /*args*/, const Comm::Bundle_t& comm) {
        if ( !RNGInstance_m )
            RNGInstance_m = RNGStream::getInstance(seed_m + comm.island_id);
    }

private:
    dist_t dist_m;

    RNGStream *RNGInstance_m;
    
    std::size_t seed_m;
};

#endif