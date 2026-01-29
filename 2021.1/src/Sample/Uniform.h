//
// Class Uniform
//   This class creates uniformly distributed samples.
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
#ifndef OPAL_UNIFORM_H
#define OPAL_UNIFORM_H

#include "Sample/SamplingMethod.h"
#include "Sample/RNGStream.h"

#include <type_traits>

template <typename T>
class Uniform : public SamplingMethod
{

public:
    typedef typename std::conditional<
                        std::is_integral<T>::value,
                        std::uniform_int_distribution<T>,
                        std::uniform_real_distribution<T>
                     >::type dist_t;

    Uniform(T lower, T upper)
        : dist_m(lower, upper)
        , RNGInstance_m(RNGStream::getInstance())
        , seed_m(RNGStream::getGlobalSeed())
    {}

    Uniform(T lower, T upper, std::size_t seed)
        : dist_m(lower, upper)
        , RNGInstance_m(nullptr)
        , seed_m(seed)
    {}

    ~Uniform() {
        if ( RNGInstance_m )
            RNGStream::deleteInstance(RNGInstance_m);
    }

    void create(boost::shared_ptr<SampleIndividual>& ind, size_t i) {
        ind->genes[i] = RNGInstance_m->getNext(dist_m);
    }
    
    void allocate(const CmdArguments_t& /*args*/, const Comm::Bundle_t& comm) {
        if ( !RNGInstance_m )
            RNGInstance_m = RNGStream::getInstance(seed_m + comm.island_id);
    }

    T getNext() {
        return RNGInstance_m->getNext(dist_m);
    }

private:
    dist_t dist_m;
    
    RNGStream *RNGInstance_m;
    
    std::size_t seed_m;
};

#endif
