//
// Class SamplingMethod
//   Base class for all sampling methods.
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
#ifndef OPAL_SAMPLING_METHOD_H
#define OPAL_SAMPLING_METHOD_H


#include "Sample/SampleIndividual.h"

#include <boost/smart_ptr.hpp>

#include "Comm/types.h"
#include "Util/CmdArguments.h"

class SamplingMethod
{

public:
    virtual ~SamplingMethod() {};
    virtual void create(boost::shared_ptr<SampleIndividual>& ind, size_t i) = 0;
    
    /*!
     * Allocate memory for sampling. Not every sampling method
     * requires that.
     * 
     * This function is used to reduce memory since only the
     * sampler ranks need these sampling methods.
     * 
     * @param args samler arguments
     * @param comm sampler communicator
     */
    virtual void allocate(const CmdArguments_t& /*args*/, const Comm::Bundle_t& /*comm*/) {
        /* Some sampling methods require a container.
         * In order to reduce memory only samplers should allocate
         * the memory
         */
    }
};

#endif
