//
// Class SampleSequence
//   This class provides a sequence of equidistant sampling points. It
//   can't be garanteed that the sampling is equidistant if
//   an integer type is chosen and the difference between
//   the upper and lower limit isn't divisible by the number
//   of sampling points.
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
#ifndef OPAL_SAMPLE_SEQUENCE_H
#define OPAL_SAMPLE_SEQUENCE_H

#include "Sample/SamplingMethod.h"

template <typename T>
class SampleSequence : public SamplingMethod
{

public:

    SampleSequence(T lower, T upper, size_t modulo, int nSample)
        : lowerLimit_m(lower)
        , stepSize_m( (upper - lower) / double(nSample - 1) )
        , numSamples_m(nSample)
        , volumeLowerDimensions_m(modulo)
    { }

    void create(boost::shared_ptr<SampleIndividual>& ind, size_t i) {
        
        unsigned int id = ind->id;
        
        int bin = int(id / volumeLowerDimensions_m) % numSamples_m;
        
        ind->genes[i] = static_cast<T>(lowerLimit_m + stepSize_m * bin);
    }
    
private:
    T lowerLimit_m;
    double stepSize_m;
    unsigned int numSamples_m; // size of this "dimension"
    size_t volumeLowerDimensions_m; // the "volume" of the sampling space of the lower "dimensions"
};

#endif
