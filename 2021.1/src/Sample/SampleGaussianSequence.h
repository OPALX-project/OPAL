//
// Class SampleGaussianSequence
//   This class provides a sequence of sampling points that have a Gaussian distribution
//   with
//         mean = 0.5 * (upper + lower)
//         sigma = (upper - lower) / 10
//   This can be achieved if the integral of the Gaussian between the sampling
//   points are all equal. The sampling points are therefore computed using
//   the inverse error function at equally distributed arguments between
//   -1 and 1.
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
#ifndef OPAL_SAMPLE_GAUSSIAN_SEQUENCE_H
#define OPAL_SAMPLE_GAUSSIAN_SEQUENCE_H

#include "Sample/SamplingMethod.h"
#include "Utilities/Util.h"

#ifdef WITH_UNIT_TESTS
#include <gtest/gtest_prod.h>
#endif

class SampleGaussianSequence : public SamplingMethod
{

public:

    SampleGaussianSequence(double lower, double upper, size_t modulo, int nSample)
        : numSamples_m(nSample)
        , volumeLowerDimensions_m(modulo)
    {
        double mean = 0.5 * (lower + upper);
        double sigma = (upper - lower) / 10; // +- 5 sigma
        double factor = sigma / sqrt(2);
        double dx = 2.0 / nSample;
        for (long i = 0; i < nSample; ++ i) {
            double x = -1.0 + (i + 0.5) * dx;
            double y = Util::erfinv(x);
            sampleChain_m.push_back(mean + factor * y);
        }
    }

    void create(boost::shared_ptr<SampleIndividual>& ind, size_t i) {
        ind->genes[i] = getNext(ind->id);
    }

    double getNext(unsigned int id) {
        int bin = int(id / volumeLowerDimensions_m) % numSamples_m;
        
        double sample = sampleChain_m[bin];
        return sample;
    }

private:
#ifdef WITH_UNIT_TESTS
    FRIEND_TEST(GaussianSampleTest, ChainTest);
#endif
    std::vector<double> sampleChain_m;
    unsigned int numSamples_m; // size of this "dimension"
    size_t volumeLowerDimensions_m; // the "volume" of the sampling space of the lower "dimensions"
};

#endif