/*
 *  Copyright (c) 2018, Chris Rogers
 *  All rights reserved.
 *  Redistribution and use in source and binary forms, with or without
 *  modification, are permitted provided that the following conditions are met:
 *  1. Redistributions of source code must retain the above copyright notice,
 *     this list of conditions and the following disclaimer.
 *  2. Redistributions in binary form must reproduce the above copyright notice,
 *     this list of conditions and the following disclaimer in the documentation
 *     and/or other materials provided with the distribution.
 *  3. Neither the name of STFC nor the names of its contributors may be used to
 *     endorse or promote products derived from this software without specific
 *     prior written permission.
 *
 *  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 *  AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 *  IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 *  ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
 *  LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 *  CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 *  SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 *  INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 *  CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 *  ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 *  POSSIBILITY OF SUCH DAMAGE.
 */

#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>

#include "Utilities/GeneralClassicException.h"
#include "Algorithms/SplineTimeDependence.h"

SplineTimeDependence::SplineTimeDependence(size_t splineOrder,
                                           std::vector<double> times,
                                           std::vector<double> values)
    : spline_m(NULL), acc_m(NULL) {
    setSpline(splineOrder, times, values);
}

SplineTimeDependence::SplineTimeDependence(const SplineTimeDependence& /*rhs*/)
    : spline_m(NULL), acc_m(NULL) {
    setSpline(splineOrder_m, times_m, values_m);
}

SplineTimeDependence::SplineTimeDependence() : spline_m(NULL), acc_m(NULL) {
}

SplineTimeDependence::~SplineTimeDependence() {
    if (spline_m != NULL) {
        gsl_spline_free(spline_m);
    }
    if (acc_m != NULL) {
        gsl_interp_accel_free(acc_m);
    }
}

SplineTimeDependence* SplineTimeDependence::clone() {
    SplineTimeDependence* timeDep = new SplineTimeDependence();
    timeDep->setSpline(splineOrder_m, times_m, values_m);
    return timeDep;
}

Inform& SplineTimeDependence::print(Inform &os) {
    if (spline_m == NULL) {
        os << "Uninitiaised SplineTimeDependence" << endl;
        return os;
    }
    os << "SplineTimeDependence of order " << splineOrder_m
       << " with " << times_m.size() << " entries" << endl;
    return os;
}

void SplineTimeDependence::setSpline(size_t splineOrder,
              std::vector<double> times,
              std::vector<double> values) {
    if (times.size() != values.size()) {
        throw GeneralClassicException(
                            "SplineTimeDependence::SplineTimeDependence",
                            "Times and values should be of equal length");
    }
    if (times.size() < splineOrder+1) {
        throw GeneralClassicException(
                        "SplineTimeDependence::SplineTimeDependence",
                        "Times and values should be of length > splineOrder+1");
    }
    if (splineOrder != 1 and splineOrder != 3) {
        throw GeneralClassicException(
                            "SplineTimeDependence::SplineTimeDependence",
                            "Only linear or cubic interpolation is supported");
    }
    for (int i = 0; i < int(times.size())-1; ++i) {
        if (times[i] >= times[i+1]) {
            throw GeneralClassicException(
                                "SplineTimeDependence::SplineTimeDependence",
                                "Times should increase monotonically");
        }
    }
    if (spline_m != NULL) {
        gsl_spline_free(spline_m);
        spline_m = NULL;
    }
    if (splineOrder == 1) {
        spline_m = gsl_spline_alloc (gsl_interp_linear, times.size());
    } else if (splineOrder == 3) {
        spline_m = gsl_spline_alloc (gsl_interp_cspline, times.size());
    }
    times_m = times;
    values_m = values;
    gsl_spline_init(spline_m, &times[0], &values[0], times.size());
    if (acc_m == NULL) {
        acc_m = gsl_interp_accel_alloc();
    } else {
        gsl_interp_accel_reset(acc_m);
    }
}

