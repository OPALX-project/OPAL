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

#ifndef _CLASSIC_SRC_ALGORITHMS_SPLINETIMEDEPENDENCE_H_
#define _CLASSIC_SRC_ALGORITHMS_SPLINETIMEDEPENDENCE_H_

#include <vector>

#include <gsl/gsl_spline.h>

#include "Utilities/GeneralClassicException.h"
#include "Algorithms/AbstractTimeDependence.h"

class Inform;

/** @class SplineTimeDependence
 * 
 *  Time dependence that follows a spline. Interpolation is supported at 
 *  linear or cubic order..
 *
 *  Interpolation is done using gsl_spline routines (1st or 3rd order).
 *
 *  Nb: I tried to use PolynomialPatch, as it is arbitrary order interpolation,
 *  but it turns out that only works for regular grids.
 */
class SplineTimeDependence : public AbstractTimeDependence {
  public:
    /** Constructor
     * 
     *  @param splineOrder the order of the spline used to fit the time 
     *         dependence; either 1 (linear interpolation) or 3 (cubic spline
     *         with quadratic smoothing) 
     *  @param times the times of successive elements in the time dependence
     *  @param values the values at each time step.
     *
     *  It is an error if times and values are not of equal length, times and
     *  values length < splineOrder + 1, or times do not increase strictly
     *  monotonically.
     */
    SplineTimeDependence(size_t splineOrder,
                         std::vector<double> times,
                         std::vector<double> values);

    /** Copy Constructor */
    SplineTimeDependence(const SplineTimeDependence& rhs);

    /** Default Constructor makes a dependence of length 2 with values 0*/
    SplineTimeDependence();
    /** Destructor does nothing */
    ~SplineTimeDependence();

    /** Return the value of the spline at a given time
     *  
     *  @param time: time at which the spline is referenced. If time is off 
     *         either end of the spline, the last few values are used to 
     *         extrapolate past the end of the spline.
     *
     */
    inline double getValue(double time);

    /** Inheritable copy constructor
     *
     *  @returns new SplineTimeDependence.
     */
    SplineTimeDependence* clone();

    /** Print summary information about the time dependence
     *  
     *  @param os "Inform" stream to which the information is printed.
     */
    Inform &print(Inform &os);
    
    /** Set the spline, deleting any existing spline data
     * 
     * 
     */
    void setSpline(size_t splineOrder,
                   std::vector<double> times,
                   std::vector<double> values);
  private:
    gsl_spline* spline_m;
    gsl_interp_accel* acc_m;
    size_t splineOrder_m;
    std::vector<double> times_m;
    std::vector<double> values_m;
  
};

inline
Inform &operator<<(Inform &os, SplineTimeDependence &p) {
  return p.print(os);
}

inline
double SplineTimeDependence::getValue(double time) {
    if (time < times_m[0] or time > times_m.back()) {
        throw GeneralClassicException(
                            "SplineTimeDependence::getValue",
                            "time out of spline range");
    }
    return gsl_spline_eval (spline_m, time, acc_m);
}

#endif

