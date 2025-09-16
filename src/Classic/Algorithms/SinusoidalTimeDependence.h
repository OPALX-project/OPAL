/*
 *  Copyright (c) 2025, Jon Thompson
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

#ifndef OPAL_SINUSOIDALTIMEDEPENDENCE_H
#define OPAL_SINUSOIDALTIMEDEPENDENCE_H

#include <vector>
#include <cmath>

#include "Algorithms/AbstractTimeDependence.h"
class Inform;

/** @class SinusoidalTimeDependence
 *
 *  Time dependence that follows sum of a set of sinusoids
 *      sigma_over_i(a[i] / 2 * sin(2 * pi * f[i] * t + p[i]) + o[i])
 *  a is the peak to peak amplitude, f is the frequency, p is the phase offset,
 *  o is the DC offset, t is the time.
 */
class SinusoidalTimeDependence final : public AbstractTimeDependence {
public:
    /** Constructor
     *  @param f the frequencies in Hz; can be of arbitrary length
     *  @param p the phase offsets in radians; can be of arbitrary length
     *  @param a the peak-to-peak amplitude; can be of arbitrary length
     *  @param o the DC offset; can be of arbitrary length
     *  (user is responsible for issues like floating point precision).
     */
    SinusoidalTimeDependence(
        const std::vector<double>& f, const std::vector<double>& p, const std::vector<double>& a,
        const std::vector<double>& o);

    /** Default Constructor makes a 0 length polynomial */
    SinusoidalTimeDependence() = default;

    /** Destructor does default */
    ~SinusoidalTimeDependence() override = default;

    /** Return the sinusoidal value
     *  @param time simulation time in seconds
     */
    double getValue(double time) override;

    /** Inheritable clone function
     *  @returns new SinusoidalTimeDependence that is a copy of this. User owns
     *  returned memory.
     */
    SinusoidalTimeDependence* clone() override;

    /** Print the sinusoidals
     *  @param os "Inform" stream to which the sinusoidals are printed.
     */
    Inform& print(Inform& os) const;

private:
    std::vector<double> f_m;
    std::vector<double> p_m;
    std::vector<double> a_m;
    std::vector<double> o_m;
};

inline Inform& operator<<(Inform& os, const SinusoidalTimeDependence& p) {
    return p.print(os);
}

#endif  // OPAL_SINUSOIDALTIMEDEPENDENCE_H
