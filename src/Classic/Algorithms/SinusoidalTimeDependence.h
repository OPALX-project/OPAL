//
// Class SinusoidalTimeDependence
//   A time dependence class that generates sine waves
//
// Copyright (c) 2025, Jon Thompson, STFC Rutherford Appleton Laboratory, Didcot, UK
//
// This file is part of OPAL.
//
// OPAL is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// You should have received a copy of the GNU General Public License
// along with OPAL.  If not, see <https://www.gnu.org/licenses/>.
//

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
