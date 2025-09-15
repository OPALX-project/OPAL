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

#include "OpalSinusoidalTimeDependence.h"
#include <string>
#include "Algorithms/SinusoidalTimeDependence.h"
#include "Attributes/Attributes.h"

const std::string OpalSinusoidalTimeDependence::doc_string =
    "The \"SINUSOIDAL_TIME_DEPENDENCE\" element defines "
    "sinusoidal coefficients for a time dependence, "
    "frequency, amplitude, phase offset, DC offset, given by "
    "f(t) = sigma_over_i(a[i] / 2 * sin(2 * pi * f[i] * t + p[i]) + o[i])";

OpalSinusoidalTimeDependence::OpalSinusoidalTimeDependence()
    : OpalElement(static_cast<int>(SIZE), "SINUSOIDAL_TIME_DEPENDENCE", doc_string.c_str()) {
    itsAttr[FREQUENCIES] = Attributes::makeRealArray(
        "FREQUENCIES",
        "Sine wave frequencies, length determines the number of size waves included.");
    itsAttr[AMPLITUDES] = Attributes::makeRealArray(
        "AMPLITUDES", "Peak-to-peak amplitude for each size wave, if missing 1.0 is assumed.");
    itsAttr[PHASE_OFFSETS] = Attributes::makeRealArray(
        "PHASE_OFFSETS", "Phase offset for each sine wave, if missing 0.0 is assumed.");
    itsAttr[DC_OFFSETS] = Attributes::makeRealArray(
        "DC_OFFSETS", "DC offset for each sine wave, if missing 0.0 is assumed.");
    registerOwnership();
}

OpalSinusoidalTimeDependence* OpalSinusoidalTimeDependence::clone(const std::string& name) {
    return new OpalSinusoidalTimeDependence(name, this);
}

void OpalSinusoidalTimeDependence::print(std::ostream& out) const {
    OpalElement::print(out);
}

OpalSinusoidalTimeDependence::OpalSinusoidalTimeDependence(
    const std::string& name, OpalSinusoidalTimeDependence* parent)
    : OpalElement(name, parent) {
}

void OpalSinusoidalTimeDependence::update() {
    AbstractTimeDependence::setTimeDependence(
        getOpalName(), std::make_shared<SinusoidalTimeDependence>(
                           Attributes::getRealArray(itsAttr[FREQUENCIES]),
                           Attributes::getRealArray(itsAttr[PHASE_OFFSETS]),
                           Attributes::getRealArray(itsAttr[AMPLITUDES]),
                           Attributes::getRealArray(itsAttr[DC_OFFSETS])));
}
