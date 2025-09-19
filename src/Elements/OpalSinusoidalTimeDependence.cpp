//
// Class OpalSinusoidalTimeDependence
//   User interface for a time dependence class that generates sine waves
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
        "AMPLITUDES", "Peak-to-peak amplitude for each size wave.  If undefined, defaults to 1.0.");
    itsAttr[PHASE_OFFSETS] = Attributes::makeRealArray(
        "PHASE_OFFSETS", "Phase offset for each sine wave.  If undefined, defaults to 0.0.");
    itsAttr[DC_OFFSETS] = Attributes::makeRealArray(
        "DC_OFFSETS", "DC offset for each sine wave.  If undefined, defaults to 0.0.");
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
