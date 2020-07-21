/*
 *  Copyright (c) 2012, Chris Rogers
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

#include "Elements/OpalPolynomialTimeDependence.h"

#include <string>

#include "Algorithms/PolynomialTimeDependence.h"

#include "Attributes/Attributes.h"

const std::string OpalPolynomialTimeDependence::doc_string =
    std::string("The \"POLYNOMIAL_TIME_DEPENDENCE\" element defines ")+\
    std::string("polynomial coefficients for time dependent RF phase, ")+\
    std::string("frequency, amplitude, etc, given by ")+\
    std::string("f(t) = P0+P1*t+P2*t^2+P3*t^3 where t is the time in ns");

// I investigated using a StringArray or RealArray here;
// Don't seem to have capacity to handle variables, so for now not implemented
OpalPolynomialTimeDependence::OpalPolynomialTimeDependence()
       : OpalElement(int(SIZE),
                     "POLYNOMIAL_TIME_DEPENDENCE",
                     doc_string.c_str()) {
    itsAttr[P0] = Attributes::makeReal("P0",
      "Constant term in the polynomial expansion.");
    itsAttr[P1] = Attributes::makeReal("P1",
      "First order (linear) term in the polynomial expansion.");
    itsAttr[P2] = Attributes::makeReal("P2",
      "Second order (quadratic) term in the polynomial expansion.");
    itsAttr[P3] = Attributes::makeReal("P3",
      "Third order (cubic) term in the polynomial expansion.");

    registerOwnership();
}

OpalPolynomialTimeDependence* OpalPolynomialTimeDependence::clone(const std::string &name) {
    return new OpalPolynomialTimeDependence(name, this);
}

void OpalPolynomialTimeDependence::print(std::ostream& out) const {
    OpalElement::print(out);
}

OpalPolynomialTimeDependence::OpalPolynomialTimeDependence(const std::string &name, OpalPolynomialTimeDependence *parent):
    OpalElement(name, parent) {
}

OpalPolynomialTimeDependence::~OpalPolynomialTimeDependence() {}

void OpalPolynomialTimeDependence::update() {
    // getOpalName() comes from AbstractObjects/Object.h
    std::vector<double> polynomial_coefficients;
    polynomial_coefficients.push_back(Attributes::getReal(itsAttr[P0]));
    polynomial_coefficients.push_back(Attributes::getReal(itsAttr[P1]));
    polynomial_coefficients.push_back(Attributes::getReal(itsAttr[P2]));
    polynomial_coefficients.push_back(Attributes::getReal(itsAttr[P3]));

    AbstractTimeDependence::setTimeDependence(getOpalName(),
                                              std::make_shared<PolynomialTimeDependence>(polynomial_coefficients));
}