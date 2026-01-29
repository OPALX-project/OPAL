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

#include <string>
#include "Attributes/Attributes.h"
#include "Utilities/OpalException.h"
#include "Elements/OpalSplineTimeDependence.h"

const std::string OpalSplineTimeDependence::doc_string =
    std::string("The \"SPLINE_TIME_DEPENDENCE\" element defines ")+\
    std::string("an array of times and corresponding values for time lookup, ")+\
    std::string("for use in time-dependent elements. Lookup is supported at ")+\
    std::string("first order or third order with quadratic smoothing.");

// I investigated using a StringArray or RealArray here;
// Don't seem to have capacity to handle variables, so for now not implemented
OpalSplineTimeDependence::OpalSplineTimeDependence()
       : OpalElement(int(SIZE),
                     "SPLINE_TIME_DEPENDENCE",
                     doc_string.c_str()) {
    itsAttr[ORDER] = Attributes::makeReal("ORDER",
      std::string("Order of the lookup - either 1 for linear interpolation, ")+
      std::string("or 3 for cubic interpolation with quadratic smoothing. ")+
      std::string("Other values make an error."));

    itsAttr[TIMES] = Attributes::makeRealArray("TIMES",
      std::string("Array of real times in ns. There must be at least \"ORDER\"+1 ")+
      std::string("elements in the array and they must be strictly monotically ")+
      std::string("increasing"));

    itsAttr[VALUES] = Attributes::makeRealArray("VALUES",
      std::string("Array of real values. The length of \"VALUES\" must be the ")+
      std::string("same as the length of \"TIMES\"."));

    registerOwnership();
}

OpalSplineTimeDependence* OpalSplineTimeDependence::clone(const std::string &name) {
    return new OpalSplineTimeDependence(name, this);
}

void OpalSplineTimeDependence::print(std::ostream& out) const {
    OpalElement::print(out);
}

OpalSplineTimeDependence::OpalSplineTimeDependence(const std::string &name,
                                                   OpalSplineTimeDependence *parent):
    OpalElement(name, parent) {
}

OpalSplineTimeDependence::~OpalSplineTimeDependence() {}

void OpalSplineTimeDependence::update() {

    double orderReal = Attributes::getReal(itsAttr[ORDER])+1e-10;
    if ((orderReal - 1.) > 1e-9 && (orderReal - 3.) > 1e-9) {
        throw OpalException("OpalSplineTimeDependence::update",
                            "SPLINE_TIME_DEPENDENCE \"ORDER\" should be 1 or 3.");
    }
    size_t order(floor(orderReal));
    std::vector<double> times = Attributes::getRealArray(itsAttr[TIMES]);
    std::vector<double> values = Attributes::getRealArray(itsAttr[VALUES]);
    std::shared_ptr<SplineTimeDependence> spline =
                                  std::make_shared<SplineTimeDependence>(order,
                                                                        times,
                                                                        values);
    AbstractTimeDependence::setTimeDependence(getOpalName(), spline);
}
