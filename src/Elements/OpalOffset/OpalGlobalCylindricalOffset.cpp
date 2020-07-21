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

#include "Elements/OpalOffset/OpalGlobalCylindricalOffset.h"

#include <string>

#include "AbsBeamline/Offset.h"
#include "Elements/OpalElement.h"

#include "Attributes/Attributes.h"

namespace OpalOffset {

const std::string OpalGlobalCylindricalOffset::doc_string =
    std::string("The \"LOCAL_CYLINDRICAL_OFFSET\" element defines an offset")+
    std::string("in cylindrical coordinates, relative to the last placed ")+
    std::string("element. All angles are defined in the midplane.");

OpalGlobalCylindricalOffset::OpalGlobalCylindricalOffset()
       : OpalElement(int(SIZE),
                     "GLOBAL_CYLINDRICAL_OFFSET",
                     doc_string.c_str()) {
    itsAttr[RADIUS] = Attributes::makeReal("RADIUS",
             "Angle between the previous element and the displacement vector.");
    itsAttr[AZIMUTHAL_ANGLE] = Attributes::makeReal("AZIMUTHAL_ANGLE",
             "Angle between the displacement vector and the next element.");
    itsAttr[TANGENTIAL_OFFSET] = Attributes::makeReal("TANGENTIAL_OFFSET",
             "Length of the offset.");

    registerOwnership();
}

OpalGlobalCylindricalOffset* OpalGlobalCylindricalOffset::clone(const std::string &name) {
    return new OpalGlobalCylindricalOffset(name, this);
}

void OpalGlobalCylindricalOffset::print(std::ostream& out) const {
    OpalElement::print(out);
}

OpalGlobalCylindricalOffset::OpalGlobalCylindricalOffset(const std::string &name, OpalGlobalCylindricalOffset *parent):
    OpalElement(name, parent) {
}

OpalGlobalCylindricalOffset::~OpalGlobalCylindricalOffset() {}

void OpalGlobalCylindricalOffset::update() {
    // getOpalName() comes from AbstractObjects/Object.h
    std::string name = getOpalName();
    double radius = Attributes::getReal(itsAttr[RADIUS]);
    double phi = Attributes::getReal(itsAttr[AZIMUTHAL_ANGLE]);
    double theta = Attributes::getReal(itsAttr[TANGENTIAL_OFFSET]);
    Offset* off = new Offset(Offset::globalCylindricalOffset(name,
                                                            radius,
                                                            phi,
                                                            theta));
    // is this a memory leak?
    setElement(off);
}
}