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

#include "Elements/OpalSBend3D.h"
#include "Physics/Physics.h"
#include "AbsBeamline/SBend3D.h"
#include "Attributes/Attributes.h"
#include "Utilities/OpalException.h"

extern Inform *gmsg;

OpalSBend3D::OpalSBend3D():
    OpalElement(SIZE, "SBEND3D",
             "The \"SBEND3D\" element defines a sector bending magnet.") {
    itsAttr[FMAPFN] = Attributes::makeString
                                       ("FMAPFN", "The name of the field map.");
    itsAttr[FIELD_UNITS] = Attributes::makeReal("FIELD_UNITS",
            "Scale the field map up or down by this factor (default 1. is T).");
    itsAttr[LENGTH_UNITS] = Attributes::makeReal("LENGTH_UNITS",
            "Units for length coordinates (default 1. is mm).");

    registerOwnership();

    setElement(new SBend3D("SBEND3D"));
}


OpalSBend3D::OpalSBend3D(const std::string &name, OpalSBend3D *parent):
    OpalElement(name, parent) {
    setElement(new SBend3D(name));
}


OpalSBend3D::~OpalSBend3D() {
}


OpalSBend3D *OpalSBend3D::clone(const std::string &name) {
    return new OpalSBend3D(name, this);
}


void OpalSBend3D::update() {
    SBend3D *bend = dynamic_cast<SBend3D*>(getElement());
    if (itsAttr[FIELD_UNITS])
        bend->setFieldUnits(Attributes::getReal(itsAttr[FIELD_UNITS]));
    if (itsAttr[LENGTH_UNITS])
        bend->setLengthUnits(Attributes::getReal(itsAttr[LENGTH_UNITS]));
    // this has to be done last as we initialise field map here
    // (need units before initialisation)
    bend->setFieldMapFileName(Attributes::getString(itsAttr[FMAPFN]));
    setElement(bend);

}