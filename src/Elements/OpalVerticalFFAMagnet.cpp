/* 
 *  Copyright (c) 2017, Chris Rogers
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

#include "Attributes/Attributes.h"  // used?
#include "Utilities/OpalException.h"  // used?

#include "AbsBeamline/EndFieldModel/Tanh.h" // classic
#include "AbsBeamline/VerticalFFAMagnet.h" // classic
#include "Elements/OpalVerticalFFAMagnet.h"

extern Inform *gmsg;

std::string OpalVerticalFFAMagnet::docstring_m = 
      std::string("The \"VerticalFFAMagnet\" element defines a vertical FFA ")+
      std::string("magnet, which has a field that increases in the vertical ")+
      std::string("direction while maintaining similar orbits.");
OpalVerticalFFAMagnet::OpalVerticalFFAMagnet() :
    OpalElement(SIZE, "VERTICALFFAMAGNET", docstring_m.c_str()) {
    itsAttr[B0] = Attributes::makeReal
          ("B0", "The nominal dipole field of the magnet at zero height [T].");
    itsAttr[FIELD_INDEX] = Attributes::makeReal("FIELD_INDEX",
      "Exponent term in the field index [m^(-1)].");
    itsAttr[WIDTH] = Attributes::makeReal("WIDTH",
      "The full width of the magnet. Particles moving more than WIDTH/2 horizontally, in either direction, are out of the aperture.");
    itsAttr[MAX_HORIZONTAL_POWER] = Attributes::makeReal("MAX_HORIZONTAL_POWER",
      "The maximum power in horizontal coordinate that will be considered in the field expansion.");
    itsAttr[END_LENGTH] = Attributes::makeReal("END_LENGTH",
      "The end length of the FFA fringe field [m].");
    itsAttr[CENTRE_LENGTH] = Attributes::makeReal("CENTRE_LENGTH",
      "The centre length of the FFA (i.e. length of the flat top) [m].");
    itsAttr[BB_LENGTH] = Attributes::makeReal("BB_LENGTH",
      "Determines the length of the bounding box. Magnet is situated symmetrically in the bounding box. [m]");
    itsAttr[HEIGHT_POS_EXTENT] = Attributes::makeReal("HEIGHT_POS_EXTENT",
      "Height of the magnet above z=0. Particles moving upwards more than HEIGHT_POS_EXTENT are out of the aperture [m].");
    itsAttr[HEIGHT_NEG_EXTENT] = Attributes::makeReal("HEIGHT_NEG_EXTENT",
      "Height of the magnet below z=0. Particles moving downwards more than HEIGHT_NEG_EXTENT are out of the aperture [m].");
    registerRealAttribute("B0");
    registerRealAttribute("FIELD_INDEX");
    registerRealAttribute("WIDTH");
    registerRealAttribute("MAX_HORIZONTAL_POWER");
    registerRealAttribute("END_LENGTH");
    registerRealAttribute("CENTRE_LENGTH");
    registerRealAttribute("BB_LENGTH");
    registerRealAttribute("HEIGHT_NEG_EXTENT");
    registerRealAttribute("HEIGHT_POS_EXTENT");
    registerOwnership();

    VerticalFFAMagnet* magnet = new VerticalFFAMagnet("VerticalFFAMagnet");
    magnet->setEndField(new endfieldmodel::Tanh(1., 1., 1));
    setElement(magnet->makeAlignWrapper());
}


OpalVerticalFFAMagnet::OpalVerticalFFAMagnet(const std::string &name,
                                             OpalVerticalFFAMagnet *parent) :
    OpalElement(name, parent) {
    VerticalFFAMagnet* magnet = new VerticalFFAMagnet(name);
    magnet->setEndField(new endfieldmodel::Tanh(1., 1., 1));
    setElement(magnet->makeAlignWrapper());
}


OpalVerticalFFAMagnet::~OpalVerticalFFAMagnet() {
}


OpalVerticalFFAMagnet *OpalVerticalFFAMagnet::clone(const std::string &name) {
    return new OpalVerticalFFAMagnet(name, this);
}


void OpalVerticalFFAMagnet::
fillRegisteredAttributes(const ElementBase &base, ValueFlag flag) {
    OpalElement::fillRegisteredAttributes(base, flag);
}


void OpalVerticalFFAMagnet::update() {
    VerticalFFAMagnet *magnet = 
              dynamic_cast<VerticalFFAMagnet*>(getElement()->removeWrappers());
    magnet->setB0(Attributes::getReal(itsAttr[B0]));
    int maxOrder = floor(Attributes::getReal(itsAttr[MAX_HORIZONTAL_POWER]));
    magnet->setMaxOrder(maxOrder);
    magnet->setFieldIndex(Attributes::getReal(itsAttr[FIELD_INDEX]));
    magnet->setNegativeVerticalExtent(Attributes::getReal(itsAttr[HEIGHT_NEG_EXTENT]));
    magnet->setPositiveVerticalExtent(Attributes::getReal(itsAttr[HEIGHT_POS_EXTENT]));
    magnet->setBBLength(Attributes::getReal(itsAttr[BB_LENGTH]));
    magnet->setWidth(Attributes::getReal(itsAttr[WIDTH]));

    // get centre length and end length in radians
    endfieldmodel::Tanh* endField = dynamic_cast<endfieldmodel::Tanh*>(magnet->getEndField());
    double mm = 1e3;
    double end_length = Attributes::getReal(itsAttr[END_LENGTH])*mm;
    double centre_length = Attributes::getReal(itsAttr[CENTRE_LENGTH])*mm/2.;
    endField->setLambda(end_length);
    // x0 is the distance between B=0.5*B0 and B=B0 i.e. half the centre length
    endField->setX0(centre_length/2.);
    endField->setTanhDiffIndices(maxOrder+2);
    magnet->initialise();
    setElement(magnet->makeAlignWrapper());
}
