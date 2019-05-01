// ------------------------------------------------------------------------
// $RCSfile: OpalMultipole.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.3.4.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: OpalMultipole
//   The class of OPAL general multipoles.
//
// ------------------------------------------------------------------------
//
// $Date: 2002/12/09 15:06:07 $
// $Author: jsberg $
//
// ------------------------------------------------------------------------

#include "Elements/OpalMultipole.h"
#include "AbstractObjects/AttributeHandler.h"
#include "AbstractObjects/Expressions.h"
#include "AbstractObjects/OpalData.h"
#include "Attributes/Attributes.h"
#include "BeamlineCore/MultipoleRep.h"
#include "ComponentWrappers/MultipoleWrapper.h"
#include "Expressions/SValue.h"
#include "Expressions/SRefExpr.h"
#include "Physics/Physics.h"
#include "Utilities/Options.h"
#include <iostream>
#include <sstream>
#include <vector>


// Class OpalMultipole
// ------------------------------------------------------------------------

OpalMultipole::OpalMultipole():
    OpalElement(SIZE, "MULTIPOLE",
                "The \"MULTIPOLE\" element defines a thick multipole.\n"
                "* If the length is non-zero, the strengths are per unit "
                "length.\n* If the length is zero, the strengths are the "
                "values integrated over the length.\n"
                "* With zero length no synchrotron radiation can be calculated.") {
    itsAttr[KN] = Attributes::makeRealArray
                  ("KN", "Normalised multipole strengths (normal) in m^(-k)");
    itsAttr[DKN] = Attributes::makeRealArray
                  ("DKN", "Normalised multipole strengths errors(normal) in m^(-k)");
    itsAttr[KS] = Attributes::makeRealArray
                  ("KS", "Normalised multipole strengths (skew) in m^(-k)");
    itsAttr[DKS] = Attributes::makeRealArray
                  ("DKS", "Normalised multipole strength errors (skew) in m^(-k)");

    registerOwnership();

    setElement((new MultipoleRep("MULTIPOLE"))->makeWrappers());
}


OpalMultipole::OpalMultipole(const std::string &name, OpalMultipole *parent):
    OpalElement(name, parent) {
    setElement((new MultipoleRep(name))->makeWrappers());
}


OpalMultipole::~OpalMultipole()
{}


OpalMultipole *OpalMultipole::clone(const std::string &name) {
    return new OpalMultipole(name, this);
}


void OpalMultipole::print(std::ostream &os) const {
    OpalElement::print(os);
}


void OpalMultipole::
fillRegisteredAttributes(const ElementBase &base, ValueFlag flag) {
    OpalElement::fillRegisteredAttributes(base, flag);
    const MultipoleWrapper *mult =
        dynamic_cast<const MultipoleWrapper *>(base.removeAlignWrapper());
    BMultipoleField field;

    // Get the desired field.
    if(flag == ERROR_FLAG) {
        field = mult->errorField();
    } else if(flag == ACTUAL_FLAG) {
        field = mult->getField();
    } else if(flag == IDEAL_FLAG) {
        field = mult->getDesign().getField();
    }

    double length = getLength();
    double scale = Physics::c / OpalData::getInstance()->getP0();
    if(length != 0.0) scale *= length;

    for(int order = 1; order <= field.order(); ++order) {
        std::ostringstream ss;
        ss << (order - 1) << std::ends;
        std::string orderString = ss.str();

        std::string normName = "K" + orderString + "L";
        registerRealAttribute(normName)->setReal(scale * field.normal(order));

        std::string skewName = "K" + orderString + "SL";
        registerRealAttribute(skewName)->setReal(scale * field.skew(order));

        scale *= double(order);
    }
}


void OpalMultipole::update() {
    OpalElement::update();

    // Magnet length.
    MultipoleRep *mult =
        dynamic_cast<MultipoleRep *>(getElement()->removeWrappers());
    double length = getLength();
    mult->setElementLength(length);

    // Field components.
    BMultipoleField field;

    const std::vector<double> norm = Attributes::getRealArray(itsAttr[KN]);
    std::vector<double> normErrors = Attributes::getRealArray(itsAttr[DKN]);
    const std::vector<double> skew = Attributes::getRealArray(itsAttr[KS]);
    std::vector<double> skewErrors = Attributes::getRealArray(itsAttr[DKS]);
    int normSize = norm.size();
    int skewSize = skew.size();
    normErrors.resize(normSize, 0.0);
    skewErrors.resize(skewSize, 0.0);
    double factor = OpalData::getInstance()->getP0() / Physics::c;
    int top = (normSize > skewSize) ? normSize : skewSize;

    for(int comp = 1; comp <= top; comp++) {
        factor /= double(comp);
        if(comp <= normSize) {
            field.setNormalComponent(comp, norm[comp-1] * factor);
            mult->setNormalComponent(comp, norm[comp-1], normErrors[comp-1]);
        }
        if(comp <= skewSize) {
            field.setSkewComponent(comp, skew[comp-1] * factor);
            mult->setSkewComponent(comp, skew[comp-1], skewErrors[comp-1]);
        }
    }

    mult->setField(field);

   // Transmit "unknown" attributes.
    OpalElement::updateUnknown(mult);
}