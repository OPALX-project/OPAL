// ------------------------------------------------------------------------
// $RCSfile: OpalElement.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.2.4.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: OpalElement
//   The base class for all OPAL beamline elements.

//   and for printing in OPAL-8 format.
//
// ------------------------------------------------------------------------
//
// $Date: 2002/12/09 15:06:07 $
// $Author: jsberg $
//
// ------------------------------------------------------------------------

#include "Elements/OpalElement.h"
#include "AbsBeamline/AlignWrapper.h"
#include "AbsBeamline/ElementImage.h"
#include "AbsBeamline/Bend.h"
#include "AbstractObjects/Attribute.h"
#include "AbstractObjects/Expressions.h"
#include "AbstractObjects/OpalData.h"
#include "Attributes/Attributes.h"
#include "Parser/Statement.h"
#include "Physics/Physics.h"
#include "Utilities/OpalException.h"
#include "Utilities/Options.h"
#include "Utilities/ParseError.h"
#include "Utilities/Round.h"
#include "Utilities/Util.h"

#include <cmath>
#include <cctype>
#include <sstream>
#include <vector>
#include <boost/regex.hpp>

extern Inform *gmsg;
// Class OpalElement
// ------------------------------------------------------------------------

std::map < std::string, OwnPtr<AttCell> > OpalElement::attributeRegistry;

OpalElement::OpalElement(int size, const char *name, const char *help):
    Element(size, name, help), itsSize(size) {
    itsAttr[TYPE]   = Attributes::makeString
                      ("TYPE", "The element design type (the project name)");
    itsAttr[LENGTH] = Attributes::makeReal
                      ("L", "The element length in m");
    itsAttr[ELEMEDGE] = Attributes::makeReal
                        ("ELEMEDGE", "The position of the element in path length (in m)");
    itsAttr[APERT]  = Attributes::makeString
                      ("APERTURE", "The element aperture");
    itsAttr[WAKEF]   = Attributes::makeString
                       ("WAKEF", "Defines the wake function");
    itsAttr[PARTICLEMATTERINTERACTION]   = Attributes::makeString
                                ("PARTICLEMATTERINTERACTION", "Defines the particle mater interaction handler");
    itsAttr[ORIGIN] = Attributes::makeRealArray
                      ("ORIGIN", "The location of the element");

    itsAttr[ORIENTATION] = Attributes::makeRealArray
                           ("ORIENTATION", "The Tait-Bryan angles for the orientation of the element");

    itsAttr[X] = Attributes::makeReal
        ("X", "The x-coordinate of the location of the element", 0);

    itsAttr[Y] = Attributes::makeReal
        ("Y", "The y-coordinate of the location of the element", 0);

    itsAttr[Z] = Attributes::makeReal
        ("Z", "The z-coordinate of the location of the element", 0);

    itsAttr[THETA] = Attributes::makeReal
        ("THETA", "The rotation about the y-axis of the element", 0);

    itsAttr[PHI] = Attributes::makeReal
        ("PHI", "The rotation about the x-axis of the element", 0);

    itsAttr[PSI] = Attributes::makeReal
        ("PSI", "The rotation about the z-axis of the element", 0);

    itsAttr[DX] = Attributes::makeReal
                  ("DX", "Misalignment in x direction",0.0);
    itsAttr[DY] = Attributes::makeReal
                  ("DY", "Misalignment in y direction",0.0);
    itsAttr[DZ] = Attributes::makeReal
                  ("DZ", "Misalignment in z direction",0.0);
    itsAttr[DTHETA] = Attributes::makeReal
                  ("DTHETA", "Misalignment in theta (Tait-Bryan angles)",0.0);
    itsAttr[DPHI] = Attributes::makeReal
                  ("DPHI", "Misalignment in theta (Tait-Bryan angles)",0.0);
    itsAttr[DPSI] = Attributes::makeReal
                  ("DPSI", "Misalignment in theta (Tait-Bryan angles)",0.0);

    const unsigned int end = COMMON;
    for (unsigned int i = 0; i < end; ++ i) {
        AttributeHandler::addAttributeOwner("Any", AttributeHandler::ELEMENT, itsAttr[i].getName());
    }

    static bool first = true;
    if(first) {
        registerStringAttribute("NAME");
        registerStringAttribute("TYPE");
        registerStringAttribute("CLASS");
        registerStringAttribute("KEYWORD");
        registerRealAttribute("L");
        registerStringAttribute("WAKEF");
        registerStringAttribute("PARTICLEMATTERINTERACTION");
        registerStringAttribute("APERT");
        registerRealAttribute("X");
        registerRealAttribute("Y");
        registerRealAttribute("Z");
        registerRealAttribute("THETA");
        registerRealAttribute("PHI");
        registerRealAttribute("PSI");
        registerRealAttribute("DX");
        registerRealAttribute("DY");
        registerRealAttribute("DZ");
        registerRealAttribute("DTHETA");
        registerRealAttribute("DPHI");
        registerRealAttribute("DPSI");
        first = false;
    }

}


OpalElement::OpalElement(const std::string &name, OpalElement *parent):
    Element(name, parent), itsSize(parent->itsSize)
{}


OpalElement::~OpalElement()
{}


void OpalElement::
fillRegisteredAttributes(const ElementBase &base, ValueFlag) {
    // Fill in the common data for all elements.
    attributeRegistry["NAME"]->setString(getOpalName());
    attributeRegistry["TYPE"]->setString(getTypeName());
    attributeRegistry["CLASS"]->setString(getParent()->getOpalName());
    attributeRegistry["KEYWORD"]->setString(getBaseObject()->getOpalName());
    attributeRegistry["L"]->setReal(base.getElementLength());

    CoordinateSystemTrafo global2local = base.getCSTrafoGlobal2Local();
    Vector_t origin = global2local.getOrigin();
    Vector_t orientation = Util::getTaitBryantAngles(global2local.getRotation().conjugate());
    attributeRegistry["X"]->setReal(origin[0]);
    attributeRegistry["Y"]->setReal(origin[1]);
    attributeRegistry["Z"]->setReal(origin[2]);
    attributeRegistry["THETA"]->setReal(orientation[0]);
    attributeRegistry["PHI"]->setReal(orientation[1]);
    attributeRegistry["PSI"]->setReal(orientation[2]);

    // Misalignments.
    const AlignWrapper *wrap = dynamic_cast<const AlignWrapper *>(&base);
    if(wrap) {
        double dx, dy, dz, dphi, dtheta, dpsi;
        wrap->offset().getAll(dx, dy, dz, dtheta, dphi, dpsi);
        attributeRegistry["DX"]->setReal(dx);
        attributeRegistry["DY"]->setReal(dy);
        attributeRegistry["DZ"]->setReal(dz);
        attributeRegistry["DTHETA"]->setReal(dtheta);
        attributeRegistry["DPHI"]->setReal(dphi);
        attributeRegistry["DPSI"]->setReal(dpsi);
    }

    CoordinateSystemTrafo misalignment = base.getMisalignment();
    Vector_t misalignmentShift = misalignment.getOrigin();
    Vector_t misalignmentAngles = Util::getTaitBryantAngles(misalignment.getRotation().conjugate());

    attributeRegistry["DX"]->setReal(misalignmentShift(0));
    attributeRegistry["DY"]->setReal(misalignmentShift(1));
    attributeRegistry["DZ"]->setReal(misalignmentShift(2));
    attributeRegistry["DTHETA"]->setReal(misalignmentAngles[0]);
    attributeRegistry["DPHI"]->setReal(misalignmentAngles[1]);
    attributeRegistry["DPSI"]->setReal(misalignmentAngles[2]);

    // Fill in the "unknown" attributes.
    ElementImage *image = base.ElementBase::getImage();
    AttributeSet::const_iterator cur = image->begin();
    AttributeSet::const_iterator end = image->end();
    for(; cur != end; ++cur) {
        attributeRegistry[cur->first]->setReal(cur->second);
    }
}


AttCell *OpalElement::findRegisteredAttribute(const std::string &name) {
    AttCell *cell = &*attributeRegistry[name];

    if(cell == 0) {
        std::string::size_type i = 0;

        if(name[i] == 'K') {
            ++i;
            while(isdigit(name[i])) ++i;
            if(name[i] == 'S') ++i;

            if(name[i] == 'L'  &&  ++i == name.length()) {
                attributeRegistry[name] = cell = new AttReal();
            } else {
                throw OpalException("OpalElement::findRegisteredAttribute()",
                                    "There is no element which has an attribute "
                                    "called \"" + name + "\".");
            }
        }
    }

    return cell;
}

std::pair<ElementBase::ApertureType, std::vector<double> > OpalElement::getApert() const {

    std::pair<ElementBase::ApertureType, std::vector<double> > retvalue(ElementBase::ELLIPTICAL,
                                                                        std::vector<double>({0.5, 0.5, 1.0}));
    if (!itsAttr[APERT]) return retvalue;

    std::string aperture = Attributes::getString(itsAttr[APERT]);

    boost::regex square("square *\\((.*)\\)", boost::regex::icase);
    boost::regex rectangle("rectangle *\\((.*)\\)", boost::regex::icase);
    boost::regex circle("circle *\\((.*)\\)", boost::regex::icase);
    boost::regex ellipse("ellipse *\\((.*)\\)", boost::regex::icase);

    boost::regex twoArguments("([^,]*),([^,]*)");
    boost::regex threeArguments("([^,]*),([^,]*),([^,]*)");

    boost::smatch match;

    const double width2HalfWidth = 0.5;

    if (boost::regex_search(aperture, match, square)) {
        std::string arguments = match[1];
        if (!boost::regex_search(arguments, match, twoArguments)) {
            retvalue.first = ElementBase::RECTANGULAR;

            try {
                retvalue.second[0] = width2HalfWidth * std::stod(arguments);
                retvalue.second[1] = retvalue.second[0];
            } catch (const std::exception &ex) {
                throw OpalException("OpalElement::getApert()",
                                    "could not convert '" + arguments + "' to double");
            }

        } else {
            retvalue.first = ElementBase::CONIC_RECTANGULAR;

            try {
                retvalue.second[0] = width2HalfWidth * std::stod(match[1]);
                retvalue.second[1] = retvalue.second[0];
                retvalue.second[2] = std::stod(match[2]);
            } catch (const std::exception &ex) {
                throw OpalException("OpalElement::getApert()",
                                    "could not convert '" + arguments + "' to doubles");
            }
        }

        return retvalue;
    }

    if (boost::regex_search(aperture, match, rectangle)) {
        std::string arguments = match[1];

        if (!boost::regex_search(arguments, match, threeArguments)) {
            retvalue.first = ElementBase::RECTANGULAR;

            try {
                size_t sz = 0;

                retvalue.second[0] = width2HalfWidth * std::stod(arguments, &sz);
                sz = arguments.find_first_of(",", sz) + 1;
                retvalue.second[1] = width2HalfWidth * std::stod(arguments.substr(sz));

            } catch (const std::exception &ex) {
                throw OpalException("OpalElement::getApert()",
                                    "could not convert '" + arguments + "' to doubles");
            }

        } else {
            retvalue.first = ElementBase::CONIC_RECTANGULAR;

            try {
                retvalue.second[0] = width2HalfWidth * std::stod(match[1]);
                retvalue.second[1] = width2HalfWidth * std::stod(match[2]);
                retvalue.second[2] = std::stod(match[3]);
            } catch (const std::exception &ex) {
                throw OpalException("OpalElement::getApert()",
                                    "could not convert '" + arguments + "' to doubles");
            }
        }

        return retvalue;
    }

    if (boost::regex_search(aperture, match, circle)) {
        std::string arguments = match[1];
        if (!boost::regex_search(arguments, match, twoArguments)) {
            retvalue.first = ElementBase::ELLIPTICAL;

            try {
                retvalue.second[0] = width2HalfWidth * std::stod(arguments);
                retvalue.second[1] = retvalue.second[0];
            } catch (const std::exception &ex) {
                throw OpalException("OpalElement::getApert()",
                                    "could not convert '" + arguments + "' to double");
            }

        } else {
            retvalue.first = ElementBase::CONIC_ELLIPTICAL;

            try {
                retvalue.second[0] = width2HalfWidth * std::stod(match[1]);
                retvalue.second[1] = retvalue.second[0];
                retvalue.second[2] = std::stod(match[2]);
            } catch (const std::exception &ex) {
                throw OpalException("OpalElement::getApert()",
                                    "could not convert '" + arguments + "' to doubles");
            }
        }

        return retvalue;
    }

    if (boost::regex_search(aperture, match, ellipse)) {
        std::string arguments = match[1];

        if (!boost::regex_search(arguments, match, threeArguments)) {
            retvalue.first = ElementBase::ELLIPTICAL;

            try {
                size_t sz = 0;

                retvalue.second[0] = width2HalfWidth * std::stod(arguments, &sz);
                sz = arguments.find_first_of(",", sz) + 1;
                retvalue.second[1] = width2HalfWidth * std::stod(arguments.substr(sz));

            } catch (const std::exception &ex) {
                throw OpalException("OpalElement::getApert()",
                                    "could not convert '" + arguments + "' to doubles");
            }

        } else {
            retvalue.first = ElementBase::CONIC_ELLIPTICAL;

            try {
                retvalue.second[0] = width2HalfWidth * std::stod(match[1]);
                retvalue.second[1] = width2HalfWidth * std::stod(match[2]);
                retvalue.second[2] = std::stod(match[3]);
            } catch (const std::exception &ex) {
                throw OpalException("OpalElement::getApert()",
                                    "could not convert '" + arguments + "' to doubles");
            }
        }

        return retvalue;
    }

    if (aperture != "")
        throw OpalException("OpalElement::getApert()",
                            "Unknown aperture type '" + aperture + "'.");

    return retvalue;
}

double OpalElement::getLength() const {
    return Attributes::getReal(itsAttr[LENGTH]);
}


const std::string OpalElement::getTypeName() const {
    const Attribute *attr = findAttribute("TYPE");
    return attr ? Attributes::getString(*attr) : std::string();
}

/**
   Functions to get the wake field parametes
*/

const std::string OpalElement::getWakeF() const {
    const Attribute *attr = findAttribute("WAKEF");
    return attr ? Attributes::getString(*attr) : std::string();
}

const std::string OpalElement::getParticleMatterInteraction() const {
    const Attribute *attr = findAttribute("PARTICLEMATTERINTERACTION");
    return attr ? Attributes::getString(*attr) : std::string();
}

void OpalElement::parse(Statement &stat) {
    while(stat.delimiter(',')) {
        std::string name = Expressions::parseString(stat, "Attribute name expected.");
        Attribute *attr = findAttribute(name);

        if(attr == 0) {
            throw OpalException("OpalElement::parse",
                                "unknown attribute \"" + name + "\"");
        }

        if(stat.delimiter('[')) {
            int index = int(Round(Expressions::parseRealConst(stat)));
            Expressions::parseDelimiter(stat, ']');

            if(stat.delimiter('=')) {
                attr->parseComponent(stat, true, index);
            } else if(stat.delimiter(":=")) {
                attr->parseComponent(stat, false, index);
            } else {
                throw ParseError("OpalElement::parse()",
                                 "Delimiter \"=\" or \":=\" expected.");
            }
        } else {
            if(stat.delimiter('=')) {
                attr->parse(stat, true);
            } else if(stat.delimiter(":=")) {
                attr->parse(stat, false);
            } else {
                attr->setDefault();
            }
        }
    }
}


void OpalElement::print(std::ostream &os) const {
    std::string head = getOpalName();

    Object *parent = getParent();
    if(parent != 0  &&  ! parent->getOpalName().empty()) {
        if(! getOpalName().empty()) head += ':';
        head += parent->getOpalName();
    }

    os << head;
    os << ';'; // << "JMJdebug OPALElement.cc" ;
    os << std::endl;
}


void OpalElement::setRegisteredAttribute
(const std::string &name, double value) {
    attributeRegistry[name]->setReal(value);
}


void OpalElement::setRegisteredAttribute
(const std::string &name, const std::string &value) {
    attributeRegistry[name]->setString(value);
}


void OpalElement::printMultipoleStrength
(std::ostream &os, int order, int &len,
 const std::string &sName, const std::string &tName,
 const Attribute &length, const Attribute &sNorm, const Attribute &sSkew) {
    // Find out which type of output is required.
    int flag = 0;
    if(sNorm) {
        if(sNorm.getBase().isExpression()) {
            flag += 2;
        } else if(Attributes::getReal(sNorm) != 0.0) {
            flag += 1;
        }
    }

    if(sSkew) {
        if(sSkew.getBase().isExpression()) {
            flag += 6;
        } else if(Attributes::getReal(sSkew) != 0.0) {
            flag += 3;
        }
    }
    //  cout << "JMJdebug, OpalElement.cc: flag=" << flag << endl ;
    // Now do the output.
    int div = 2 * (order + 1);

    switch(flag) {

        case 0:
            // No component at all.
            break;

        case 1:
        case 2:
            // Pure normal component.
        {
            std::string normImage = sNorm.getImage();
            if(length) {
                normImage = "(" + normImage + ")*(" + length.getImage() + ")";
            }
            printAttribute(os, sName, normImage, len);
        }
        break;

        case 3:
        case 6:
            // Pure skew component.
        {
            std::string skewImage = sSkew.getImage();
            if(length) {
                skewImage = "(" + skewImage + ")*(" + length.getImage() + ")";
            }
            printAttribute(os, sName, skewImage, len);
            double tilt = Physics::pi / double(div);
            printAttribute(os, tName, tilt, len);
        }
        break;

        case 4:
            // Both components are non-zero constants.
        {
            double sn = Attributes::getReal(sNorm);
            double ss = Attributes::getReal(sSkew);
            double strength = sqrt(sn * sn + ss * ss);
            if(strength) {
                std::ostringstream ts;
                ts << strength;
                std::string image = ts.str();
                if(length) {
                    image = "(" + image + ")*(" + length.getImage() + ")";
                }
                printAttribute(os, sName, image, len);
                double tilt = - atan2(ss, sn) / double(div);
                if(tilt) printAttribute(os, tName, tilt, len);
            }
        }
        break;

        case 5:
        case 7:
        case 8:
            // One or both components is/are expressions.
        {
            std::string normImage = sNorm.getImage();
            std::string skewImage = sSkew.getImage();
            std::string image =
                "SQRT((" + normImage + ")^2+(" + skewImage + ")^2)";
            printAttribute(os, sName, image, len);
            if(length) {
                image = "(" + image + ")*(" + length.getImage() + ")";
            }
            std::string divisor;
            if(div < 9) {
                divisor = "0";
                divisor[0] += div;
            } else {
                divisor = "00";
                divisor[0] += div / 10;
                divisor[1] += div % 10;
            }
            image = "-ATAN2(" + skewImage + ',' + normImage + ")/" + divisor;
            printAttribute(os, tName, image, len);
            break;
        }
    }
}

void OpalElement::update() {
    ElementBase *base = getElement()->removeWrappers();

    auto apert = getApert();
    base->setAperture(apert.first, apert.second);

    if (itsAttr[ORIGIN] || itsAttr[ORIENTATION]) {
        std::vector<double> ori = Attributes::getRealArray(itsAttr[ORIGIN]);
        std::vector<double> dir = Attributes::getRealArray(itsAttr[ORIENTATION]);
        Vector_t origin(0.0);
        Quaternion rotation;

        if (dir.size() == 3) {
            Quaternion rotTheta(cos(0.5 * dir[0]), 0,                 sin(0.5 * dir[0]), 0);
            Quaternion rotPhi(cos(0.5 * dir[1]),   sin(0.5 * dir[1]), 0,                 0);
            Quaternion rotPsi(cos(0.5 * dir[2]),   0,                 0,                 sin(0.5 * dir[2]));
            rotation = rotTheta * (rotPhi * rotPsi);
        } else {
            if (itsAttr[ORIENTATION]) {
                throw OpalException("Line::parse","Parameter orientation is array of 3 values (theta, phi, psi);\n" +
                                    std::to_string(dir.size()) + " values provided");
            }
        }

        if (ori.size() == 3) {
            origin = Vector_t(ori[0], ori[1], ori[2]);
        } else {
            if (itsAttr[ORIGIN]) {
                throw OpalException("Line::parse","Parameter origin is array of 3 values (x, y, z);\n" +
                                    std::to_string(ori.size()) + " values provided");
            }
        }

        CoordinateSystemTrafo global2local(origin,
                                           rotation.conjugate());
        base->setCSTrafoGlobal2Local(global2local);
        base->fixPosition();

    } else if (!itsAttr[PSI].defaultUsed() &&
               itsAttr[X].defaultUsed() &&
               itsAttr[Y].defaultUsed() &&
               itsAttr[Z].defaultUsed() &&
               itsAttr[THETA].defaultUsed() &&
               itsAttr[PHI].defaultUsed()) {
        base->setRotationAboutZ(Attributes::getReal(itsAttr[PSI]));
    } else if (!itsAttr[X].defaultUsed() ||
               !itsAttr[Y].defaultUsed() ||
               !itsAttr[Z].defaultUsed() ||
               !itsAttr[THETA].defaultUsed() ||
               !itsAttr[PHI].defaultUsed() ||
               !itsAttr[PSI].defaultUsed()) {
        const Vector_t origin(Attributes::getReal(itsAttr[X]),
                              Attributes::getReal(itsAttr[Y]),
                              Attributes::getReal(itsAttr[Z]));

        const double theta = Attributes::getReal(itsAttr[THETA]);
        const double phi = Attributes::getReal(itsAttr[PHI]);
        const double psi = Attributes::getReal(itsAttr[PSI]);

        Quaternion rotTheta(cos(0.5 * theta), 0,              sin(0.5 * theta), 0);
        Quaternion rotPhi(cos(0.5 * phi),     sin(0.5 * phi), 0,                0);
        Quaternion rotPsi(cos(0.5 * psi),     0,              0,                sin(0.5 * psi));
        Quaternion rotation = rotTheta * (rotPhi * rotPsi);

        CoordinateSystemTrafo global2local(origin,
                                           rotation.conjugate());
        base->setCSTrafoGlobal2Local(global2local);
        base->fixPosition();
        base->setRotationAboutZ(Attributes::getReal(itsAttr[PSI]));
    }

    Vector_t misalignmentShift(Attributes::getReal(itsAttr[DX]),
                               Attributes::getReal(itsAttr[DY]),
                               Attributes::getReal(itsAttr[DZ]));
    double dtheta = Attributes::getReal(itsAttr[DTHETA]);
    double dphi = Attributes::getReal(itsAttr[DPHI]);
    double dpsi = Attributes::getReal(itsAttr[DPSI]);
    Quaternion rotationY(cos(0.5 * dtheta), 0,               sin(0.5 * dtheta), 0);
    Quaternion rotationX(cos(0.5 * dphi),   sin(0.5 * dphi), 0,                 0);
    Quaternion rotationZ(cos(0.5 * dpsi),   0,               0,                 sin(0.5 * dpsi));
    Quaternion misalignmentRotation = rotationY * rotationX * rotationZ;
    CoordinateSystemTrafo misalignment(misalignmentShift,
                                       misalignmentRotation.conjugate());

    base->setMisalignment(misalignment);

    if (itsAttr[ELEMEDGE])
        base->setElementPosition(Attributes::getReal(itsAttr[ELEMEDGE]));
}

void OpalElement::updateUnknown(ElementBase *base) {
    for(std::vector<Attribute>::size_type i = itsSize;
        i < itsAttr.size(); ++i) {
        Attribute &attr = itsAttr[i];
        base->setAttribute(attr.getName(), Attributes::getReal(attr));

    }
}


void OpalElement::printAttribute
(std::ostream &os, const std::string &name, const std::string &image, int &len) {
    len += name.length() + image.length() + 2;
    if(len > 74) {
        os << ",&\n  ";
        len = name.length() + image.length() + 3;
    } else {
        os << ',';
    }
    os << name << '=' << image;
}

void OpalElement::printAttribute
(std::ostream &os, const std::string &name, double value, int &len) {
    std::ostringstream ss;
    ss << value << std::ends;
    printAttribute(os, name, ss.str(), len);
}


AttCell *OpalElement::registerRealAttribute(const std::string &name) {
    OwnPtr<AttCell> &cell = attributeRegistry[name];
    if(! cell.isValid()) {
        cell = new AttReal();
    }
    return &*cell;
}


AttCell *OpalElement::registerStringAttribute(const std::string &name) {
    OwnPtr<AttCell> &cell = attributeRegistry[name];
    if(! cell.isValid()) {
        cell = new AttString();
    }
    return &*cell;
}

void OpalElement::registerOwnership() const {
    if (getParent() != 0) return;

    const unsigned int end = itsSize;
    const std::string name = getOpalName();
    for (unsigned int i = COMMON; i < end; ++ i) {
        AttributeHandler::addAttributeOwner(name, AttributeHandler::ELEMENT, itsAttr[i].getName());
    }
}