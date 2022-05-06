//
// Class OpalElement
//   Base class for all beam line elements.
//
// Copyright (c) 200x - 2020, Paul Scherrer Institut, Villigen PSI, Switzerland
// All rights reserved
//
// This file is part of OPAL.
//
// OPAL is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// You should have received a copy of the GNU General Public License
// along with OPAL. If not, see <https://www.gnu.org/licenses/>.
//
#include "Elements/OpalElement.h"

#include "AbsBeamline/Bend2D.h"
#include "AbstractObjects/Attribute.h"
#include "AbstractObjects/Expressions.h"
#include "AbstractObjects/OpalData.h"
#include "Attributes/Attributes.h"
#include "Parser/Statement.h"
#include "Physics/Physics.h"
#include "Utilities/OpalException.h"
#include "Utilities/Options.h"
#include "Utilities/ParseError.h"
#include "Utilities/Util.h"

#include <boost/regex.hpp>

#include <cmath>
#include <cctype>
#include <sstream>
#include <vector>


OpalElement::OpalElement(int size, const char* name, const char* help):
    Element(size, name, help), itsSize(size) {
    itsAttr[TYPE]   = Attributes::makePredefinedString
        ("TYPE", "The element design type.",
         {"RING",
          "CARBONCYCL",
          "CYCIAE",
          "AVFEQ",
          "FFA",
          "BANDRF",
          "SYNCHROCYCLOTRON",
          "SINGLEGAP",
          "STANDING",
          "TEMPORAL",
          "SPATIAL"});

    itsAttr[LENGTH] = Attributes::makeReal
        ("L", "The element length in m");

    itsAttr[ELEMEDGE] = Attributes::makeReal
        ("ELEMEDGE", "The position of the element in path length (in m)");

    itsAttr[APERT] = Attributes::makeString
        ("APERTURE", "The element aperture");

    itsAttr[WAKEF] = Attributes::makeString
        ("WAKEF", "Defines the wake function");

    itsAttr[PARTICLEMATTERINTERACTION] = Attributes::makeString
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

    itsAttr[OUTFN] = Attributes::makeString
        ("OUTFN", "Output filename");

    itsAttr[DELETEONTRANSVERSEEXIT] = Attributes::makeBool
        ("DELETEONTRANSVERSEEXIT", "Flag controlling if particles should be deleted if they exit "
                                   "the element transversally", true);

    const unsigned int end = COMMON;
    for (unsigned int i = 0; i < end; ++ i) {
        AttributeHandler::addAttributeOwner("Any", AttributeHandler::ELEMENT, itsAttr[i].getName());
    }
}


OpalElement::OpalElement(const std::string& name, OpalElement* parent):
    Element(name, parent), itsSize(parent->itsSize)
{}


OpalElement::~OpalElement()
{}


std::pair<ApertureType, std::vector<double> > OpalElement::getApert() const {

    std::pair<ApertureType, std::vector<double> > retvalue(ApertureType::ELLIPTICAL,
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
            retvalue.first = ApertureType::RECTANGULAR;

            try {
                retvalue.second[0] = width2HalfWidth * std::stod(arguments);
                retvalue.second[1] = retvalue.second[0];
            } catch (const std::exception &ex) {
                throw OpalException("OpalElement::getApert()",
                                    "could not convert '" + arguments + "' to double");
            }

        } else {
            retvalue.first = ApertureType::CONIC_RECTANGULAR;

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
            retvalue.first = ApertureType::RECTANGULAR;

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
            retvalue.first = ApertureType::CONIC_RECTANGULAR;

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
            retvalue.first = ApertureType::ELLIPTICAL;

            try {
                retvalue.second[0] = width2HalfWidth * std::stod(arguments);
                retvalue.second[1] = retvalue.second[0];
            } catch (const std::exception &ex) {
                throw OpalException("OpalElement::getApert()",
                                    "could not convert '" + arguments + "' to double");
            }

        } else {
            retvalue.first = ApertureType::CONIC_ELLIPTICAL;

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
            retvalue.first = ApertureType::ELLIPTICAL;

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
            retvalue.first = ApertureType::CONIC_ELLIPTICAL;

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

    if (!aperture.empty()) {
        throw OpalException("OpalElement::getApert()",
                            "Unknown aperture type '" + aperture + "'.");
    }

    return retvalue;
}

double OpalElement::getLength() const {
    return Attributes::getReal(itsAttr[LENGTH]);
}


const std::string OpalElement::getTypeName() const {
    const Attribute* attr = findAttribute("TYPE");
    return attr ? Attributes::getString(*attr) : std::string();
}

/**
   Functions to get the wake field parametes
*/
const std::string OpalElement::getWakeF() const {
    const Attribute* attr = findAttribute("WAKEF");
    return attr ? Attributes::getString(*attr) : std::string();
}

const std::string OpalElement::getParticleMatterInteraction() const {
    const Attribute* attr = findAttribute("PARTICLEMATTERINTERACTION");
    return attr ? Attributes::getString(*attr) : std::string();
}

void OpalElement::parse(Statement& stat) {
    while (stat.delimiter(',')) {
        std::string name = Expressions::parseString(stat, "Attribute name expected.");
        Attribute *attr = findAttribute(name);

        if (attr == 0) {
            throw OpalException("OpalElement::parse",
                                "unknown attribute \"" + name + "\"");
        }

        if (stat.delimiter('[')) {
            int index = int(std::round(Expressions::parseRealConst(stat)));
            Expressions::parseDelimiter(stat, ']');

            if (stat.delimiter('=')) {
                attr->parseComponent(stat, true, index);
            } else if (stat.delimiter(":=")) {
                attr->parseComponent(stat, false, index);
            } else {
                throw ParseError("OpalElement::parse()",
                                 "Delimiter \"=\" or \":=\" expected.");
            }
        } else {
            if (stat.delimiter('=')) {
                attr->parse(stat, true);
            } else if (stat.delimiter(":=")) {
                attr->parse(stat, false);
            } else {
                attr->setDefault();
            }
        }
    }
}


void OpalElement::print(std::ostream& os) const {
    std::string head = getOpalName();

    Object* parent = getParent();
    if (parent != 0  &&  ! parent->getOpalName().empty()) {
        if (! getOpalName().empty()) head += ':';
        head += parent->getOpalName();
    }

    os << head;
    os << ';'; // << "JMJdebug OPALElement.cc" ;
    os << std::endl;
}


void OpalElement::printMultipoleStrength(std::ostream& os,
                                         int order, int& len,
                                         const std::string& sName,
                                         const std::string& tName,
                                         const Attribute& length,
                                         const Attribute& sNorm,
                                         const Attribute& sSkew) {
    // Find out which type of output is required.
    int flag = 0;
    if (sNorm) {
        if (sNorm.getBase().isExpression()) {
            flag += 2;
        } else if (Attributes::getReal(sNorm) != 0.0) {
            flag += 1;
        }
    }

    if (sSkew) {
        if (sSkew.getBase().isExpression()) {
            flag += 6;
        } else if (Attributes::getReal(sSkew) != 0.0) {
            flag += 3;
        }
    }
    //  cout << "JMJdebug, OpalElement.cc: flag=" << flag << endl ;
    // Now do the output.
    int div = 2 * (order + 1);

    switch (flag) {

        case 0:
            // No component at all.
            break;

        case 1:
        case 2:
            // Pure normal component.
        {
            std::string normImage = sNorm.getImage();
            if (length) {
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
            if (length) {
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
            double strength = std::sqrt(sn * sn + ss * ss);
            if (strength) {
                std::ostringstream ts;
                ts << strength;
                std::string image = ts.str();
                if (length) {
                    image = "(" + image + ")*(" + length.getImage() + ")";
                }
                printAttribute(os, sName, image, len);
                double tilt = - std::atan2(ss, sn) / double(div);
                if (tilt) printAttribute(os, tName, tilt, len);
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
            if (length) {
                image = "(" + image + ")*(" + length.getImage() + ")";
            }
            std::string divisor;
            if (div < 9) {
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
    ElementBase* base = getElement();

    auto apert = getApert();
    base->setAperture(apert.first, apert.second);

    if (itsAttr[ORIGIN] || itsAttr[ORIENTATION]) {
        std::vector<double> ori = Attributes::getRealArray(itsAttr[ORIGIN]);
        std::vector<double> dir = Attributes::getRealArray(itsAttr[ORIENTATION]);
        Vector_t origin(0.0);
        Quaternion rotation;

        if (dir.size() == 3) {
            Quaternion rotTheta(std::cos(0.5 * dir[0]), 0,                 std::sin(0.5 * dir[0]), 0);
            Quaternion rotPhi(std::cos(0.5 * dir[1]),   std::sin(0.5 * dir[1]), 0,                 0);
            Quaternion rotPsi(std::cos(0.5 * dir[2]),   0,                 0,                 std::sin(0.5 * dir[2]));
            rotation = rotTheta * (rotPhi * rotPsi);
        } else {
            if (itsAttr[ORIENTATION]) {
                throw OpalException("OpalElement::update",
                                    "Parameter orientation is array of 3 values (theta, phi, psi);\n" +
                                    std::to_string(dir.size()) + " values provided");
            }
        }

        if (ori.size() == 3) {
            origin = Vector_t(ori[0], ori[1], ori[2]);
        } else {
            if (itsAttr[ORIGIN]) {
                throw OpalException("OpalElement::update",
                                    "Parameter origin is array of 3 values (x, y, z);\n" +
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

        Quaternion rotTheta(std::cos(0.5 * theta), 0,              std::sin(0.5 * theta), 0);
        Quaternion rotPhi(std::cos(0.5 * phi),     std::sin(0.5 * phi), 0,                0);
        Quaternion rotPsi(std::cos(0.5 * psi),     0,              0,                std::sin(0.5 * psi));
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
    Quaternion rotationY(std::cos(0.5 * dtheta), 0,               std::sin(0.5 * dtheta), 0);
    Quaternion rotationX(std::cos(0.5 * dphi),   std::sin(0.5 * dphi), 0,                 0);
    Quaternion rotationZ(std::cos(0.5 * dpsi),   0,               0,                 std::sin(0.5 * dpsi));
    Quaternion misalignmentRotation = rotationY * rotationX * rotationZ;
    CoordinateSystemTrafo misalignment(misalignmentShift,
                                       misalignmentRotation.conjugate());

    base->setMisalignment(misalignment);

    if (itsAttr[ELEMEDGE])
        base->setElementPosition(Attributes::getReal(itsAttr[ELEMEDGE]));

    base->setFlagDeleteOnTransverseExit(Attributes::getBool(itsAttr[DELETEONTRANSVERSEEXIT]));
}

void OpalElement::updateUnknown(ElementBase* base) {
    for (std::vector<Attribute>::size_type i = itsSize;
        i < itsAttr.size(); ++i) {
        Attribute &attr = itsAttr[i];
        base->setAttribute(attr.getName(), Attributes::getReal(attr));

    }
}

void OpalElement::printAttribute(std::ostream& os, const std::string& name,
                                 const std::string& image, int& len) {
    len += name.length() + image.length() + 2;
    if (len > 74) {
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


void OpalElement::registerOwnership() const {
    if (getParent() != 0) return;

    const unsigned int end = itsSize;
    const std::string name = getOpalName();
    for (unsigned int i = COMMON; i < end; ++ i) {
        AttributeHandler::addAttributeOwner(name, AttributeHandler::ELEMENT, itsAttr[i].getName());
    }
}
