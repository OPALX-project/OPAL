
// ------------------------------------------------------------------------
// $RCSfile: Line.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.3.4.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: Line
//   The class for OPAL beam lines.
//
// ------------------------------------------------------------------------
//
// $Date: 2002/12/09 15:06:08 $
// $Author: jsberg $
//
// ------------------------------------------------------------------------

#include "Lines/Line.h"
#include "AbsBeamline/ElementBase.h"
#include "AbstractObjects/BeamSequence.h"
#include "AbstractObjects/Expressions.h"
#include "AbstractObjects/OpalData.h"
#include "Attributes/Attributes.h"
#include "Beamlines/Beamline.h"
#include "Beamlines/FlaggedElmPtr.h"
#include "Errors/AlignHandler.h"
#include "Errors/AlignWriter.h"
#include "Expressions/SAutomatic.h"
#include "Expressions/SBinary.h"
#include "Expressions/SRefExpr.h"
#include "Expressions/TFunction2.h"
#include "Lines/LineTemplate.h"
#include "Lines/Replacer.h"
#include "Parser/Statement.h"
#include "Utilities/Options.h"
#include "Utilities/ParseError.h"
#include <algorithm>
#include <iostream>

using namespace Expressions;

// Class Line
// ------------------------------------------------------------------------

// The attributes of class Line.
namespace {
    enum {
        TYPE,        // The type attribute.
        LENGTH,      // The line length.
        ORIGIN,      // The location of the particle source
        ORIENTATION,      // The orientation of the particle source
        X,
        Y,
        Z,
        THETA,
        PHI,
        PSI,
        SIZE
    };

    // Auxiliary function for adding two element lengths.
    double AttAdd(double a, double b)
    { return a + b; }
}


Line::Line():
    BeamSequence(SIZE, "LINE",
                 "The \"LINE\" statement defines a beamline list.\n"
                 "\t<name> : line = (<list>)") {
    itsAttr[TYPE] = Attributes::makeString
                    ("TYPE", "Design type");

    itsAttr[LENGTH] = Attributes::makeReal
                      ("L", "Total length of line in m");
    itsAttr[LENGTH].setReadOnly(true);

    itsAttr[ORIGIN] = Attributes::makeRealArray
        ("ORIGIN", "The location of the particle source");

    itsAttr[ORIENTATION] = Attributes::makeRealArray
        ("ORIENTATION", "The Tait-Bryan angles for the orientation of the particle source");

    itsAttr[X] = Attributes::makeReal
        ("X", "The x-coordinate of the location of the particle source", 0);

    itsAttr[Y] = Attributes::makeReal
        ("Y", "The y-coordinate of the location of the particle source", 0);

    itsAttr[Z] = Attributes::makeReal
        ("Z", "The z-coordinate of the location of the particle source", 0);

    itsAttr[THETA] = Attributes::makeReal
        ("THETA", "The rotation about the y-axis of the particle source", 0);

    itsAttr[PHI] = Attributes::makeReal
        ("PHI", "The rotation about the x-axis of the particle source", 0);

    itsAttr[PSI] = Attributes::makeReal
        ("PSI", "The rotation about the z-axis of the particle source", 0);

    setElement((new FlaggedBeamline("LINE"))->makeAlignWrapper());

    registerOwnership(AttributeHandler::STATEMENT);
}


Line::Line(const std::string &name, Line *parent):
    // Do not copy list members; they will be filled at parse time.
    BeamSequence(name, parent) {
    setElement((new FlaggedBeamline(name))->makeAlignWrapper());
}


Line::~Line()
{}


Line *Line::clone(const std::string &name) {
    return new Line(name, this);
}


Line *Line::copy(const std::string &name) {
    Line *clone = new Line(name, this);
    FlaggedBeamline *oldLine = fetchLine();
    FlaggedBeamline *newLine = clone->fetchLine();
    std::copy(oldLine->begin(), oldLine->end(), std::back_inserter(*newLine));
    return &*clone;
}

double Line::getLength() const {
    return Attributes::getReal(itsAttr[LENGTH]);
}


Object *Line::makeTemplate
(const std::string &name, TokenStream &is, Statement &statement) {
    LineTemplate *macro = new LineTemplate(name, this);

    try {
        macro->parseTemplate(is, statement);
        return macro;
    } catch(...) {
        delete macro;
        macro = 0;
        throw;
    }
}


void Line::parse(Statement &stat) {
    static const TFunction2<double, double> plus = { "+", 4, AttAdd };

    // Check for delimiters.
    parseDelimiter(stat, '=');
    parseDelimiter(stat, '(');

    // Clear all reference counts.
    OpalData::getInstance()->apply(OpalData::ClearReference());

    // Parse the line list.
    parseList(stat);

    // Insert the begin and end markers.
    FlaggedBeamline *line = fetchLine();
    ElementBase *markS = Element::find("#S")->getElement();
    FlaggedElmPtr first(markS, false);
    line->push_front(first);

    ElementBase *markE = Element::find("#E")->getElement();
    FlaggedElmPtr last(markE, false);
    line->push_back(last);

    // Construct expression for length, and fill in occurrence counters.
    PtrToScalar<double> expr;
    for(FlaggedBeamline::iterator i = line->begin(); i != line->end(); ++i) {
        // Accumulate length.
        PtrToScalar<double> temp =
            new SRefExpr<double>(i->getElement()->getName(), "L");
        if(expr.isValid() &&  temp.isValid()) {
            expr = SBinary<double, double>::make(plus, expr, temp);
        } else {
            expr = temp;
        }

        // Do the required instantiations.
        ElementBase *base = i->getElement();
        i->setElement(base->copyStructure());
        Element *elem = Element::find(base->getName());
        i->setCounter(elem->increment());
    }
    if(expr.isValid()) itsAttr[LENGTH].set(new SAutomatic<double>(expr));

    while(stat.delimiter(',')) {
        std::string name = Expressions::parseString(stat, "Attribute name expected.");
        Attribute *attr = findAttribute(name);

        if (attr != 0) {
            if(stat.delimiter('=')) {
                attr->parse(stat, true);
            } else if(stat.delimiter(":=")) {
                attr->parse(stat, false);
            } else {
                attr->setDefault();
            }
        }
    }

    if (itsAttr[ORIGIN] || itsAttr[ORIENTATION]) {
        std::vector<double> origin = Attributes::getRealArray(itsAttr[ORIGIN]);
        if (origin.size() == 3) {
            line->setOrigin3D(Vector_t(origin[0],
                                       origin[1],
                                       origin[2]));
        } else {
            line->setOrigin3D(Vector_t(0.0));
            if (itsAttr[ORIGIN]) {
                throw OpalException("Line::parse","Parameter origin is array of 3 values (x, y, z);\n" +
                                    std::to_string(origin.size()) + " values provided");
            }
        }

        std::vector<double> direction = Attributes::getRealArray(itsAttr[ORIENTATION]);
        if (direction.size() == 3) {
            const double &theta = direction[0];
            const double &phi = direction[1];
            const double &psi = direction[2];

            Quaternion rotTheta(cos(0.5 * theta), 0, -sin(0.5 * theta), 0);
            Quaternion rotPhi(cos(0.5 * phi), -sin(0.5 * phi), 0, 0);
            Quaternion rotPsi(cos(0.5 * psi), 0, 0, -sin(0.5 * psi));
            line->setInitialDirection(rotPsi * rotPhi * rotTheta);
        } else {
            line->setInitialDirection(Quaternion(1, 0, 0, 0));
            if (itsAttr[ORIENTATION]) {
                throw OpalException("Line::parse","Parameter orientation is array of 3 values (theta, phi, psi);\n" +
                                    std::to_string(direction.size()) + " values provided");
            }
        }

        line->setRelativeFlag(true);
    } else {
        const Vector_t origin(Attributes::getReal(itsAttr[X]),
                              Attributes::getReal(itsAttr[Y]),
                              Attributes::getReal(itsAttr[Z]));

        const double theta = Attributes::getReal(itsAttr[THETA]);
        const double phi = Attributes::getReal(itsAttr[PHI]);
        const double psi = Attributes::getReal(itsAttr[PSI]);

        Quaternion rotTheta(cos(0.5 * theta), 0, -sin(0.5 * theta), 0);
        Quaternion rotPhi(cos(0.5 * phi), -sin(0.5 * phi), 0, 0);
        Quaternion rotPsi(cos(0.5 * psi), 0, 0, -sin(0.5 * psi));

        line->setOrigin3D(origin);
        line->setInitialDirection(rotPsi * rotPhi * rotTheta);

        line->setRelativeFlag(!itsAttr[X].defaultUsed() ||
                              !itsAttr[Y].defaultUsed() ||
                              !itsAttr[Z].defaultUsed() ||
                              !itsAttr[THETA].defaultUsed() ||
                              !itsAttr[PHI].defaultUsed() ||
                              !itsAttr[PSI].defaultUsed());
    }
}


void Line::print(std::ostream &os) const {
    os << getOpalName() << ":LINE=(";
    const FlaggedBeamline *line = fetchLine();
    bool seen = false;

    for(FlaggedBeamline::const_iterator i = line->begin();
        i != line->end(); ++i) {
        const std::string name = i->getElement()->getName();
        if(name[0] != '#') {
            if(seen) os << ',';
            if(i->getReflectionFlag()) os << '-';
            os << name;
            seen = true;
        }
    }

    os << ')';
    os << ';';
    os << std::endl;
}


FlaggedBeamline *Line::fetchLine() const {
    return dynamic_cast<FlaggedBeamline *>(getElement()->removeWrappers());
}


void Line::parseList(Statement &stat) {
    FlaggedBeamline *line = fetchLine();

    do {
        // Reversed member ?
        bool rev = stat.delimiter('-');

        // Repetition count.
        int repeat = 1;
        if(stat.integer(repeat)) parseDelimiter(stat, '*');

        // List member.
        if(stat.delimiter('(')) {
            // Anonymous sub-line is expanded immediately.
            Line nestedLine;
            nestedLine.parseList(stat);
            FlaggedBeamline *subLine = nestedLine.fetchLine();

            while(repeat-- > 0) {
                if(rev) {
                    for(FlaggedBeamline::reverse_iterator i = subLine->rbegin();
                        i != subLine->rend(); ++i) {
                        FlaggedElmPtr fep(*i);
                        fep.setReflectionFlag(false);
                        line->push_back(fep);
                    }
                } else {
                    for(FlaggedBeamline::iterator i = subLine->begin();
                        i != subLine->end(); ++i) {
                        FlaggedElmPtr fep(*i);
                        line->push_back(fep);
                    }
                }
            }
        } else {
            // Identifier.
            std::string name = parseString(stat, "Line member expected.");
            Pointer<Object> obj = OpalData::getInstance()->find(name);

            if(! obj.isValid()) {
                throw ParseError("Line::parseList()",
                                 "Element \"" + name + "\" is undefined.");
            }

            if(stat.delimiter('(')) {
                // Line or sequence macro.
                // This instance will always be an anonymous object.
                obj = obj->makeInstance(name, stat, 0);
            }

            if(Element *elem = dynamic_cast<Element *>(&*obj)) {
                while(repeat-- > 0) {
                    ElementBase *base = elem->getElement();
                    FlaggedElmPtr member(base, rev);
                    line->push_back(member);
                }
            } else {
                throw ParseError("Line::parseList()",
                                 "Object \"" + name + "\" cannot be a line member.");
            }
        }
    } while(stat.delimiter(','));

    parseDelimiter(stat, ')');
}


void Line::replace(Object *oldObject, Object *newObject) {
    Element *oldElement = dynamic_cast<Element *>(oldObject);
    Element *newElement = dynamic_cast<Element *>(newObject);

    if(oldElement != 0  &&  newElement != 0) {
        Replacer rep(*fetchLine(),
                     oldElement->getOpalName(),
                     newElement->getElement());
        rep.execute();
    }
}