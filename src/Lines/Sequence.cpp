// ------------------------------------------------------------------------
// $RCSfile: Sequence.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.3.4.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: Sequence
//   The class of OPAL element sequences.
//
// ------------------------------------------------------------------------
//
// $Date: 2002/12/09 15:06:08 $
// $Author: jsberg $
//
// ------------------------------------------------------------------------

#include "Lines/Sequence.h"
#include "AbstractObjects/Attribute.h"
#include "AbstractObjects/OpalData.h"
#include "Attributes/Attributes.h"
#include "BeamlineCore/DriftRep.h"
#include "Beamlines/Beamline.h"
#include "Elements/OpalDrift.h"
#include "Lines/FlatWriter.h"

#include "Lines/Replacer.h"
#include "Lines/SequenceParser.h"
#include "Lines/SequenceTemplate.h"
#include "Utilities/OpalException.h"
#include "Utilities/Options.h"

#include <cmath>
#include <iostream>
#include <sstream>
#include <vector>

// Class Sequence
// ------------------------------------------------------------------------

// The attributes of class Sequence.
namespace {
    enum {
        TYPE,     // The design type.
        LENGTH,   // The total sequence length.
        REFER,    // The reference position for members.
        REFPOS,   // The reference position for the sequence.
        SIZE
    };
}


Sequence::Sequence():
    BeamSequence(SIZE, "SEQUENCE",
                 "The \"SEQUENCE\" statement initiates parsing of an "
                 "element sequence.\n"
                 "\t<label>: SEQUENCE,L=<length>,REFER=<reference>\n"
                 "\t\t...\n"
                 "\t\t<object>: <class>,AT=<real>{,<attribute>=<value>}\n"
                 "\t\t<object>: <class>,DRIFT=<real>{,<attribute>=<value>}\n"
                 "\t\t...\n"
                 "\tEND;") {
    itsAttr[TYPE] = Attributes::makeString
                    ("TYPE", "The design type");
    itsAttr[LENGTH] = Attributes::makeReal
                      ("L", "Total length of sequence in m");

    itsAttr[REFER] = Attributes::makeString
                     ("REFER",
                      "Reference position for members:\n"
                      "\tENTRY | EXIT | CENTRE (default is CENTRE)", "CENTRE");
    itsAttr[REFPOS] = Attributes::makeString
                      ("REFPOS",
                       "Element giving reference position for this sequence"
                       "\t(if given, this position is used instead of the centre, when the"
                       "\tsequence is nested in another sequence with \"REFER=CENTRE\")");

    setElement((new TLine("SEQUENCE"))->makeAlignWrapper());

    registerOwnership(AttributeHandler::STATEMENT);
}


Sequence::Sequence(const std::string &name, Sequence *parent):
    BeamSequence(name, parent)
    // Do not copy parent's line, it will be filled in at parse time,
    // In case of a clone within a sequence, it is filled by the method
    // SequenceParser::parseMember().
{
    setElement((new TLine(name))->makeAlignWrapper());
}


Sequence::~Sequence()
{}


Sequence *Sequence::clone(const std::string &name) {
    return new Sequence(name, this);
}


Sequence *Sequence::copy(const std::string &name) {
    TLine *oldLine = fetchLine();
    TLine *newLine = new TLine(name);

    for(TLine::iterator i = oldLine->begin(); i != oldLine->end(); ++i) {
        SequenceMember member(*i);

        if(i->itsType == SequenceMember::GENERATED) {
            member.setElement(member.getElement()->clone());
        } else {
            member.setElement(member.getElement()->copyStructure());
        }

        newLine->push_back(member);
    }

    Sequence *newSeq = dynamic_cast<Sequence *>(clone(name));
    newSeq->storeLine(*newLine);
    return newSeq;
}

Sequence::TLine::iterator
Sequence::findNamedPosition(TLine &line, const std::string &name) const {
    TLine::iterator first = line.begin();
    TLine::iterator last  = line.end();

    for(TLine::iterator iter = first; iter != last; ++iter) {
        if(iter->itsType != SequenceMember::GENERATED) {
            if(iter->getElement()->getName() == name) {
                TLine::iterator test = iter;

                while(++test != line.end()) {
                    if(test->getElement()->getName() == name) {
                        throw OpalException("Sequence::findNamedPosition()", "Element \"" +
                                            name + "\" is not unique in sequence.");
                    }
                }

                return iter;
            }
        }
    }

    throw OpalException("Sequence::findNamedPosition()",
                        "Element \"" + name + "\" not found in sequence.");
}


double Sequence::getLength() const {
    return Attributes::getReal(itsAttr[LENGTH]);
}


double Sequence::getEntrance(ReferenceType ref) const {
    if(itsRefPoint.empty()) {
        return Element::getEntrance(ref);
    } else {
        TLine *line = fetchLine();
        TLine::iterator iter = findNamedPosition(*line, itsRefPoint);
        return (- iter->itsPosition);
    }
}


double Sequence::getExit(ReferenceType ref) const {
    if(itsRefPoint.empty()) {
        return Element::getExit(ref);
    } else {
        TLine *line = fetchLine();
        TLine::iterator iter = findNamedPosition(*line, itsRefPoint);
        return (getLength() - iter->itsPosition);
    }
}


Sequence::ReferenceType Sequence::getReference() const {
    std::string ref = Attributes::getString(itsAttr[REFER]);
    ReferenceType code = IS_CENTRE;

    if(ref == "ENTRY") {
        code = IS_ENTRY;
    } else if(ref == "EXIT") {
        code = IS_EXIT;
    }

    return code;
}


Object *Sequence::makeTemplate
(const std::string &name, TokenStream &is, Statement &statement) {
    SequenceTemplate *macro = new SequenceTemplate(name, this);

    try {
        macro->parseTemplate(is, statement);
        return macro;
    } catch(...) {
        delete macro;
        macro = 0;
        throw;
    }
}


void Sequence::parse(Statement &statement) {
    // Look for "REFER" and "L" attributes.
    Object::parse(statement);

    // Set up to parse members.
    SequenceParser parser(this);
    parser.run();
    if(itsAttr[REFPOS]) itsRefPoint = Attributes::getString(itsAttr[REFPOS]);
}



void Sequence::print(std::ostream &os) const {
    TLine *line = fetchLine();
    std::streamsize old_prec = os.precision(12);

    if(isShared()) os << "SHARED ";
    os << getOpalName() << ":SEQUENCE";

    for(const Attribute& i : itsAttr) {
        if(i) {
            os << ',' << i.getName()
               << (i.isExpression() ? ":=" : "=") << i;
        }
    }

    os << ";" << std::endl;

    for(TLine::const_iterator iter = line->begin();
        iter != line->end(); ++iter) {
        const SequenceMember &member = *iter;

        if(member.itsType != SequenceMember::GENERATED) {
            const std::string &name = member.getElement()->getName();
            if(name[0] != '#') {
                // Name of current element.
                os << "  ";
                if(member.getReflectionFlag()) os << '-';
                os << name;

                // Position of current element.
                os << ",AT=" << member.itsPosition;
                os << ';' << std::endl;
            }
        }
    }

    os << "ENDSEQUENCE;" << std::endl;

    os.precision(old_prec);
}


void Sequence::replace(Object *oldObject, Object *newObject) {
    Element *oldElement = dynamic_cast<Element *>(oldObject);
    Element *newElement = dynamic_cast<Element *>(newObject);

    if(oldElement != 0  &&  newElement != 0) {
        Replacer rep(*fetchLine(),
                     oldElement->getOpalName(),
                     newElement->getElement());
        rep.execute();
    }
}


void Sequence::update() {
    TLine *line = fetchLine();
    if(! line->empty()) updateList(this, line);
}


Sequence::TLine *Sequence::fetchLine() const {
    return dynamic_cast<TLine *>(getElement()->removeWrappers());
}


void Sequence::storeLine(TLine &newLine) {
    // Remove any old line and assign new one.
    TLine *line = fetchLine();
    line->erase(line->begin(), line->end());
    line->swap(newLine);

    // Store the reference code.
    itsCode = getReference();

    // Clear all occurrence counts.
    OpalData::getInstance()->apply(OpalData::ClearReference());

    // Set occurrence counts.
    for(TLine::iterator i = line->begin(); i != line->end(); ++i) {
        if(i->itsType != SequenceMember::GENERATED) {
            const std::string &name = i->getElement()->getName();
            Element *elem = Element::find(name);
            // ada 4.5 2000 to speed up matching, add a pointer to
            // opal elements in order to avoid serching the opal elements
            i->OpalElement = elem;
            i->setCounter(elem->increment());
        }
    }

    // Fill in the drift spaces.
    // update();
}


void Sequence::addEndMarkers(TLine &line) const {
    SequenceMember member;
    member.setElement(Element::find("#S")->getElement());
    member.itsPosition = 0.0;
    member.itsFlag = SequenceMember::ABSOLUTE;
    member.itsType = SequenceMember::GLOBAL;
    line.push_front(member);
    member.setElement(Element::find("#E")->getElement());
    member.itsPosition = getLength();
    line.push_back(member);
}


double Sequence::findDriftLength(TLine::iterator drift) const {
    // It is assumed that the elements preceding and following the drift exist.
    TLine::iterator prev(drift);
    --prev;
    const std::string prev_name = prev->getElement()->getName();

    // ada 4.5 2000 to speed up matching, add a pointer to
    // opal elements in order to avoid serching the opal elements
    // Element *prev_elem = Element::find(prev_name);
    Pointer<Element> prev_elem = prev->OpalElement;
    double prev_exit = prev->itsPosition + prev_elem->getExit(itsCode);

    TLine::iterator next(drift);
    ++next;
    const std::string next_name = next->getElement()->getName();

    // ada 4.5 2000 to speed up matching, add a pointer to
    // opal elements in order to avoid serching the opal elements
    // Element *next_elem = Element::find(next_name);
    Pointer<Element> next_elem = next->OpalElement;
    double next_entry = next->itsPosition + next_elem->getEntrance(itsCode);

    double driftLength = next_entry - prev_exit;

    static double tolerance = 1.0e-8;
    if(std::abs(driftLength) < tolerance) {
        driftLength = 0.0;
    } else if(driftLength < 0.0) {
        driftLength = 0.0;
        std::ostringstream os;
        os << "Inconsistent positions in sequence \"" << getOpalName() + "\";\n"
           << "previous element:  \"" << prev_name + "\" at = "
           << prev->itsPosition << ",\n"
           << "following element: \"" << next_name + "\" at = "
           << next->itsPosition << "." << std::ends;
        throw OpalException("Sequence::findDriftLength()", os.str());
        // ada 15-6-2000 move driftLength = 0.0; bevore throw
    }

    return driftLength;
}


void Sequence::insertDrifts(Sequence::TLine &line) {
    TLine::iterator i = line.begin();
    while(true) {
        ++i;
        if(i == line.end()) break;
        SequenceMember member;
        DriftRep *drift = new DriftRep();
        drift->setName("[DRIFT]");
        member.setElement(drift);
        member.itsType = SequenceMember::GENERATED;
        line.insert(i, member);
    }
}


void Sequence::updateList(Sequence *seq, TLine *line) {
    TLine::iterator iter = line->begin();
    TLine::iterator last = line->end();

    while(true) {
        // Recursive call for nested beam non-shared sequence.
        if(iter == last) break;
        ElementBase *base = iter->getElement()->removeWrappers();
        if(! base->isSharable()) {
            TLine *sub_line = dynamic_cast<TLine *>(base);
            if(sub_line != 0) {
                const std::string &sub_name = sub_line->getName();
                Sequence *sub_seq = dynamic_cast<Sequence *>(OpalData::getInstance()->find(sub_name));
                updateList(sub_seq, sub_line);
            }
        }
        ++iter;

        // Fill in drift length.
        if(iter == last) break;
        double driftLength = seq->findDriftLength(iter);
        iter->setLength(driftLength);
        ++iter;
    }
}