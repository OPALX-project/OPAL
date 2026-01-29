// ------------------------------------------------------------------------
// $RCSfile: SequenceParser.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.2 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Definitions for class: SequenceParser
//   This class defines the OPAL sequence parser.
//
// ------------------------------------------------------------------------
//
// $Date: 2001/08/13 15:16:16 $
// $Author: jowett $
//
// ------------------------------------------------------------------------

#include "Lines/SequenceParser.h"
#include "AbstractObjects/Attribute.h"
#include "AbstractObjects/BeamSequence.h"
#include "AbstractObjects/Expressions.h"
#include "AbstractObjects/OpalData.h"
#include "Attributes/Attributes.h"
#include "Lines/Sequence.h"
#include "OpalParser/OpalParser.h"
#include "MemoryManagement/Pointer.h"
#include "Parser/Statement.h"
#include "Utilities/ClassicException.h"
#include "Utilities/OpalException.h"
#include "Utilities/Options.h"
#include "Utilities/ParseError.h"


// Class SequenceParser
// ------------------------------------------------------------------------

SequenceParser::SequenceParser(Sequence *seq):
    itsSequence(seq), okFlag(true), itsLine()
{}


SequenceParser::~SequenceParser()
{}


void SequenceParser::parse(Statement &stat) const {
    if(stat.keyword("ENDSEQUENCE")) {
        // Finish setting up of sequence.
        if(okFlag) {
            try {
                // Finish parsing the ENDSEQUENCE statement.
                itsSequence->Object::parse(stat);
                parseEnd(stat);

                // Finish setting up of sequence.
                itsSequence->addEndMarkers(itsLine);
                fillPositions();
                itsSequence->insertDrifts(itsLine);
                itsSequence->storeLine(itsLine);
                OpalData::getInstance()->makeDirty(itsSequence);
            } catch(...) {
                // In case of error, erase sequence and return to previous parser.
                itsLine.clear();
                stop();
                throw;
            }
        }

        // Return to previous parser.
        itsLine.clear();
        stop();
    } else {
        // Backtrack to start and parse new member.
        try {
            stat.start();
            parseMember(stat);
        } catch(...) {
            okFlag = false;
            throw;
        }
    }
}


void SequenceParser::fillPositions() const {
    findFromPositions();
    findNeighbourPositions();

    while(! references.empty()) {
        int count = 0;

        for(RefList::iterator i = references.begin();
            i != references.end();) {
            if(i->fromPosition != 0  &&
               i->fromPosition->itsFlag == SequenceMember::ABSOLUTE) {
                for(std::list<SequenceMember *>::iterator j = i->itsList.begin();
                    j != i->itsList.end(); ++j) {
                    (*j)->itsPosition += i->fromPosition->itsPosition;
                    (*j)->itsFlag = SequenceMember::ABSOLUTE;
                }

                ++count;
                i = references.erase(i);
            } else {
                ++i;
            }
        }

        if(count == 0) {
            throw OpalException("SequenceParser::fillPositions()",
                                "Circular \"FROM\" clauses.");
        }
    }
}


void SequenceParser::findFromPositions() const {
    for(RefList::iterator i = references.begin(); i != references.end(); ++i) {
        TLine::iterator pos = itsSequence->findNamedPosition(itsLine, i->fromName);

        if(pos != itsLine.end()) {
            i->fromPosition = &(*pos);
        } else {
            throw OpalException("SequenceParser::findFromPositions()",
                                "Not all FROM clauses resolved.");
        }
    }
}


void SequenceParser::findNeighbourPositions() const {
    Sequence::ReferenceType type = Sequence::IS_CENTRE;
    Attribute *attr = itsSequence->findAttribute("REFER");
    std::string ref = Attributes::getString(*attr);

    if(ref == "ENTRY") {
        type = Sequence::IS_ENTRY;
    } else if(ref == "CENTRE"  ||  ref == "") {
        type = Sequence::IS_CENTRE;
    } else if(ref == "EXIT") {
        type = Sequence::IS_EXIT;
    } else {
        throw OpalException("SequenceParser::findNeighbourPositions()",
                            "Unknown reference flag \"" + ref +
                            "\" in <sequence> \"" +
                            itsSequence->getOpalName() + "\".");
    }

    for(TLine::iterator i = itsLine.begin(); i != itsLine.end(); ++i) {
        TLine::iterator j = i;
        switch(i->itsFlag) {

            case SequenceMember::DRIFT:
            case SequenceMember::IMMEDIATE:
                // Add scaled length of current.
                if(type != Sequence::IS_ENTRY) {
                    std::string name = i->getElement()->getName();
                    Element *elem = Element::find(name);
                    double length = elem->getLength();

                    if(type == Sequence::IS_CENTRE) {
                        i->itsPosition += length / 2.0;
                    } else if(type == Sequence::IS_EXIT) {
                        i->itsPosition += length;
                    }
                }

                // Add scaled length of previous.
                if(j != itsLine.begin()) {
                    --j;
                    if(type != Sequence::IS_EXIT) {
                        std::string name = j->getElement()->getName();
                        Element *elem = Element::find(name);
                        double length = elem->getLength();

                        if(type == Sequence::IS_ENTRY) {
                            i->itsPosition += length;
                        } else if(type == Sequence::IS_CENTRE) {
                            i->itsPosition += length / 2.0;
                        }
                    }

                    // Current position is relative to previous.
                    Reference ref("DRIFT", &(*i));
                    ref.fromPosition = &(*j);
                    references.push_back(ref);
                } else {
                    // Current position is relative to sequence begin.
                    i->itsFlag = SequenceMember::ABSOLUTE;
                }
                break;

            case SequenceMember::PREVIOUS:
                if(j != itsLine.begin()) {
                    --j;
                    Reference ref("PREVIOUS", &(*i));
                    ref.fromPosition = &(*j);
                    references.push_back(ref);
                } else {
                    i->itsFlag = SequenceMember::ABSOLUTE;
                }
                break;

            case SequenceMember::NEXT:
                ++j;
                if(j != itsLine.end()) {
                    // Current position is relative to *j;
                    Reference ref("NEXT", &(*i));
                    ref.fromPosition = &(*j);
                    references.push_back(ref);
                } else {
                    // Current position is relative to end of sequence.
                    Attribute *attr = itsSequence->findAttribute("L");
                    i->itsPosition += Attributes::getReal(*attr);
                    i->itsFlag = SequenceMember::ABSOLUTE;
                }
                break;

            case SequenceMember::END: {
                // Current position is relative to end of sequence.
                Attribute *attr = itsSequence->findAttribute("L");
                i->itsPosition += Attributes::getReal(*attr);
                i->itsFlag = SequenceMember::ABSOLUTE;
            }
            break;


            default:
                ;
        }
    }
}


void SequenceParser::parseMember(Statement &stat) const {
    std::string objName;
    std::string clsName;
    SequenceMember member;

    // Read object identifier.
    if(stat.delimiter('-')) {
        member.setReflectionFlag(true);
        clsName = Expressions::parseString(stat, "Class name expected.");
        member.itsType = SequenceMember::GLOBAL;
    } else {
        clsName = Expressions::parseString(stat, "Object or class name expected.");

        if(stat.delimiter(':')) {
            member.itsType = SequenceMember::LOCAL;
            member.setReflectionFlag(stat.delimiter('-'));
            objName = clsName;
            clsName = Expressions::parseString(stat, "Class name expected.");
        } else {
            member.itsType = SequenceMember::GLOBAL;
        }
    }

    // Find exemplar object.
    if(Object *obj = OpalData::getInstance()->find(clsName)) {
        Pointer<Object> copy;
        if(obj) {
            if(member.itsType == SequenceMember::LOCAL) {
                // We must make a copy or a macro instance.
                if(objName == itsSequence->getOpalName()) {
                    throw ParseError("SequenceParser::parseMember()",
                                     "You cannot redefine \"" + objName +
                                     "\" within the sequence of the same name.");
                } else if(clsName == itsSequence->getOpalName()) {
                    throw ParseError("SequenceParser::parseMember()",
                                     "You cannot refer to \"" + objName +
                                     "\" within the sequence of the same name.");
                } else if(OpalData::getInstance()->find(objName)) {
                    throw ParseError("SequenceParser::parseMember()",
                                     "You cannot redefine \"" + objName +
                                     "\" within a sequence.");
                }

                if(stat.delimiter('(')) {
                    copy = obj->makeInstance(objName, stat, this);
                } else if(BeamSequence *line = dynamic_cast<BeamSequence *>(obj)) {
                    copy = line->copy(objName);
                } else {
                    copy = obj->clone(objName);
                }
            } else {
                // We can use the global object.
                copy = obj;
            }
        }

        // Test for correct object type.
        if(Element *elem = dynamic_cast<Element *>(&*copy)) {
            ElementBase *base = elem->getElement();
            member.setElement(base->copyStructure());

            // ada 4.5 2000 to speed up matching, add a pointer to
            // opal elements in order to avoid serching the opal elements
            member.OpalElement = elem;

            itsLine.push_back(member);

            // Parse position and attributes.
            parsePosition(stat, *copy, member.itsType == SequenceMember::LOCAL);
            parseEnd(stat);

            // For local instantiation define element.
            if(member.itsType == SequenceMember::LOCAL) OpalData::getInstance()->define(&*copy);
        } else {
            throw ParseError("SequenceParser::parseMember()", "Object \"" +
                             objName + "\" is not an element.");
        }
    } else {
        throw ParseError("SequenceParser::parseMember()",
                         "Element class name \"" + clsName + "\" is unknown.");
    }
}


void SequenceParser::parsePosition
(Statement &stat, Object &elem, bool defined) const {
    SequenceMember &member = itsLine.back();
    double position = 0.0;
    std::string fromName;
    int flag = 0;

    while(stat.delimiter(',')) {
        std::string name =
            Expressions::parseString(stat, "Attribute name expected.");

        if(name == "AT") {
            Expressions::parseDelimiter(stat, '=');
            position = Expressions::parseRealConst(stat);
            flag |= (flag & 0x0001) ? 0x0008 : 0x0001;
        } else if(name == "DRIFT") {
            Expressions::parseDelimiter(stat, '=');
            position = Expressions::parseRealConst(stat);
            flag |= (flag & 0x0002) ? 0x0008 : 0x0002;
        } else if(name == "FROM") {
            Expressions::parseDelimiter(stat, '=');
            if(stat.delimiter('#')) {
                if(stat.keyword("S")) {
                    fromName = "#S";
                } else if(stat.keyword("E")) {
                    fromName = "#E";
                } else {
                    throw ParseError("SequenceParser::parsePosition()",
                                     "Expected 'S' or 'E'.");
                }
            } else {
                fromName = Expressions::parseString
                           (stat, "Expected <name>, \"#S\", or \"#E\".");
            }
            flag |= (flag & 0x0004) ? 0x0008 : 0x0004;
        } else if(defined) {
            // May be element attribute.
            if(Attribute *attr = elem.findAttribute(name)) {
                if(stat.delimiter('=')) {
                    attr->parse(stat, true);
                } else if(stat.delimiter(":=")) {
                    attr->parse(stat, false);
                } else {
                    throw ParseError("SequenceParser::parsePosition()",
                                     "Delimiter \"=\" or \":=\" expected.");
                }
            } else {
                throw ParseError("SequenceParser::parsePosition()",
                                 "Element \"" + elem.getOpalName() +
                                 "\" has no attribute \"" + name + "\".");
            }
        } else {
            throw ParseError("SequenceParser::parsePosition()",
                             "Overriding attributes not permitted here.");
        }
    }

    // Check for valid data and store.
    member.itsPosition = 0.0;

    switch(flag) {

        case 0x0001:
            // "at=...".
            member.itsPosition = position;
            break;

        case 0x0002:
            // "drift=...".  Set distance to drift length.
            member.itsPosition = position;
            member.itsFlag = SequenceMember::DRIFT;
            break;

        case 0x0000:
            // No position data.
            member.itsPosition = 0.0;
            member.itsFlag = SequenceMember::IMMEDIATE;
            break;

        case 0x0005:
            // "at=...,from=...".
            member.itsPosition = position;

            if(fromName == "PREVIOUS") {
                // Position refers to previous element.
                member.itsFlag = SequenceMember::PREVIOUS;
            } else if(fromName == "NEXT") {
                // Position refers to next element.
                member.itsFlag = SequenceMember::NEXT;
            } else if(fromName == "#S") {
                // Position refers to begin of sequence.
                member.itsFlag = SequenceMember::ABSOLUTE;
            } else if(fromName == "#E") {
                // Position refers to end of sequence.
                member.itsFlag = SequenceMember::END;
            } else {
                // Position refers to named element.
                member.itsFlag = SequenceMember::FROM;

                // Find reference to fromName.
                for(RefList::iterator i = references.begin();
                    i != references.end(); ++i) {
                    if(i->fromName == fromName) {
                        // Build reference to a new "FROM" name.
                        i->itsList.push_back(&member);
                        return;
                    }
                }

                // Add reference to an already known "FROM" name.
                Reference ref;
                ref.fromName = fromName;
                ref.itsList.push_back(&member);
                references.push_back(ref);
            }

            break;

        default:
            throw ParseError("SequenceParser:parsePosition()",
                             "Conflicting keywords.");
    }
}
