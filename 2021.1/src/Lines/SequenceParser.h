#ifndef OPAL_SequenceParser_HH
#define OPAL_SequenceParser_HH

// ------------------------------------------------------------------------
// $RCSfile: SequenceParser.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: SequenceParser
//
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:33:43 $
// $Author: Andreas Adelmann $
//
// ------------------------------------------------------------------------

#include "OpalParser/OpalParser.h"
#include "Lines/Sequence.h"
#include <list>

class Element;
class Sequence;
class SequenceMember;
class Statement;


// Class SequenceParser
// ------------------------------------------------------------------------
/// The parser for SEQUENCE members.
//  Recognises the format for all possible SEQUENCE members.
//  The SequenceParser is instantiated when a SEQUENCE header is seen.
//  It takes over parsing for the members of the sequence and returns
//  when it sees ENDSEQUENCE.  Before returning it makes sure that all
//  sequence members are placed correctly.

class SequenceParser: public OpalParser {

public:

    /// Constructor.
    //  Assign the current sequence being parsed.
    SequenceParser(Sequence *);

    virtual ~SequenceParser();

    /// Parse sequence member.
    virtual void parse(Statement &) const;

private:

    // The data structure for relative position definitions.
    struct Reference {
        Reference(): fromName()
        { fromPosition = 0; }

        Reference(std::string name, SequenceMember *member): fromName(name)
        { fromPosition = 0; itsList.push_back(member); }

        ~Reference()
        { }

        std::string fromName;
        SequenceMember *fromPosition;
        std::list<SequenceMember *> itsList;
    };

    typedef std::list<Reference> RefList;

    // The type of CLASSIC beamline embedded in a Sequence.
    typedef Sequence::TLine TLine;

    // Not implemented.
    SequenceParser();
    SequenceParser(const SequenceParser &);
    void operator=(const SequenceParser &);

    // Fill in relative positions and convert to absolute.
    void fillPositions() const;

    // Find positions of elements named in "FROM" clauses.
    void findFromPositions() const;

    // Find positions of neighbour elements in sequence line.
    void findNeighbourPositions() const;

    // Parse single member.
    void parseMember(Statement &) const;

    // Parse element position.
    void parsePosition(Statement &, Object &, bool defined) const;

    // The sequence being parsed.
    Sequence* itsSequence;

    // Flag for success.
    mutable bool okFlag;

    // The relative positions to be filled in.
    mutable RefList references;

    // The sequence being parsed.
    mutable TLine itsLine;
};

#endif // OPAL_SequenceParser_HH
