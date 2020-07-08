//
// Class Sequence
//   The SEQUENCE definition.
//   A Sequence contains a CLASSIC TBeamline<SequenceMember> which represents
//   the sequence of elements in the line and their positions.  The sequence
//   definition is parsed in collaboration with a SequenceParser.
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
#ifndef OPAL_Sequence_HH
#define OPAL_Sequence_HH

#include "AbstractObjects/BeamSequence.h"
#include "AbsBeamline/ElementBase.h"
#include "Beamlines/TBeamline.h"
#include "Lines/SequenceMember.h"
#include "MemoryManagement/Pointer.h"
#include <list>

class TokenStream;

class Sequence: public BeamSequence {

    friend class Edit;
    friend class SequenceParser;
    friend class SequenceTemplate;

public:

    /// The type of a sequence line.
    typedef TBeamline<SequenceMember> TLine;

    /// Exemplar constructor.
    Sequence();

    virtual ~Sequence();

    /// Make clone.
    // The new object is an empty sequence, it will be filled in by the parser.
    virtual Sequence *clone(const std::string &name);

    /// Make copy of the sequence line.
    virtual Sequence *copy(const std::string &name);

    /// Return sequence length.
    virtual double getLength() const;

    /// Return the arc length from origin to entrance.
    virtual double getEntrance(ReferenceType) const;

    /// Return the arc length from origin to exit.
    virtual double getExit(ReferenceType) const;

    /// Return the reference type flag.
    //  "ENTRY", "CENTRE", or "EXIT".
    ReferenceType getReference() const;

    /// Make a sequence template.
    //  Arguments:
    //  [ol]
    //  [li] Name to be given to the template.
    //  [li] Token stream to be read for the sequence line.
    //  [li] Statement to be read for the arguments.
    //  [/ol]
    virtual Object *makeTemplate(const std::string &, TokenStream &, Statement &);

    /// Parse sequence.
    virtual void parse(Statement &);

    /// Print sequence.
    virtual void print(std::ostream &) const;

    /// Replace references to elements.
    virtual void replace(Object *oldObject, Object *newObject);

    /// Update the embedded CLASSIC beam line.
    //  Recompute drift lengths.
    virtual void update();

    /// Return the embedded CLASSIC beam line.
    //  The result it the ideal line.
    virtual TLine *fetchLine() const;

    /// Store sequence line.
    //  Assign to the underlying ideal line and re-insert the drifts.
    void storeLine(TLine &line);

private:

    // Not implemented.
    Sequence(const Sequence &);
    void operator=(const Sequence &);

    // Clone constructor.
    Sequence(const std::string &name, Sequence *parent);

    // Add the top-level begin and end markers.
    void addEndMarkers(TLine &line) const;

    // Compute drift length for a generated drift.
    double findDriftLength(TLine::iterator drift) const;

    // Find named position in sequence line.
    TLine::iterator findNamedPosition(TLine &, const std::string &) const;

    // Insert the top-level drift spaces.
    void insertDrifts(TLine &line);

    //  Recompute drift lengths for a given sequence or sub-sequence.
    static void updateList(Sequence *, TLine *);

    // The reference code.
    ReferenceType itsCode;

    // The reference point.
    std::string itsRefPoint;
};

#endif // OPAL_Sequence_HH
