#ifndef OPAL_BeamSequence_HH
#define OPAL_BeamSequence_HH

// ------------------------------------------------------------------------
// $RCSfile: BeamSequence.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: BeamSequence
//
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:33:34 $
// $Author: Andreas Adelmann $
//
// ------------------------------------------------------------------------

#include "AbstractObjects/Element.h"

class Beamline;


// Class BeamSequence
// ------------------------------------------------------------------------
/// The base class for all OPAL beam lines and sequences.
//  It implements the common behaviour of sequences, it can also be used via
//  dynamic casting to determine whether an object represents a sequence.

class BeamSequence: public Element {

public:

    virtual ~BeamSequence();

    /// Make complete copy.
    //  Copy also the line list.
    virtual BeamSequence *copy(const std::string &name) = 0;

    /// Find a BeamSequence by name.
    static BeamSequence *find(const std::string &name);

    /// Return the object category as a string.
    //  Return the string "SEQUENCE".
    virtual const std::string getCategory() const;

    /// Return the embedded CLASSIC beam line.
    //  The result it the ideal line.
    virtual Beamline *fetchLine() const = 0;

protected:

    /// Constructor for exemplars.
    BeamSequence(int size, const char *name, const char *help);

    /// Constructor for clones.
    //  The clone will be [b]empty[/b].
    //  It has to be filled in by the corresponding parser.
    BeamSequence(const std::string &name, BeamSequence *parent);

private:

    // Not implemented.
    BeamSequence();
    BeamSequence(const BeamSequence &);
    void operator=(const BeamSequence &);
};

#endif // OPAL_BeamSequence_HH
