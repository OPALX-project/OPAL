#ifndef OPAL_AttWriter_HH
#define OPAL_AttWriter_HH 1

// ------------------------------------------------------------------------
// $RCSfile: AttWriter.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.2 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: AttWriter
//
// ------------------------------------------------------------------------
//
// $Date: 2001/08/13 15:25:21 $
// $Author: jowett $
//
// ------------------------------------------------------------------------

#include "AbstractObjects/Element.h"
#include "Algorithms/DefaultVisitor.h"
#include "Beamlines/FlaggedElmPtr.h"
#include "Elements/OpalElement.h"
#include "Utilities/Timer.h"
#include <iosfwd>
#include <string>
#include <vector>

class AttCell;
class Beamline;
class CCollimator;
class Corrector;
class Drift;
class ElementBase;
class Multipole;
class Patch;
class RBend;
class RFCavity;
class TravelingWave;
class SBend;
class Separator;
class Solenoid;


// Class AttWriter
// ------------------------------------------------------------------------
/// The worker class for ATTLIST commands.
//  A ``Visitor'' which looks at all elements in turn.  For each element
//  it asks to store all defined attributes in a registry, and then
//  requests the desired values from that registry to build a print line.

class AttWriter: public DefaultVisitor {

public:

    /// Constructor.
    AttWriter(const Beamline &,
              std::ostream &,
              const std::vector<AttCell *> &);

    virtual ~AttWriter();


    /// Apply the algorithm to a FlaggedElmPtr.
    virtual void visitFlaggedElmPtr(const FlaggedElmPtr &);

private:

    // Not implemented.
    AttWriter();
    AttWriter(const AttWriter &);
    void operator=(const AttWriter &);

    // The output stream.
    std::ostream &itsStream;

    // The current output line representation.
    const std::vector<AttCell *> &itsBuffer;
};

#endif // OPAL_AttWriter_HH
