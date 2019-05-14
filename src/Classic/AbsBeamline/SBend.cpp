// ------------------------------------------------------------------------
// $RCSfile: SBend.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Definitions for class: SBend
//   Defines the abstract interface for a sector bend magnet.
//
// ------------------------------------------------------------------------
// Class category: AbsBeamline
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:32:31 $
// $Author: fci $
//
// ------------------------------------------------------------------------

#include "AbsBeamline/SBend.h"
#include "Algorithms/PartBunchBase.h"
#include "AbsBeamline/BeamlineVisitor.h"
#include "Utilities/Options.h"
#include "Fields/Fieldmap.h"
#include "AbstractObjects/OpalData.h"
#include <iostream>
#include <fstream>

extern Inform *gmsg;

// Class SBend
// ------------------------------------------------------------------------

SBend::SBend():
    Bend2D()
{
    setMessageHeader("SBend ");
}

SBend::SBend(const SBend &right):
    Bend2D(right)
{
    setMessageHeader("SBend ");
}

SBend::SBend(const std::string &name):
    Bend2D(name)
{
    setMessageHeader("SBend ");
}

SBend::~SBend() {
}

void SBend::accept(BeamlineVisitor &visitor) const {
    visitor.visitSBend(*this);
}

/*
 * OPAL-MAP methods
 * ================
 */
double SBend::getNormalComponent(int n) const {
    return getField().getNormalComponent(n);
}

double SBend::getSkewComponent(int n) const {
    return getField().getSkewComponent(n);
}

void SBend::setNormalComponent(int n, double v) {
    getField().setNormalComponent(n, v);
}

void SBend::setSkewComponent(int n, double v) {
    getField().setSkewComponent(n, v);
}


/*
 * OPAL-T Methods.
 * ===============
 */

/*
 *  This function merely repackages the field arrays as type Vector_t and calls
 *  the equivalent method but with the Vector_t data types.
 */

ElementBase::ElementType SBend::getType() const {
    return SBEND;
}


bool SBend::findChordLength(Inform &msg,
                            double &chordLength) {

    /*
     * Find bend chord length. If this was not set by the user using the
     * L (length) attribute, infer it from the field map.
     */
    chordLength = getLength();
    return true;
}