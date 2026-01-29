// ------------------------------------------------------------------------
// $RCSfile: BeamSequence.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: BeamSequence
//   The abstract base class for all OPAL beam sequences (LINE and SEQUENCE).
//
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:33:34 $
// $Author: Andreas Adelmann $
//
// ------------------------------------------------------------------------

#include "AbstractObjects/BeamSequence.h"
#include "AbstractObjects/OpalData.h"
#include "Utilities/OpalException.h"


// Class BeamSequence
// ------------------------------------------------------------------------

BeamSequence::~BeamSequence()
{}


const std::string BeamSequence::getCategory() const {
    return "SEQUENCE";
}


BeamSequence *BeamSequence::find(const std::string &name) {
    OpalData *opal = OpalData::getInstance();
    BeamSequence *bs = dynamic_cast<BeamSequence *>(opal->find(name));
    if(bs == 0) {
        throw OpalException("BeamSequence::find()",
                            "Beam line or sequence \"" + name + "\" not found.");
    }
    return bs;
}


BeamSequence::BeamSequence(int size, const char *name, const char *help):
    Element(size, name, help)
{}


BeamSequence::BeamSequence(const std::string &name, BeamSequence *parent):
    Element(name, parent)
{}
