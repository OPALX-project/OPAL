// ------------------------------------------------------------------------
// $RCSfile: AttWriter.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.4 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: AttWriter
//   The worker class for OPAL ATTLIST commands.
//
// ------------------------------------------------------------------------
//
// $Date: 2002/03/28 21:27:54 $
// $Author: jsberg $
//
// ------------------------------------------------------------------------

#include "Tables/AttWriter.h"

#include "AbsBeamline/CCollimator.h"
#include "AbsBeamline/FlexibleCollimator.h"
#include "AbsBeamline/Corrector.h"
#include "AbsBeamline/Drift.h"
#include "AbsBeamline/ElementBase.h"
#include "AbsBeamline/Multipole.h"
#include "AbsBeamline/Patch.h"
#include "AbsBeamline/RBend.h"
#include "AbsBeamline/RFCavity.h"
#include "AbsBeamline/TravelingWave.h"
#include "AbsBeamline/SBend.h"
#include "AbsBeamline/Separator.h"
#include "AbsBeamline/Solenoid.h"

#include "AbstractObjects/BeamSequence.h"
#include "AbstractObjects/Element.h"
#include "AbstractObjects/OpalData.h"
#include "Algorithms/DefaultVisitor.h"
#include "Beamlines/Beamline.h"
#include "Beamlines/FlaggedElmPtr.h"
#include "Elements/AttCell.h"
#include "Elements/OpalElement.h"
#include "Physics/Physics.h"
#include "Utilities/Timer.h"
#include <fstream>
#include <vector>


// Class AttWriter
// ------------------------------------------------------------------------

AttWriter::AttWriter(const Beamline &line,
                     std::ostream &os,
                     const std::vector<AttCell *> &buffer):
    DefaultVisitor(line, false, false),
    itsStream(os),
    itsBuffer(buffer)
{}


AttWriter::~AttWriter()
{}


void AttWriter::visitFlaggedElmPtr(const FlaggedElmPtr &fep) {
    ElementBase *base = fep.getElement();
    const std::string &nam = base->getName();
    if(dynamic_cast<Beamline *>(base)) {
        DefaultVisitor::visitFlaggedElmPtr(fep);
    } else if(fep.getSelectionFlag()) {
        // Fill the line buffer.
        if(nam[0] == '[') {
            OpalElement::setRegisteredAttribute("L", base->getElementLength());
            OpalElement::setRegisteredAttribute("KEYWORD", "DRIFT");
        } else {
            OpalElement *elem = dynamic_cast<OpalElement *>(OpalData::getInstance()->find(nam));
            elem->fillRegisteredAttributes(*fep.getElement());
        }

        // Write the current output line and clear it.
        itsStream << ' ';
        std::vector<AttCell *>::size_type n = itsBuffer.size();
        std::vector<AttCell *>::size_type i;
        for(i = 0; i < n; ++i) {
            itsStream << ' ';
            itsBuffer[i]->printValue(itsStream);
            itsBuffer[i]->clearValue();
        }
        itsStream << '\n';
    }
}