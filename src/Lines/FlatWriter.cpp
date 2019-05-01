// ------------------------------------------------------------------------
// $RCSfile: FlatWriter.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.2.4.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: FlatWriter
//   The class for the OPAL FlatWriter command.
//
// ------------------------------------------------------------------------
//
// Revision History:
// $Date: 2002/12/09 15:06:08 $
// $Author: jsberg $
//
// JMJ: removing dirty C-style hacking where splitName is constructed.
//      13/10/2000
// ------------------------------------------------------------------------

#include "Lines/FlatWriter.h"
#include "AbsBeamline/Drift.h"
#include "AbsBeamline/ElementBase.h"
#include "AbstractObjects/Attribute.h"
#include "AbstractObjects/Element.h"
#include "AbstractObjects/OpalData.h"
#include "AbstractObjects/Object.h"
#include "Algorithms/MPSplitIntegrator.h"
#include "Attributes/Attributes.h"
#include "Beamlines/Beamline.h"
#include "Beamlines/FlaggedElmPtr.h"
#include <vector>
#include <sstream>

extern Inform *gmsg;

// Class FlatWriter.
// ------------------------------------------------------------------------
// Visitor class for writing a line or sequence as a flat sequence in
// OPAL-8 format.

FlatWriter::FlatWriter(const Beamline &bl, Sequence::TLine &seq,
                       std::ostream &os):
    DefaultVisitor(bl, false, false),
    itsPosition(0.0), itsSequence(seq), itsStream(os)
{}


FlatWriter::~FlatWriter()
{}


const Sequence::TLine &FlatWriter::getSequence() const {
    return itsSequence;
}


void FlatWriter::visitDrift(const Drift &drf) {
    // All drifts are skipped, but their length is accumulated.
    itsPosition += drf.getElementLength();
}


void FlatWriter::visitMapIntegrator(const MapIntegrator &i) {
    ElementBase *base = i.getElement();

    if(const MPSplitIntegrator *mpi =
           dynamic_cast<const MPSplitIntegrator *>(&i)) {
        // For a MPSliceIntegrator store the equivalent thin lenses.
        double length = base->getElementLength();
        std::vector<double> itsSlices;
        mpi->getSlices(itsSlices);
        int slices = itsSlices.size() - 1;

        if(slices > 1  &&  length != 0.0) {
            // Construct internal element name.
            // JMJ modified following to get rid of char stuff 13/10/2000
            // This removes leading zeros in slice number and also changes
            // underscore character to a pair of dots.
            const std::string name = base->getName();

            // Store the slices in the temporary sequence.
            // JMJ 18/12/2000, put name generation inside loop
            // JMJ 19/12/2000, put declaration of splitNameStream inside loop so
            //                 that it would be re-initialised each pass.
            for(int k = 0; k < slices; ++k) {

                std::ostringstream splitNameStream;
                splitNameStream << name << ".." << k + 1 << std::ends ;
                std::string splitName = splitNameStream.str();

                // Find or construct internal element.
                Object *split = OpalData::getInstance()->find(splitName);
                if(split == 0) {
                    split = OpalData::getInstance()->find(name)->clone(splitName);
                    Attribute *attr = split->findAttribute("L");
                    Attributes::setReal(*attr, length / double(slices));
                    OpalData::getInstance()->define(split);
                    itsStream << *split;
                }
                SequenceMember member;
                member.setElement(static_cast<Element *>(split)->getElement());
                member.itsType = SequenceMember::GLOBAL;
                member.itsPosition = itsPosition + itsSlices[k] * length;
                itsSequence.push_back(member);

            }

            itsPosition += length;
        } else {
            FlatWriter::applyDefault(*base);
        }
    } else {
        FlatWriter::applyDefault(*base);
    }
}


void FlatWriter::applyDefault(const ElementBase &elem) {
    const std::string name = elem.getName();

    if(name[0] != '#') {
        SequenceMember member;
        member.setElement(const_cast<ElementBase *>(&elem));
        member.itsType = SequenceMember::GLOBAL;
        double length = elem.getElementLength();
        member.itsPosition = itsPosition + length / 2.0;
        itsSequence.push_back(member);
        itsPosition += length;
    }
}
