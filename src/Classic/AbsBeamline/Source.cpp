
#include "AbsBeamline/Source.h"
#include "Algorithms/PartBunchBase.h"
#include "AbsBeamline/BeamlineVisitor.h"
#include "Fields/Fieldmap.h"
#include "Physics/Physics.h"
#include "Elements/OpalBeamline.h"

#include <iostream>
#include <fstream>

extern Inform *gmsg;

// Class Source
// ------------------------------------------------------------------------

const std::string Source::defaultDistribution("DISTRIBUTION");

Source::Source():
    Source("")
{}


Source::Source(const Source &right):
    Component(right),
    ElementEdge_m(right.ElementEdge_m),
    startField_m(right.startField_m),
    endField_m(right.endField_m)
{
    setElType(isSource);
}


Source::Source(const std::string &name):
    Component(name),
    ElementEdge_m(0.0),
    startField_m(0.0),
    endField_m(0.0)
{
    setElType(isSource);
}

Source::~Source() {

}


void Source::accept(BeamlineVisitor &visitor) const {
    visitor.visitSource(*this);
}

/**
 * \brief ENVELOPE COMPONENT for radial focussing of the beam
 * Calculates the transverse envelope component for the
 * solenoid element and adds it to the K vector
*/
void Source::addKR(int i, double t, Vector_t &K) {
}

/**
 * ENVELOPE COMPONENT for transverse kick (only important for x0, y0)
 * Calculates the transverse kick component for the solenoid element and adds it to
 * the K vector, only important for off track tracking. Otherwise KT = 0.
*/
void Source::addKT(int i, double t, Vector_t &K) {
}

bool Source::apply(const double &t) {
    Vector_t externalE = Vector_t(0.0);
    Vector_t externalB = Vector_t(0.0);
    beamline_m->getFieldAt(csTrafoGlobal2Local_m.getOrigin(),
                           Vector_t(0.0),
                           RefPartBunch_m->getT() + 0.5 * RefPartBunch_m->getdT(),
                           externalE,
                           externalB);
    for (Distribution* dist: distrs_m) {
        dist->emitParticles(RefPartBunch_m, externalE(2));
    }

    return false;
}

void Source::initialise(PartBunchBase<double, 3> *bunch, double &startField, double &endField) {
    // Inform msg("Source ", *gmsg);
    // double zBegin = 0.0, zEnd = 0.0, rBegin = 0.0, rEnd = 0.0;

    RefPartBunch_m = bunch;
}

void Source::finalise()
{}

bool Source::bends() const {
    return false;
}


void Source::goOnline(const double &) {
    online_m = true;
}

void Source::goOffline() {
    online_m = false;
}

void Source::getDimensions(double &zBegin, double &zEnd) const {
    zBegin = startField_m;
    zEnd = endField_m;
}


ElementBase::ElementType Source::getType() const {
    return SOURCE;
}

void Source::setDistribution(std::vector<std::string> distNames) {
    const size_t numberOfDistributions = distNames.size();
    if (numberOfDistributions == 0) {
        distrs_m.push_back(Distribution::find(defaultDistribution));
    } else {
        for (std::string name: distNames) {
            distrs_m.push_back(Distribution::find(name));
            distrs_m.back()->setNumberOfDistributions(numberOfDistributions);
        }
    }
}