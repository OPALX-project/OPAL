
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

bool Source::apply(const size_t &i, const double &t, Vector_t &E, Vector_t &B) {
    const Vector_t &R = RefPartBunch_m->R[i];
    const Vector_t &P = RefPartBunch_m->P[i];
    const double &dt = RefPartBunch_m->dt[i];
    const double recpgamma = Physics::c * dt / Util::getGamma(P);
    const double end = getElementLength();
    if (online_m && dt < 0.0) {
        if (R(2) > end &&
            (R(2) + P(2) * recpgamma) < end) {
            double frac = (end - R(2)) / (P(2) * recpgamma);

            lossDs_m->addParticle(R + frac * recpgamma * P,
                                  P, RefPartBunch_m->ID[i], t + frac * dt, 0);

            return true;
        }
    }

    return false;
}

void Source::initialise(PartBunchBase<double, 3> *bunch, double &startField, double &endField) {
    RefPartBunch_m = bunch;
    endField = startField;
    startField -= getElementLength();

    std::string filename = getName();
    lossDs_m = std::unique_ptr<LossDataSink>(new LossDataSink(filename,
                                                              !Options::asciidump,
                                                              getType()));

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
    lossDs_m->save();
}

void Source::getDimensions(double &zBegin, double &zEnd) const {
    zBegin = startField_m;
    zEnd = endField_m;
}


ElementBase::ElementType Source::getType() const {
    return SOURCE;
}