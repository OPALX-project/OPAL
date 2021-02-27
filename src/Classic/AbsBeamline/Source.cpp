
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
    startField_m(right.startField_m),
    endField_m(right.endField_m),
    isTransparent_m(right.isTransparent_m)
{}


Source::Source(const std::string &name):
    Component(name),
    startField_m(0.0),
    endField_m(0.0),
    isTransparent_m(false)
{}

Source::~Source() {
}


void Source::accept(BeamlineVisitor &visitor) const {
    visitor.visitSource(*this);
}

bool Source::apply(const size_t &i, const double &t, Vector_t &/*E*/, Vector_t &/*B*/) {
    if (isTransparent_m) {
        return false;
    }

    const Vector_t &R = RefPartBunch_m->R[i];
    const Vector_t &P = RefPartBunch_m->P[i];
    const double &dt = RefPartBunch_m->dt[i];
    const double recpgamma = Physics::c * dt / Util::getGamma(P);
    if (online_m && R(2) <= 0.0 && P(2) < 0.0) {
        double frac = -R(2) / (P(2) * recpgamma);

        lossDs_m->addParticle(OpalParticle(RefPartBunch_m->ID[i],
                                           R + frac * recpgamma * P, P,
                                           t + frac * dt,
                                           RefPartBunch_m->Q[i], RefPartBunch_m->M[i],
                                           OpalParticle::SPATIAL),
                              std::make_pair(0, RefPartBunch_m->bunchNum[i]));

        return true;
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

void Source::setTransparent() {
    isTransparent_m = true;
}