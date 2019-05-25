// ------------------------------------------------------------------------
// $RCSfile: Degrader.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: Degrader
//   Defines the abstract interface for a beam Degrader.
//
// ------------------------------------------------------------------------
// Class category: AbsBeamline
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:32:31 $
// $Author: fci $
//
// ------------------------------------------------------------------------

#include "AbsBeamline/Degrader.h"
#include "Algorithms/PartBunchBase.h"
#include "AbsBeamline/BeamlineVisitor.h"
#include "Fields/Fieldmap.h"
#include "Structure/LossDataSink.h"
#include "Utilities/Options.h"
#include "Solvers/ParticleMatterInteractionHandler.hh"
#include "Physics/Physics.h"
#include "Utilities/Util.h"
#include <memory>

extern Inform *gmsg;

using namespace std;

// Class Degrader
// ------------------------------------------------------------------------

Degrader::Degrader():
    Component(),
    filename_m(""),
    PosX_m(0),
    PosY_m(0),
    PosZ_m(0),
    MomentumX_m(0),
    MomentumY_m(0),
    MomentumZ_m(0),
    time_m(0),
    id_m(0)
{}

Degrader::Degrader(const Degrader &right):
    Component(right),
    filename_m(right.filename_m),
    PosX_m(right.PosX_m),
    PosY_m(right.PosY_m),
    PosZ_m(right.PosZ_m),
    MomentumX_m(right.MomentumX_m),
    MomentumY_m(right.MomentumY_m),
    MomentumZ_m(right.MomentumZ_m),
    time_m(right.time_m),
    id_m(right.id_m)
{}

Degrader::Degrader(const std::string &name):
    Component(name),
    filename_m(""),
    PosX_m(0),
    PosY_m(0),
    PosZ_m(0),
    MomentumX_m(0),
    MomentumY_m(0),
    MomentumZ_m(0),
    time_m(0),
    id_m(0)
{}


Degrader::~Degrader() {

  if(online_m)
    goOffline();
}


void Degrader::accept(BeamlineVisitor &visitor) const {
    visitor.visitDegrader(*this);
}


inline bool Degrader::isInMaterial(double z ) {
 /**
     check if the particle is in the degarder material

  */
    return ((z > 0.0) && (z <= getElementLength()));
}

bool Degrader::apply(const size_t &i, const double &t, Vector_t &E, Vector_t &B) {

    const Vector_t &R = RefPartBunch_m->R[i];
    const Vector_t &P = RefPartBunch_m->P[i];
    const double recpgamma = Physics::c * RefPartBunch_m->dt[i] / sqrt(1.0  + dot(P, P));

    if(isInMaterial(R(2))) {
        //if particle was allready marked as -1 (it means it should have gone into degrader but didn't)
        //set the label to -2 (will not go into degrader and will be deleted when particles per core > 2)
        if (RefPartBunch_m->Bin[i] < 0)
            RefPartBunch_m->Bin[i] = -2;
        else
            RefPartBunch_m->Bin[i] = -1;

        double frac = -R(2) / P(2) * recpgamma;
        PosX_m.push_back(R(0));
        PosY_m.push_back(R(1));
        PosZ_m.push_back(R(2));
        MomentumX_m.push_back(P(0));
        MomentumY_m.push_back(P(1));
        MomentumZ_m.push_back(P(2));
        time_m.push_back(t + frac * RefPartBunch_m->dt[i]);
        id_m.push_back(RefPartBunch_m->ID[i]);
    }

    return false;
}

bool Degrader::applyToReferenceParticle(const Vector_t &R,
                                        const Vector_t &P,
                                        const double &t,
                                        Vector_t &E,
                                        Vector_t &B) {
    if (!isInMaterial(R(2))) return false;

    Vector_t updatedP = P;
    bool isDead = getParticleMatterInteraction()->computeEnergyLoss(updatedP, RefPartBunch_m->getdT(), false);
    double deltaP = euclidean_norm(updatedP) - euclidean_norm(P);
    E(2) += deltaP * RefPartBunch_m->getM() / (RefPartBunch_m->getdT() * RefPartBunch_m->getQ() * Physics::c);

    return isDead;
}

void Degrader::initialise(PartBunchBase<double, 3> *bunch, double &startField, double &endField) {
    RefPartBunch_m = bunch;
    endField = startField + getElementLength();

    if (filename_m == std::string(""))
        lossDs_m = std::unique_ptr<LossDataSink>(new LossDataSink(getName(), !Options::asciidump));
    else
        lossDs_m = std::unique_ptr<LossDataSink>(new LossDataSink(filename_m.substr(0, filename_m.rfind(".")), !Options::asciidump));
}

void Degrader::initialise(PartBunchBase<double, 3> *bunch) {
    RefPartBunch_m = bunch;
    if (filename_m == std::string(""))
        lossDs_m = std::unique_ptr<LossDataSink>(new LossDataSink(getName(), !Options::asciidump));
    else
        lossDs_m = std::unique_ptr<LossDataSink>(new LossDataSink(filename_m.substr(0, filename_m.rfind(".")), !Options::asciidump));
}


void Degrader::finalise()
{
  *gmsg << "* Finalize Degrader" << endl;
}

void Degrader::goOnline(const double &) {
    Inform msg("Degrader::goOnline ");

    PosX_m.reserve((int)(1.1 * RefPartBunch_m->getLocalNum()));
    PosY_m.reserve((int)(1.1 * RefPartBunch_m->getLocalNum()));
    PosZ_m.reserve((int)(1.1 * RefPartBunch_m->getLocalNum()));
    MomentumX_m.reserve((int)(1.1 * RefPartBunch_m->getLocalNum()));
    MomentumY_m.reserve((int)(1.1 * RefPartBunch_m->getLocalNum()));
    MomentumZ_m.reserve((int)(1.1 * RefPartBunch_m->getLocalNum()));
    time_m.reserve((int)(1.1 * RefPartBunch_m->getLocalNum()));
    id_m.reserve((int)(1.1 * RefPartBunch_m->getLocalNum()));
    online_m = true;
}

void Degrader::goOffline() {
    Inform msg("Degrader::goOffline ");
    online_m = false;
    lossDs_m->save();
    msg << " done..." << endl;
}

bool Degrader::bends() const {
    return false;
}

void Degrader::setOutputFN(std::string fn) {
    filename_m = fn;
}

string Degrader::getOutputFN() {
    if (filename_m == std::string(""))
        return getName();
    else
        return filename_m.substr(0, filename_m.rfind("."));
}

void Degrader::getDimensions(double &zBegin, double &zEnd) const {
    zBegin = 0.0;
    zEnd = getElementLength();

}

ElementBase::ElementType Degrader::getType() const {
    return DEGRADER;
}

string Degrader::getDegraderShape() {
    return "DEGRADER";

}