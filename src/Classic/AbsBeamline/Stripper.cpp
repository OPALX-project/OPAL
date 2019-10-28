#include "AbsBeamline/Stripper.h"

#include "AbsBeamline/BeamlineVisitor.h"
#include "Algorithms/PartBunchBase.h"
#include "Physics/Physics.h"
#include "Structure/LossDataSink.h"

extern Inform *gmsg;

Stripper::Stripper():Stripper("")
{}

Stripper::Stripper(const std::string &name):
    PluginElement(name),
    opcharge_m(0.0),
    opmass_m(0.0),
    opyield_m(1.0),
    stop_m(true)
{}

Stripper::Stripper(const Stripper &right):
    PluginElement(right),
    opcharge_m(right.opcharge_m),
    opmass_m(right.opmass_m),
    opyield_m(right.opyield_m),
    stop_m(right.stop_m)
{}

Stripper::~Stripper() {}

void Stripper::accept(BeamlineVisitor &visitor) const {
    visitor.visitStripper(*this);
}

void Stripper::doFinalise() {
    *gmsg << "* Finalize stripper " << getName() << endl;
}

void Stripper::setOPCharge(double charge) {
    opcharge_m = charge;
}

void Stripper::setOPMass(double mass) {
    opmass_m = mass;
}

void Stripper::setOPYield(double yield) {
    opyield_m = yield;
}

void Stripper::setStop(bool stopflag) {
    stop_m = stopflag;
}

double Stripper::getOPCharge() const {
    return opcharge_m;
}

double Stripper::getOPMass() const {
    return opmass_m;
}

double  Stripper::getOPYield() const {
    return opyield_m;
}

bool  Stripper::getStop () const {
    return stop_m;
}

bool Stripper::doPreCheck(PartBunchBase<double, 3> *bunch) {
    Vector_t rmin, rmax, strippoint;
    bunch->get_bounds(rmin, rmax);
    // interested in absolute maximum
    double xmax = std::max(std::abs(rmin(0)), std::abs(rmax(0)));
    double ymax = std::max(std::abs(rmin(1)), std::abs(rmax(1)));
    double rbunch_max = std::hypot(xmax, ymax);

    if(rbunch_max > rmin_m - 10.0) {
        return true;
    }
    return false;
}

//change the stripped particles to outcome particles
bool Stripper::doCheck(PartBunchBase<double, 3> *bunch, const int turnnumber, const double t, const double tstep) {
    bool flagNeedUpdate = false;
    Vector_t strippoint;

    size_t count = 0;
    size_t tempnum = bunch->getLocalNum();

    Inform gmsgALL("OPAL", INFORM_ALL_NODES);
    for(unsigned int i = 0; i < tempnum; ++i) {
        if(bunch->PType[i] != ParticleType::REGULAR) continue;

        double tangle = calculateIncidentAngle(bunch->P[i](0), bunch->P[i](1));
        changeWidth(bunch, i, tstep, tangle);
        int pflag = checkPoint(bunch->R[i](0), bunch->R[i](1));
        if(pflag == 0) continue;

        // dist1 > 0, right hand, dt > 0; dist1 < 0, left hand, dt < 0
        double dist1 = (A_m * bunch->R[i](0) + B_m * bunch->R[i](1) + C_m) / R_m * 1.0e-3; // [m]
        double dist2 = dist1 * sqrt(1.0 + 1.0 / tangle / tangle);
        double dt = dist2 / (sqrt(1.0 - 1.0 / (1.0 + dot(bunch->P[i], bunch->P[i]))) * Physics::c) * 1.0e9; // [ns]
        strippoint(0) = (B_m * B_m * bunch->R[i](0) - A_m * B_m* bunch->R[i](1) - A_m * C_m) / (R_m * R_m);
        strippoint(1) = (A_m * A_m * bunch->R[i](1) - A_m * B_m* bunch->R[i](0) - B_m * C_m) / (R_m * R_m);
        strippoint(2) = bunch->R[i](2);
        lossDs_m->addParticle(strippoint, bunch->P[i], bunch->ID[i], t+dt,
                              turnnumber, bunch->bunchNum[i]);

        flagNeedUpdate = true;
        if (stop_m) {
            bunch->Bin[i] = -1;
            gmsgALL << level4 << getName() << ": Particle " << bunch->ID[i] << " is deleted by stripper " << getName() << endl;
        } else {
            gmsgALL << level4 << getName() << ": Particle " << bunch->ID[i] << " collide in stripper " << getName() << endl;
            // change charge and mass of PartData when the reference particle hits the stripper.
            if(bunch->ID[i] == 0)
                bunch->setPType(ParticleType::STRIPPED);

            // change the mass and charge
            bunch->M[i] = opmass_m;
            bunch->Q[i] = opcharge_m * Physics::q_e;
            bunch->PType[i] = ParticleType::STRIPPED;

            int j = 1;
            //create new particles
            while (j < opyield_m){
                bunch->create(1);
                size_t index = tempnum + count;
                bunch->R[index] = bunch->R[i];
                bunch->P[index] = bunch->P[i];
                bunch->Q[index] = bunch->Q[i];
                bunch->M[index] = bunch->M[i];
                // once the particle is stripped, change PType from 0 to 1 as a flag so as to avoid repetitive stripping.
                bunch->PType[index] = ParticleType::STRIPPED;
                if(bunch->weHaveBins())
                    bunch->Bin[index] = bunch->Bin[i];

                count++;
                j++;
            }
        }
    }
    return flagNeedUpdate;
}

bool Stripper::doFinaliseCheck(PartBunchBase<double, 3> *bunch, bool flagNeedUpdate) {
    reduce(&flagNeedUpdate, &flagNeedUpdate + 1, &flagNeedUpdate, OpBitwiseOrAssign());

    if(!stop_m){
        // change charge and mass of PartData when the reference particle hits the stripper.
        if (bunch->getPType() == ParticleType::STRIPPED) {
            bunch->resetM(opmass_m * 1.0e9); // GeV -> eV
            bunch->resetQ(opcharge_m);       // elementary charge
        }
    }

    return flagNeedUpdate;
}

ElementBase::ElementType Stripper::getType() const {
    return STRIPPER;
}
