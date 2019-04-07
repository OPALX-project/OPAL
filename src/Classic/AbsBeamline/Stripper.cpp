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

//change the stripped particles to outcome particles
bool Stripper::doCheck(PartBunchBase<double, 3> *bunch, const int turnnumber, const double t, const double tstep) {

    bool flagNeedUpdate = false;
    bool flagresetMQ = false;
    Vector_t rmin, rmax, strippoint;
    bunch->get_bounds(rmin, rmax);
    // interested in absolute maximum
    double xmax = std::max(std::abs(rmin(0)), std::abs(rmax(0)));
    double ymax = std::max(std::abs(rmin(1)), std::abs(rmax(1)));
    double rbunch_max = std::hypot(xmax, ymax);

    if(rbunch_max > rstart_m - 10.0 ){

        size_t count = 0;
        size_t tempnum = bunch->getLocalNum();
        int pflag = 0;

        Vector_t meanP(0.0, 0.0, 0.0);
        for(unsigned int i = 0; i < bunch->getLocalNum(); ++i) {
          for(int d = 0; d < 3; ++d) {
              meanP(d) += bunch->P[i](d);
          }
        }
        reduce(meanP, meanP, OpAddAssign());
        meanP = meanP / Vector_t(bunch->getTotalNum());

        double sk1, sk2, stangle = 0.0;
        if ( B_m == 0.0 ){
            sk1 = meanP(1)/meanP(0);
            if(sk1 == 0.0)
                stangle =1.0e12;
            else
                stangle = std::abs(1/sk1);
        }else if (meanP(0) == 0.0 ){
            sk2 = - A_m/B_m;
            if(sk2 == 0.0)
              stangle =1.0e12;
            else
              stangle = std::abs(1/sk2);
        }else {
            sk1 = meanP(1)/meanP(0);
            sk2 = - A_m/B_m;
            stangle = std::abs(( sk1-sk2 )/(1 + sk1*sk2));
        }
        double lstep = (sqrt(1.0-1.0/(1.0+dot(meanP, meanP))) * Physics::c) * tstep*1.0e-6; // [mm]
        double Swidth = lstep /  sqrt( 1+1/stangle/stangle ) * 1.2;
        setGeom(Swidth);

        for(unsigned int i = 0; i < tempnum; ++i) {
            if(bunch->PType[i] == ParticleType::REGULAR) {
                pflag = checkPoint(bunch->R[i](0), bunch->R[i](1));
                if(pflag != 0) {
                    // dist1 > 0, right hand, dt > 0; dist1 < 0, left hand, dt < 0
                    double dist1 = (A_m*bunch->R[i](0)+B_m*bunch->R[i](1)+C_m)/R_m/1000.0;
                    double k1, k2, tangle = 0.0;
                    if ( B_m == 0.0 ){
                        k1 = bunch->P[i](1)/bunch->P[i](0);
                        if (k1 == 0.0)
                            tangle = 1.0e12;
                        else
                            tangle = std::abs(1/k1);
                    }else if (bunch->P[i](0) == 0.0 ){
                        k2 = -A_m/B_m;
                        if (k2 == 0.0)
                            tangle = 1.0e12;
                        else
                            tangle = std::abs(1/k2);
                    }else {
                        k1 = bunch->P[i](1)/bunch->P[i](0);
                        k2 = -A_m/B_m;
                        tangle = std::abs(( k1-k2 )/(1 + k1*k2));
                    }
                    double dist2 = dist1 * sqrt( 1+1/tangle/tangle );
                    double dt = dist2/(sqrt(1.0-1.0/(1.0 + dot(bunch->P[i], bunch->P[i]))) * Physics::c)*1.0e9;
                    strippoint(0) = (B_m*B_m*bunch->R[i](0) - A_m*B_m*bunch->R[i](1)-A_m*C_m)/(R_m*R_m);
                    strippoint(1) = (A_m*A_m*bunch->R[i](1) - A_m*B_m*bunch->R[i](0)-B_m*C_m)/(R_m*R_m);
                    strippoint(2) = bunch->R[i](2);
                    lossDs_m->addParticle(strippoint, bunch->P[i], bunch->ID[i], t+dt, turnnumber);

                    if (stop_m) {
                        bunch->Bin[i] = -1;
                        flagNeedUpdate = true;
                    }else{

                        flagNeedUpdate = true;
                        // change charge and mass of PartData when the reference particle hits the stripper.
                        if(bunch->ID[i] == 0)
                          flagresetMQ = true;

                        // change the mass and charge
                        bunch->M[i] = opmass_m;
                        bunch->Q[i] = opcharge_m * Physics::q_e;
                        bunch->PType[i] = ParticleType::STRIPPED;

                        int j = 1;
                        //create new particles
                        while (j < opyield_m){
                          bunch->create(1);
                          bunch->R[tempnum+count] = bunch->R[i];
                          bunch->P[tempnum+count] = bunch->P[i];
                          bunch->Q[tempnum+count] = bunch->Q[i];
                          bunch->M[tempnum+count] = bunch->M[i];
                          // once the particle is stripped, change PType from 0 to 1 as a flag so as to avoid repetitive stripping.
                          bunch->PType[tempnum+count] = ParticleType::STRIPPED;
                          count++;
                          j++;
                        }

                        if(bunch->weHaveBins())
                            bunch->Bin[bunch->getLocalNum()-1] = bunch->Bin[i];

                    }
                }
            }
        }
    }
    reduce(&flagNeedUpdate, &flagNeedUpdate + 1, &flagNeedUpdate, OpBitwiseOrAssign());

    if(!stop_m){
        reduce(&flagresetMQ, &flagresetMQ + 1, &flagresetMQ, OpBitwiseOrAssign());
        if(flagresetMQ){
            bunch->resetM(opmass_m * 1.0e9); // GeV -> eV
            bunch->resetQ(opcharge_m);     // elementary charge
        }
    }

    return flagNeedUpdate;
}

ElementBase::ElementType Stripper::getType() const {
    return STRIPPER;
}
