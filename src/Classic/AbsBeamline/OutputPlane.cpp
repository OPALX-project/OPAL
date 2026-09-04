// Copyright (c) 2023, Chris Rogers
// All rights reserved
//
// This file is part of OPAL.
//
// OPAL is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// You should have received a copy of the GNU General Public License
// along with OPAL. If not, see <https://www.gnu.org/licenses/>.
//

#include "AbsBeamline/OutputPlane.h"

#include "AbsBeamline/BeamlineVisitor.h"
#include "Algorithms/PartBunchBase.h"
#include "Physics/Physics.h"
#include "Structure/LossDataSink.h"
#include "Physics/Units.h"
#include "Utilities/GeneralClassicException.h"

extern Inform *gmsg;

OutputPlane::OutputPlane() : OutputPlane("")
{}

OutputPlane::OutputPlane(const std::string& name):
    PluginElement(name),
    maxIterations_m(10),
    algorithm_m(algorithm::INTERPOLATION),
    verbose_m(0) {
}

OutputPlane::OutputPlane(const OutputPlane& right):
    PluginElement(right),
    maxIterations_m(10) {

    field_m = right.field_m;
    normal_m = right.normal_m;
    centre_m = right.centre_m;
    radialExtent_m = right.radialExtent_m;
    verticalExtent_m = right.verticalExtent_m;
    horizontalExtent_m = right.horizontalExtent_m;
    maxIterations_m = right.maxIterations_m;
    tolerance_m = right.tolerance_m;
    nullfield_m = right.nullfield_m;
    geom_m = right.geom_m;
    recentre_m = right.recentre_m;
    algorithm_m = right.algorithm_m;
    verbose_m = right.verbose_m;
}

OutputPlane::~OutputPlane()
{}

ElementBase* OutputPlane::clone() const {
    ElementBase* element = dynamic_cast<ElementBase*>(new OutputPlane(*this));
    return element;
}

void OutputPlane::accept(BeamlineVisitor& visitor) const {
    visitor.visitOutputPlane(*this);
}

void OutputPlane::doInitialise(PartBunchBase<double, 3>* /*bunch*/) {
    //void OutputPlane::doInitialise() {
    *gmsg << "* Initialize OutputPlane at " << centre_m << " with normal " << normal_m << "\n*";
    if (radialExtent_m >= 0) {
        *gmsg << " Radial extent " << radialExtent_m;
    }
    if (horizontalExtent_m >= 0) {
        *gmsg << " Horizontal extent " << horizontalExtent_m;
    }
    if (verticalExtent_m >= 0) {
        *gmsg << " Vertical extent " << verticalExtent_m;
    }
    if (field_m != nullptr) {
        *gmsg << " Using field map " << field_m->getName();
    } else {
        *gmsg << " Using empty field map";
    }
    if (recentre_m >= 0) {
        *gmsg << " Recentre event " << recentre_m;
    }
    *gmsg << endl;
}

void OutputPlane::doGoOffline() {
    *gmsg << "* OutputPlane goes offline " << getName() << endl;
}

void OutputPlane::RK4Step(
                     const double& tStep,
                     const double& chargeToMass,
                     const double& t,
                     Vector_t& R1,
                     Vector_t& P1) const {
    double thalfStep = tStep / 2.;
    double tPlusHalf = t + thalfStep;
    double tPlusStep = t + tStep;
    // f = dy/dt
    // f1 = f(x,t)
    double  deriv1[6];

    getDerivatives(R1, P1, t, chargeToMass, deriv1);
    // f2 = f(x+dt*f1/2, t+dt/2 )
    double  deriv2[6];
    Vector_t R2({R1[0] + thalfStep * deriv1[0], R1[1] + thalfStep * deriv1[1], R1[2] + thalfStep * deriv1[2]});
    Vector_t P2({P1[0] + thalfStep * deriv1[3], P1[1] + thalfStep * deriv1[4], P1[2] + thalfStep * deriv1[5]});
    getDerivatives(R2, P2, tPlusHalf, chargeToMass, deriv2);

    // f3 = f( x+dt*f2/2, t+dt/2 )
    double  deriv3[6];
    Vector_t R3({R1[0] + thalfStep * deriv2[0], R1[1] + thalfStep * deriv2[1], R1[2] + thalfStep * deriv2[2]});
    Vector_t P3({P1[0] + thalfStep * deriv2[3], P1[1] + thalfStep * deriv2[4], P1[2] + thalfStep * deriv2[5]});
    getDerivatives(R3, P3, tPlusHalf, chargeToMass, deriv3);

    // f4 = f(x+dt*f3, t+dt ).
    double  deriv4[6];
    Vector_t R4({R1[0] + tStep * deriv3[0], R1[1] + tStep * deriv3[1], R1[2] + tStep * deriv3[2]});
    Vector_t P4({P1[0] + tStep * deriv3[3], P1[1] + tStep * deriv3[4], P1[2] + tStep * deriv3[5]});
    getDerivatives(R4, P4, tPlusStep, chargeToMass, deriv4);

    // Return x(t+dt) computed from fourth-order R-K.
    R1[0] += (deriv1[0] + deriv4[0] + 2. * (deriv2[0] + deriv3[0])) * tStep / 6.;
    R1[1] += (deriv1[1] + deriv4[1] + 2. * (deriv2[1] + deriv3[1])) * tStep / 6.;
    R1[2] += (deriv1[2] + deriv4[2] + 2. * (deriv2[2] + deriv3[2])) * tStep / 6.;
    P1[0] += (deriv1[3] + deriv4[3] + 2. * (deriv2[3] + deriv3[3])) * tStep / 6.;
    P1[1] += (deriv1[4] + deriv4[4] + 2. * (deriv2[4] + deriv3[4])) * tStep / 6.;
    P1[2] += (deriv1[5] + deriv4[5] + 2. * (deriv2[5] + deriv3[5])) * tStep / 6.;
}

void OutputPlane::getDerivatives(const Vector_t& R,
                            const Vector_t& P,
                            const double& t,
                            const double& chargeToMass,
                            double* yp) const {

    double gamma = std::sqrt(1 + (P[0] * P[0] + P[1] * P[1] + P[2] * P[2]));

    Vector_t beta = P / gamma;
    double betax = beta[0];
    double betay = beta[1];
    double betaz = beta[2];

    yp[0] = Physics::c * betax;
    yp[1] = Physics::c * betay;
    yp[2] = Physics::c * betaz;

    if (field_m == nullptr) {
        throw GeneralClassicException("OutputPlane::getDerivatives",
                                      "Field was null");
    }
    Vector_t externalB, externalE;
    field_m->apply(R, P, t, externalE, externalB);
    //double kiloGaussToTesla = 0.1;
    externalB *= Units::kG2T;
    externalE *= Units::kV2V / Units::mm2m / Physics::c;

    yp[3] = chargeToMass * (externalB(2) * betay - externalB(1) * betaz+externalE(0));
    yp[4] = chargeToMass * (externalB(0) * betaz - externalB(2) * betax + externalE(1));
    yp[5] = chargeToMass * (externalB(1) * betax - externalB(0) * betay + externalE(2));
}

bool OutputPlane::checkOne(const int index, const double tstep, double chargeToMass, 
                           double& t, Vector_t& R, Vector_t& P) {

    // distance from particle to the output plane
    // time units are ns
    Vector_t delta = R-centre_m;
    double distance = dot(normal_m, delta);

    // maximum step, assuming no curvature, rough guess for relativistic beta
    double betaEstimate = euclidean_norm(P);
    if (betaEstimate > 1) {
        betaEstimate = 1.0;
    }
    double sStep = tstep*betaEstimate*Physics::c;
    if (verbose_m > 3) {
        *gmsg << "* First check crossing of plane " << getName() << " at " << centre_m << " with normal " << normal_m << endl;
        *gmsg << "    Particle " << index << " with R " << R << " P " << P << " t0 " << t << " tstep " << tstep << endl;
        *gmsg << "    Distance prestep " << distance << " compared to s step length " << sStep << endl;
    }
    if (std::abs(distance) > sStep) {
        // we can't cross the plane
        return false;
    }

    Vector_t rTest(R);
    Vector_t pTest(P);
    RK4Step(tstep, chargeToMass, t, rTest, pTest);
    double distanceTest = dot(normal_m, (rTest-centre_m));
    if (verbose_m > 2) {
        *gmsg << "* Second check crossing of plane " << getName() << " at " << centre_m << " with normal " << normal_m << endl;
        *gmsg << "    Particle " << index << " with R " << R << " P " << P << " tstep " << tstep << endl;
        *gmsg << "    After RK4 " << rTest << " " << pTest << endl;
        *gmsg << "    Step distance " << distanceTest << " compared to " << distance << endl;
    }

    if (distance != 0 && distanceTest/distance > 0) {
        // step does not cross the plane - give up
        // note that particle could cross and cross back the plane in a single
        // time-step; in this case the particle is not registered
        return false;
    }

    if (algorithm_m == algorithm::RK4STEP) {
        rk4Test(tstep, chargeToMass, t, R, P);
    } else if (algorithm_m == algorithm::INTERPOLATION) {
        interpolation(t, R, P);
    }

    delta = R-centre_m;
    if (verbose_m > 1) {
        *gmsg << "* Track estimate RK4? " << bool(algorithm_m == algorithm::RK4STEP)
              << " INTERPOLATION? " << bool(algorithm_m == algorithm::INTERPOLATION)
              << " R " << R << " P " << P << " t " << t << endl
              << " delta " << delta << endl;
    }
    if (horizontalExtent_m > 0 && delta[0] * delta[0] + delta[1] * delta[1] > horizontalExtent_m*horizontalExtent_m) {
        // out of horizontal extent. Again, defined in global coordinates
        return false;
    }
    if (verticalExtent_m > 0 && abs(delta[2]) > verticalExtent_m) {
        // out of vertical extent; note that this is defined in global coordinates
        // not in the coordinates of the plane
        return false;
    }
    if (radialExtent_m > 0 &&
        delta[0] * delta[0] + delta[1] * delta[1] + delta[2] * delta[2] > radialExtent_m * radialExtent_m) {
        // out of radial extent
        return false;
    }
    if (recentre_m == index) {
        recentre(R, P);
        *gmsg << "* Recentred output plane to " << centre_m
            << " with normal " << normal_m << " by event " << index << endl;
    }
    if (verbose_m > 0) {
        *gmsg << "* Found track" << endl;
    }
    return true;
}

void OutputPlane::rk4Test(double tstep, double chargeToMass,
                          double& t, Vector_t& R, Vector_t& P) {
    double preStepDistance = 0.0;
    double postStepDistance = dot(normal_m, (R-centre_m));
    size_t iteration = 0;
    while (std::abs(postStepDistance) > tolerance_m && iteration < maxIterations_m) {
        preStepDistance = postStepDistance;
        RK4Step(tstep, chargeToMass, t, R, P); // this updates R and P
        if (verbose_m > 2) {
            *gmsg << "    RK4 iteration " << iteration << " step distance " << preStepDistance << " R " << R << " P " << P
                  << " centre " << centre_m << " d: " << R-centre_m << " t: " << t << " dt: " << tstep << endl;
            Vector_t externalB, externalE;
            field_m->apply(R, P, t, externalE, externalB);
            *gmsg << "        B " << externalB << " [kG] E " << externalE << " [MV/m] " << endl;
        }
        t += tstep;
        postStepDistance = dot(normal_m, (R - centre_m));
        if (postStepDistance/preStepDistance < 0) { // straddling the plane
            // we stepped too far; step in opposite direction
            // step length in ratio of distance to plane
            tstep *= -abs(postStepDistance)/abs(postStepDistance - preStepDistance);
        } else {
            // we didn't step far enough; step in same direction
            tstep *= abs(postStepDistance)/abs(postStepDistance - preStepDistance);
        }
        iteration++;
    }
}

void OutputPlane::interpolation(double& t, Vector_t& R, Vector_t& P) {
    // trajectory R = R0 + V dt 
    // plane = (X-X0).n
    // intersection time t0 = (X0-R0).N / (V.N)
    // intersection position = R0 + V dt
    // relativistic gamma:
    double gamma = std::sqrt(1 + P[0] * P[0] + P[1] * P[1] + P[2] * P[2]);
    Vector_t velocity = P/gamma*Physics::c; // m/ns
    double dt = dot((centre_m-R), normal_m) / dot(velocity, normal_m);
    R += velocity*dt;
    t += dt;
}


bool OutputPlane::doCheck(PartBunchBase<double, 3> *bunch, const int turnnumber,
                          const double t, const double tstep) {
    size_t tempnum = bunch->getLocalNum();
    for(unsigned int i = 0; i < tempnum; ++i) {
        if (verbose_m > 2) {
            *gmsg << "OutputPlane checking at time " << t
                  << " turn number " << turnnumber << " track id " << i << endl;
        }
        Vector_t R(bunch->R[i]);
        Vector_t P(bunch->P[i]);
        double t0(t);
        double chargeToMass = bunch->Q[i] / Physics::q_e * Physics::c * Physics::c
                             / bunch->M[i]/Units::GeV2eV; // electron charge cancels
        bool crossing = checkOne(i, tstep, chargeToMass, t0, R, P);
        if (crossing && lossDs_m) {
            nHits_m += 1;
            lossDs_m->addParticle(OpalParticle(bunch->ID[i], R, P,
                                               t0, bunch->Q[i], bunch->M[i]),
                                  std::make_pair(turnnumber, bunch->bunchNum[i]));
        }
    }
    return false;
}

void OutputPlane::recentre(Vector_t R, Vector_t P) {
    setCentre(R);
    setNormal(P);
    recentre_m = -1; // don't recentre again
}

ElementType OutputPlane::getType() const {
    return ElementType::OUTPUTPLANE;
}
