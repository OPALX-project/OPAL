//
// Class Undulator
//   Defines all the methods used by the Undulator element.
//   The Undulator element uses a full wave solver from the
//   MITHRA library, see <https://github.com/aryafallahi/mithra/>.
//
// Copyright (c) 2020, Arnau Albà, Paul Scherrer Institut, Villigen PSI, Switzerland
// All rights reserved.
//
// Implemented as part of the MSc thesis
// "Start-to-End Modelling of the AWA Micro-Bunched Electron Cooling POP-Experiment"
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
#include "AbsBeamline/Undulator.h"

#include "AbsBeamline/BeamlineVisitor.h"
#include "Algorithms/PartBunchBase.h"
#include "Physics/Physics.h"
#include "Physics/Units.h"

#include <sys/time.h>
#include <cmath>
#include "mithra/classes.h"
#include "mithra/database.h"
#include "mithra/datainput.h"
#include "mithra/fdtdSC.h"
#include "mithra/fieldvector.h"
#include "mithra/readdata.h"
#include "mithra/stdinclude.h"

extern Inform* gmsg;

Undulator::Undulator() : Undulator("") {
}

Undulator::Undulator(const Undulator& right)
    : Component(right),
      k_m(right.k_m),
      lambda_m(right.lambda_m),
      numPeriods_m(right.numPeriods_m),
      angle_m(right.angle_m),
      fname_m(right.fname_m),
      meshLength_m(right.meshLength_m),
      meshResolution_m(right.meshResolution_m),
      truncationOrder_m(right.truncationOrder_m),
      totalTime_m(right.totalTime_m),
      dtBunch_m(right.dtBunch_m),
      hasBeenSimulated_m(right.hasBeenSimulated_m) {
}

Undulator::Undulator(const std::string& name)
    : Component(name),
      k_m(0.0),
      lambda_m(0.0),
      numPeriods_m(0),
      angle_m(0.0),
      fname_m(""),
      meshLength_m(3, 0.0),
      meshResolution_m(3, 0.0),
      truncationOrder_m(2),
      totalTime_m(0.0),
      dtBunch_m(0.0),
      hasBeenSimulated_m(false) {
}

Undulator::~Undulator() {
}

void Undulator::accept(BeamlineVisitor& visitor) const {
    visitor.visitUndulator(*this);
}

void Undulator::initialise(
    PartBunchBase<double, 3>* bunch, double& /*startField*/, double& /*endField*/) {
    RefPartBunch_m = bunch;
}

void Undulator::apply(
    PartBunchBase<double, 3>* itsBunch, CoordinateSystemTrafo const& refToLocalCSTrafo) {
    Inform msg("MITHRA FW solver ", *gmsg);

    // Get local coordinates w.r.t. undulator.
    const unsigned int localNum = itsBunch->getLocalNum();
    for (unsigned int i = 0; i < localNum; ++i) {
        itsBunch->R[i] = refToLocalCSTrafo.transformTo(itsBunch->R[i]);
        itsBunch->P[i] = refToLocalCSTrafo.rotateTo(itsBunch->P[i]);
    }

    itsBunch->calcBeamParameters();
    msg << "Bunch before undulator in local coordinate system: " << endl;
    itsBunch->print(msg);

    MITHRA::helloMessage();

    // Prepare parameters for full wave solver.
    MITHRA::BunchInitialize bunchInit;
    bunchInit.bunchType_         = "other";
    bunchInit.numberOfParticles_ = itsBunch->getTotalNum();
    bunchInit.cloudCharge_ =
        itsBunch->getTotalNum() * itsBunch->getChargePerParticle() / (-Physics::q_e);
    bunchInit.initialGamma_ = itsBunch->get_gamma();
    for (unsigned int d = 0; d < 3; ++d)
        bunchInit.initialDirection_[d] = itsBunch->get_pmean()[d];
    bunchInit.initialDirection_ /= euclidean_norm(itsBunch->get_pmean());
    MITHRA::Bunch bunch;
    bunch.bunchInit_.push_back(bunchInit);
    bunch.timeStep_ = getDtBunch();
    msg << "Bunch parameters have been transferred to the full-wave solver." << endl;

    MITHRA::Undulator undulator;
    undulator.k_      = getK();
    undulator.lu_     = getLambda();
    undulator.length_ = getNumPeriods();
    undulator.theta_  = getAngle();
    double lFringe    = 2 * undulator.lu_;                     // Fringe field length is 2*lu.
    undulator.dist_ = lFringe - itsBunch->get_maxExtent()[2];  // Bunch-head to undulator distance.
    std::vector<MITHRA::Undulator> undulators;
    undulators.push_back(undulator);
    msg << "Undulator parameters have been transferred to the full-wave solver." << endl;

    MITHRA::Mesh mesh;
    mesh.initialize();
    mesh.lengthScale_    = 1.0;  // OPAL uses metres.
    mesh.timeScale_      = 1.0;  // OPAL uses seconds.
    mesh.meshCenter_     = MITHRA::FieldVector<double>(0.0);
    mesh.meshLength_     = getMeshLength();
    mesh.meshResolution_ = getMeshResolution();
    mesh.totalTime_      = getTotalTime();
    // If simulation time is not given, run the full-wave solver until the end of the undulator.
    if (mesh.totalTime_ == 0)
        mesh.totalDist_ = lFringe + undulator.lu_ * undulator.length_;
    mesh.truncationOrder_  = getTruncationOrder();
    mesh.spaceCharge_      = true;
    mesh.optimizePosition_ = true;
    msg << "Mesh parameters have been transferred to the full-wave solver." << endl;

    MITHRA::Seed seed;
    seed.a0_ = 0.0; // initialise unitialised member (see #658)
    std::vector<MITHRA::ExtField> externalFields;
    std::vector<MITHRA::FreeElectronLaser> FELs;

    // Get filename with desired output data.
    if (!fname_m.empty()) {
        std::list<std::string> jobFile = MITHRA::read_file(fname_m.c_str());
        MITHRA::cleanJobFile(jobFile);
        MITHRA::ParseDarius parser(jobFile, mesh, bunch, seed, undulators, externalFields, FELs);
        parser.setJobParameters();
    }

    MITHRA::FdTdSC solver(mesh, bunch, seed, undulators, externalFields, FELs);

    // Transfer particles to full wave solver and destroy them from itsBunch.
    MITHRA::Charge charge;
    charge.q = itsBunch->getChargePerParticle() / (-Physics::q_e);
    for (unsigned int i = 0; i < localNum; ++i) {
        for (unsigned int d = 0; d < 3; ++d) {
            charge.rnp[d]  = itsBunch->R[i][d];
            charge.gbnp[d] = itsBunch->P[i][d];
        }
        solver.chargeVectorn_.push_back(charge);
    }
    itsBunch->destroy(localNum, 0, true);
    msg << "Particles have been transferred to the full-wave solver." << endl;

    // Print the parameters for the simulation.
    mesh.show();
    bunch.show();
    seed.show();
    for (unsigned int i = 0; i < undulators.size(); i++) {
        undulators[i].show();
    }
    for (unsigned int i = 0; i < externalFields.size(); i++) {
        externalFields[i].show();
    }

    // Run the full-wave solver.
    timeval simulationStart;
    gettimeofday(&simulationStart, nullptr);
    solver.solve();

    // Get total computational time of the full wave simulation.
    timeval simulationEnd;
    gettimeofday(&simulationEnd, nullptr);
    double deltaTime = (simulationEnd.tv_usec - simulationStart.tv_usec) * Units::us2s;
    deltaTime += (simulationEnd.tv_sec - simulationStart.tv_sec);
    msg << "::: Total full wave simulation time [seconds] = " << deltaTime << endl;

    // Lorentz Transformation back to undulator local coordinates.
    // First you need to get the position of the bunch tail.
    double zMin = 1e100;
    for (auto iter = solver.chargeVectorn_.begin(); iter != solver.chargeVectorn_.end(); iter++) {
        zMin = std::min(zMin, iter->rnp[2]);
    }
    allreduce(&zMin, 1, std::less<double>());

    const double gammaBeta = solver.gamma_ * solver.beta_;
    const double factor    = gammaBeta * solver.c0_ * (solver.timeBunch_ + solver.dt_) + lFringe;
    for (auto iter = solver.chargeVectorn_.begin(); iter != solver.chargeVectorn_.end(); iter++) {
        double dist = zMin - iter->rnp[2];
        // Lorentz transform.
        iter->rnp[2] = solver.gamma_ * iter->rnp[2] + factor;
        iter->gbnp[2] =
            solver.gamma_ * iter->gbnp[2] + gammaBeta * std::sqrt(1 + iter->gbnp.norm2());
        // Shift to bring all particles to same time in lab frame.
        double gammaParticle = std::sqrt(1 + iter->gbnp.norm2());
        iter->rnp[0] += iter->gbnp[0] / gammaParticle * dist * gammaBeta;
        iter->rnp[1] += iter->gbnp[1] / gammaParticle * dist * gammaBeta;
        iter->rnp[2] += iter->gbnp[2] / gammaParticle * dist * gammaBeta;
    }

    // Get total time elapsed in laboratory frame.
    mesh.totalTime_ =
        solver.gamma_ * (solver.time_ + solver.beta_ / solver.c0_ * (zMin - bunch.zu_));

    // Return particles to itsBunch in local coordinates.
    msg << "Transferring particles back to OPAL bunch." << endl;
    itsBunch->create(solver.chargeVectorn_.size());
    const double dt                          = itsBunch->getdT();
    const unsigned int newLocalNum           = itsBunch->getLocalNum();
    std::list<MITHRA::Charge>::iterator iter = solver.chargeVectorn_.begin();
    for (unsigned int i = 0; i < newLocalNum; ++i) {
        for (unsigned int d = 0; d < 3; ++d) {
            itsBunch->R[i][d] = iter->rnp[d];
            itsBunch->P[i][d] = iter->gbnp[d];
        }
        itsBunch->Q[i]  = iter->q * (-Physics::q_e);
        itsBunch->dt[i] = dt;
        iter++;
    }
    itsBunch->setT(itsBunch->getT() + mesh.totalTime_);

    // Transform back to reference coordinate system.
    CoordinateSystemTrafo localToRefCSTrafo = refToLocalCSTrafo.inverted();
    for (unsigned int i = 0; i < newLocalNum; ++i) {
        itsBunch->R[i] = localToRefCSTrafo.transformTo(itsBunch->R[i]);
        itsBunch->P[i] = localToRefCSTrafo.rotateTo(itsBunch->P[i]);
    }
    itsBunch->calcBeamParameters();

    // Update reference particle.
    // The reference particle becomes the bunch-centroid after the undulator.
    itsBunch->RefPartR_m = itsBunch->toLabTrafo_m.transformTo(itsBunch->get_centroid());
    itsBunch->RefPartP_m = itsBunch->toLabTrafo_m.rotateTo(itsBunch->get_pmean());

    msg << "Bunch after undulator in reference coordinate system: " << endl;
    itsBunch->print(msg);

    setHasBeenSimulated(true);
}

void Undulator::finalise() {
}

bool Undulator::bends() const {
    return false;
}

void Undulator::getDimensions(double& /*zBegin*/, double& /*zEnd*/) const {
}

ElementType Undulator::getType() const {
    return ElementType::UNDULATOR;
}

void Undulator::setK(double k) {
    k_m = k;
}
double Undulator::getK() const {
    return k_m;
}

void Undulator::setLambda(double lambda) {
    lambda_m = lambda;
}
double Undulator::getLambda() const {
    return lambda_m;
}

void Undulator::setNumPeriods(unsigned int np) {
    numPeriods_m = np;
}
unsigned int Undulator::getNumPeriods() const {
    return numPeriods_m;
}

void Undulator::setAngle(double theta) {
    angle_m = theta;
}
double Undulator::getAngle() const {
    return angle_m;
}

void Undulator::setFilename(const std::string& fname) {
    fname_m = fname;
}
const std::string& Undulator::getFilename() const {
    return fname_m;
}

void Undulator::setMeshLength(const std::vector<double>& ml) {
    meshLength_m = ml;
}
std::vector<double> Undulator::getMeshLength() const {
    return meshLength_m;
}

void Undulator::setMeshResolution(const std::vector<double>& mr) {
    meshResolution_m = mr;
}
std::vector<double> Undulator::getMeshResolution() const {
    return meshResolution_m;
}

void Undulator::setTruncationOrder(unsigned int trunOrder) {
    truncationOrder_m = trunOrder;
}
unsigned int Undulator::getTruncationOrder() const {
    return truncationOrder_m;
}

void Undulator::setTotalTime(double tt) {
    totalTime_m = tt;
}
double Undulator::getTotalTime() const {
    return totalTime_m;
}

void Undulator::setDtBunch(double dtb) {
    dtBunch_m = dtb;
}
double Undulator::getDtBunch() const {
    return dtBunch_m;
}

void Undulator::setHasBeenSimulated(bool hbs) {
    hasBeenSimulated_m = hbs;
}
bool Undulator::getHasBeenSimulated() const {
    return hasBeenSimulated_m;
}
