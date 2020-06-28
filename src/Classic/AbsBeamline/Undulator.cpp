// ------------------------------------------------------------------------
// $RCSfile: Undulator.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: Undulator
//   Defines the abstract interface for a drift space.
//
// ------------------------------------------------------------------------
// Class category: AbsBeamline
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:32:31 $
// $Author: fci $
//
// ------------------------------------------------------------------------

#include "AbsBeamline/Undulator.h"
#include "AbsBeamline/BeamlineVisitor.h"
#include "Algorithms/PartBunchBase.h"

#ifdef OPAL_FEL
#include "mithra/fieldvector.h"
#include "mithra/stdinclude.h"
#include "mithra/readdata.h"
#include "mithra/database.h"
#include "mithra/classes.h"
#include "mithra/datainput.h"
#include "mithra/readdata.h"
#include "mithra/solver.h"
#include "mithra/fdtdSC.h"
#include <ctime>
#include <cmath>
#endif

extern Inform *gmsg;

// Class Undulator
// ------------------------------------------------------------------------

Undulator::Undulator():
    Component(),
    nSlices_m(1),
    k_m(0.0),
    lambda_m(0.0),
    numPeriods_m(0.0),
    angle_m(0.0),
    fname_m(""),
    meshLength_m(3, 0.0),
    meshResolution_m(3, 0.0),
    truncationOrder_m(2),
    totalTime_m(0.0),
    dtBunch_m(0.0),
    hasBeenSimulated_m(0)
{ }


Undulator::Undulator(const Undulator &right):
    Component(right),
    nSlices_m(right.nSlices_m),
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
    hasBeenSimulated_m(right.hasBeenSimulated_m)
{ }


Undulator::Undulator(const std::string &name):
    Component(name),
    nSlices_m(1),
    k_m(0.0),
    lambda_m(0.0),
    numPeriods_m(0.0),
    angle_m(0.0),
    fname_m(""),
    meshLength_m(3, 0.0),
    meshResolution_m(3, 0.0),
    truncationOrder_m(2),
    totalTime_m(0.0),
    dtBunch_m(0.0),
    hasBeenSimulated_m(0)
{ }


Undulator::~Undulator()
{ }


void Undulator::accept(BeamlineVisitor &visitor) const {
    visitor.visitUndulator(*this);
}

void Undulator::initialise(PartBunchBase<double, 3> *bunch, double &startField, double &endField) {
    endField = startField + getElementLength();
    RefPartBunch_m = bunch;
    startField_m = startField;
}

#ifdef OPAL_FEL
void Undulator::apply(PartBunchBase<double, 3> *itsBunch, CoordinateSystemTrafo const& refToLocalCSTrafo) {
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

    // Mithra full wave solver initial message.
    msg << " :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::" << endl;
    msg << " MITHRA-2.0: Completely Numerical Calculation of Free Electron Laser Radiation" << endl;
    msg << " Version 2.0, Copyright 2019, Arya Fallahi" << endl;
    msg << " Written by Arya Fallahi, IT'IS Foundation, Zurich, Switzerland" << endl;
    msg << " :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::" << endl;

    // Prepare parameters for full wave solver.
    MITHRA::BunchInitialize bunchInit;
    bunchInit.bunchType_ = "other";
    bunchInit.numberOfParticles_ = itsBunch->getTotalNum();
    bunchInit.cloudCharge_ = itsBunch->getTotalNum() * itsBunch->getChargePerParticle() / (-Physics::q_e);
    bunchInit.initialGamma_ = itsBunch->get_gamma();
    for (unsigned int d = 0; d < 3; ++d) 
        bunchInit.initialDirection_[d] = itsBunch->get_pmean()[d];
    bunchInit.initialDirection_ /= euclidean_norm(itsBunch->get_pmean());
    MITHRA::Bunch bunch;
    bunch.bunchInit_.push_back(bunchInit);
    bunch.timeStep_ = getDtBunch();
    msg << "Bunch parameters have been transferred to the full-wave solver." << endl;

    MITHRA::Undulator undulator;
    undulator.k_ = getK();
    undulator.lu_ = getLambda();
    undulator.length_ = getNumPeriods();
    undulator.theta_ = getAngle() * Physics::pi / 180.0;
    lFringe_m = 2 * undulator.lu_;  // Default fringe field length is 2*lu.
    undulator.dist_ = lFringe_m - itsBunch->get_maxExtent()[2];  // Bunch-head to undulator distance.
    std::vector<MITHRA::Undulator> undulators;
    undulators.push_back(undulator);
    msg << "Undulator parameters have been transferred to the full-wave solver." << endl;

    MITHRA::Mesh mesh;
    mesh.initialize();
    mesh.lengthScale_ = 1.0;  // OPAL uses metres.
    mesh.timeScale_ = 1.0;  // OPAL uses seconds.
    mesh.meshCenter_ = MITHRA::FieldVector<double> (0.0);
    mesh.meshLength_ = getMeshLength();
    mesh.meshResolution_ = getMeshResolution();
    mesh.totalTime_ = getTotalTime();
    mesh.totalDist_ = lFringe_m + undulator.lu_ * undulator.length_;
    mesh.truncationOrder_ = getTruncationOrder();
    mesh.spaceCharge_ = true;
    mesh.optimizePosition_ = true;
    msg << "Mesh parameters have been transferred to the full-wave solver." << endl;

    MITHRA::Seed seed;
    std::vector<MITHRA::ExtField> externalFields;
    std::vector<MITHRA::FreeElectronLaser> FELs;
    
    // Get filename with desired output data.    
    if (!fname_m.empty()) {
        std::list<std::string> jobFile = MITHRA::read_file(fname_m.c_str());
        MITHRA::cleanJobFile(jobFile);
        MITHRA::ParseDarius parser (jobFile, mesh, bunch, seed, undulators, externalFields, FELs);
        parser.setJobParameters();
    }
    
    MITHRA::FdTdSC solver (mesh, bunch, seed, undulators, externalFields, FELs);
    
    // Transfer particles to full wave solver and destroy them from itsBunch.
    MITHRA::Charge charge;
    charge.q = itsBunch->getChargePerParticle() / (-Physics::q_e);
    for (unsigned int i = 0; i < localNum; ++i) {
        for (unsigned int d = 0; d < 3; ++d) {
            charge.rnp[d] = itsBunch->R[i][d];
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
    for (unsigned int i = 0; i < externalFields.size();  i++) {
        externalFields[i].show();
    }

    // Run the full-wave solver.
    timeval simulationStart;
    gettimeofday(&simulationStart, NULL);
    solve(solver, mesh, bunch, seed);

    // Get total computational time of the full wave simulation.
    timeval simulationEnd;
    gettimeofday(&simulationEnd, NULL);
    double deltaTime = ( simulationEnd.tv_usec - simulationStart.tv_usec ) / 1.0e6;
    deltaTime += ( simulationEnd.tv_sec - simulationStart.tv_sec );
    msg << "::: Total full wave simulation time [seconds] = " << deltaTime << endl;
    
    // Lorentz Transformation back to undulator local coordinates.
    double zMin = 1e100;
    for (auto iter = solver.chargeVectorn_.begin(); iter != solver.chargeVectorn_.end(); iter++) {
        zMin = std::min(zMin, iter->rnp[2]);
    }
    allreduce(&zMin, 1, std::less<double>());

    const double gb = solver.gamma_ * solver.beta_;
    const double factor = solver.gamma_ * (solver.beta_ * solver.c0_ * (solver.timeBunch_ + solver.dt_)) + lFringe_m;
    for (auto iter = solver.chargeVectorn_.begin(); iter != solver.chargeVectorn_.end(); iter++) {
        double dist = zMin - iter->rnp[2];
        // Lorentz transform.
        iter->rnp[2] = solver.gamma_ * iter->rnp[2] + factor;
        iter->gbnp[2] = solver.gamma_ * iter->gbnp[2] + gb * std::sqrt(1 + iter->gbnp.norm2());
        // Shift to bring all particles to same time in lab frame.
        double g = std::sqrt(1 + iter->gbnp.norm2());
        iter->rnp[0] += iter->gbnp[0] / g * dist * gb;
        iter->rnp[1] += iter->gbnp[1] / g * dist * gb;
        iter->rnp[2] += iter->gbnp[2] / g * dist * gb;
    }
    
    // Get total time ellapsed in laboratory frame.
    mesh.totalTime_ = solver.gamma_ * (solver.time_ + solver.beta_ / solver.c0_ * (zMin - bunch.zu_));

    
    // Return particles to itsBunch in local coordinates.
    msg << "Transferring particles back to OPAL bunch." << endl;
    itsBunch->create(solver.chargeVectorn_.size());
    const double dt = itsBunch->getdT();
    const unsigned int newLocalNum = itsBunch->getLocalNum();
    std::list<MITHRA::Charge>::iterator iter = solver.chargeVectorn_.begin();
    for (unsigned int i = 0; i < newLocalNum; ++i) {
        for (unsigned int d = 0; d < 3; ++d) {
            itsBunch->R[i][d] = iter->rnp[d];
            itsBunch->P[i][d] = iter->gbnp[d];
        }
        itsBunch->Q[i] = iter->q * (-Physics::q_e);
        itsBunch->dt[i] = dt;
        iter++;
    }
    itsBunch->setT(itsBunch->getT() + mesh.totalTime_);

    // Transfrom back to reference coordinate system.
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

    setHasBeenSimlated(true);
}

void Undulator::solve(MITHRA::FdTdSC & solver, MITHRA::Mesh& mesh, MITHRA::Bunch& bunch, MITHRA::Seed& seed) const {
    Inform msg("MITHRA FW solver ", *gmsg); 

    // Remark: This function is almost entirely copied from mithra/src/solver.cpp, but is adapted for OPAL.

    solver.initialize();

    
    // New total time try
    double Lu = solver.undulator_[0].lu_ * solver.undulator_[0].length_ / solver.gamma_;
    double zEnd = Lu + lFringe_m / solver.gamma_;
    double zMin = 1e100;
    double bz = 0;
    for (auto iter = solver.chargeVectorn_.begin(); iter != solver.chargeVectorn_.end(); iter++) {
        zMin = std::min(zMin, iter->rnp[2]);
        bz += iter->gbnp[2] / std::sqrt(1 + iter->gbnp.norm2());
    }
    allreduce(&zMin, 1, std::less<double>());
    allreduce(&bz, 1, std::plus<double>());
    unsigned int Nq = solver.chargeVectorn_.size();
    allreduce(&Nq, 1, std::plus<double>());
    bz /= Nq;

    if (mesh.totalTime_ == 0) {
        mesh.totalTime_ = 1 / (solver.c0_ * (bz + solver.beta_)) * (zEnd - solver.beta_ * solver.c0_ * solver.dt_ - zMin + bz / solver.beta_* Lu);
        msg << "Total time of the full wave simulation has been set to " << mesh.totalTime_ * solver.gamma_ << endl;
    }
    

    timeval simulationStart;
    gettimeofday(&simulationStart, NULL);

    msg << std::fixed << std::setprecision(3);
    msg << "-> Run the time domain simulation ..." << endl;
    double percentTime = 0.0;
    
    while (solver.time_ < mesh.totalTime_) {
        // Advance the fields for one time step using the FDTD algorithm.
        solver.fieldUpdate();

        /* If sampling of the field is enabled and the rhythm for sampling is achieved. Sample the
         * field at the given position and save them into the file.	*/
        if ( seed.sampling_ && fmod(solver.time_, seed.samplingRhythm_) < mesh.timeStep_ && solver.time_ > 0.0 )
            solver.fieldSample();

        /* If visualization of the field is enabled and the rhythm for visualization is achieved,
         * visualize the fields and save the vtk data in the given file name. */
        for (unsigned int i = 0; i < seed.vtk_.size(); i++) {
            if ( seed.vtk_[i].sample_ && fmod(solver.time_, seed.vtk_[i].rhythm_) < mesh.timeStep_ && solver.time_ > 0.0 ) {
                if ( seed.vtk_[i].type_ == MITHRA::ALLDOMAIN )
                    solver.fieldVisualizeAllDomain(i);
                else if (seed.vtk_[i].type_ == MITHRA::INPLANE)
                    solver.fieldVisualizeInPlane(i);
	        }
        }

        /* If profiling of the field is enabled and the time for profiling is achieved, write the
         * field profile and save the data in the given file name. */
        if (seed.profile_) {
            for (unsigned int i = 0; i < seed.profileTime_.size(); i++)
                if ( solver.time_ - seed.profileTime_[i] < mesh.timeStep_ && solver.time_ > seed.profileTime_[i] )
                    solver.fieldProfile();
            if ( fmod(solver.time_, seed.profileRhythm_) < mesh.timeStep_ && solver.time_ > 0.0 && seed.profileRhythm_ != 0 )
                solver.fieldProfile();
        }

        // Reset the charge and current values to zero.
        solver.currentReset();

        // Update the position and velocity parameters.
        for (auto iter = solver.chargeVectorn_.begin(); iter != solver.chargeVectorn_.end(); iter++) {
            iter->rnm  = iter->rnp;
            iter->gbnm = iter->gbnp;
        }

        /* Advance the particles till the time of the bunch properties reaches the time instant of the
         * field. */
        for (double t = 0.0; t < solver.nUpdateBunch_; t += 1.0) {
            solver.bunchUpdate();
            solver.timeBunch_ += bunch.timeStep_;
            ++solver.nTimeBunch_;
        }

        // Get particles going through a monitor.
        solver.screenProfile();

        // Deposit current and charge density on grid.
        solver.currentUpdate();
        solver.currentCommunicate();

        /* If sampling of the bunch is enabled and the rhythm for sampling is achieved. Sample the
         * bunch and save them into the file. */
        if ( bunch.sampling_ && fmod(solver.time_, bunch.rhythm_) < mesh.timeStep_ && solver.time_ > 0.0 )
            solver.bunchSample();

        /* If visualization of the bunch is enabled and the rhythm for visualization is achieved,
         * visualize the bunch and save the vtk data in the given file name. */
        if ( bunch.bunchVTK_ && fmod(solver.time_, bunch.bunchVTKRhythm_) < mesh.timeStep_ && solver.time_ > 0.0 )
            solver.bunchVisualize();

        /* If profiling of the bunch is enabled and the time for profiling is achieved, write the bunch
         * profile and save the data in the given file name. */
        if (bunch.bunchProfile_ > 0) {
            for (unsigned int i = 0; i < (bunch.bunchProfileTime_).size(); i++)
                if ( solver.time_ - bunch.bunchProfileTime_[i] < mesh.timeStep_ && solver.time_ > bunch.bunchProfileTime_[i] )
                    solver.bunchProfile();
            if ( fmod(solver.time_, bunch.bunchProfileRhythm_) < mesh.timeStep_ && solver.time_ > 0.0 && bunch.bunchProfileRhythm_ != 0.0 )
                solver.bunchProfile();
        }

        /* If radiation power of the FEL output is enabled and the rhythm for sampling is achieved.
         * Sample the radiation power at the given position and save them into the file. */
        solver.powerSample();
        solver.powerVisualize();

        /* If radiation energy of the FEL output is enabled and the rhythm for sampling is achieved.
         * Sample the radiation energy at the given position and save them into the file. */
        solver.energySample();

        /* Shift the computed fields and the time points for the fields. */
        solver.fieldShift();

        solver.timem1_ += mesh.timeStep_;
        solver.time_   += mesh.timeStep_;
        solver.timep1_ += mesh.timeStep_;
        ++solver.nTime_;

        timeval simulationEnd;
        gettimeofday(&simulationEnd, NULL);
        double deltaTime  = ( simulationEnd.tv_usec - simulationStart.tv_usec ) / 1.0e6;
        deltaTime += ( simulationEnd.tv_sec - simulationStart.tv_sec );

        if ( solver.rank_ == 0 && solver.time_ / mesh.totalTime_ * 1000.0 > percentTime ) {
            msg << " Percentage of the simulation completed (%)      = "
                  << solver.time_ / mesh.totalTime_ * 100.0 << endl;
            msg << " Average calculation time for each time step (s) = "
                  << deltaTime / (double)(solver.nTime_) << endl;
            msg << " Estimated remaining time (min)                  = "
                  << (mesh.totalTime_ / solver.time_ - 1) * deltaTime / 60 << endl;
            percentTime += 1.0;
        }
    }

    solver.finalize();
}

#endif

//set the number of slices for map tracking
void Undulator::setNSlices(const std::size_t& nSlices) { 
    nSlices_m = nSlices;
}

//get the number of slices for map tracking
std::size_t Undulator::getNSlices() const {
    return nSlices_m;
}

void Undulator::finalise() {
}

bool Undulator::bends() const {
    return false;
}

void Undulator::getDimensions(double &zBegin, double &zEnd) const {
    zBegin = startField_m;
    zEnd = startField_m + getElementLength();
}

ElementBase::ElementType Undulator::getType() const {
    return UNDULATOR;
}

void Undulator::setK(double k) { k_m = k; }
double Undulator::getK() const { return k_m; }

void Undulator::setLambda(double lambda) { lambda_m = lambda; }
double Undulator::getLambda() const { return lambda_m; }

void Undulator::setNumPeriods(int np) { numPeriods_m = np; }
int Undulator::getNumPeriods() const { return numPeriods_m; }

void Undulator::setAngle(double theta) { angle_m = theta; }
double Undulator::getAngle() const { return angle_m; }

void Undulator::setFilename(const std::string&fname) { fname_m = fname; }
const std::string& Undulator::getFilename() const { return fname_m;}

void Undulator::setMeshLength(std::vector<double> ml) { meshLength_m = ml; }
std::vector<double> Undulator::getMeshLength() const { return meshLength_m; }

void Undulator::setMeshResolution(std::vector<double> mr) { meshResolution_m = mr; }
std::vector<double> Undulator::getMeshResolution() const { return meshResolution_m; }

void Undulator::setTruncationOrder(unsigned int trunOrder) { truncationOrder_m = trunOrder; } 
unsigned int Undulator::getTruncationOrder() const { return truncationOrder_m; }

void Undulator::setTotalTime(double tt) { totalTime_m = tt; }
double Undulator::getTotalTime() const { return totalTime_m; }

void Undulator::setDtBunch(double dtb) { dtBunch_m = dtb; }
double Undulator::getDtBunch() const { return dtBunch_m; }

void Undulator::setHasBeenSimlated(bool hbs) { hasBeenSimulated_m = hbs; }
bool Undulator::getHasBeenSimulated() const { return hasBeenSimulated_m; }
