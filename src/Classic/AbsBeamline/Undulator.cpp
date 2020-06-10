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

/** Include header files for Mithra full-wave solver.  */
#include "mithra/fieldvector.h"
#include "mithra/stdinclude.h"
#include "mithra/readdata.h"
#include "mithra/database.h"
#include "mithra/classes.h"
#include "mithra/datainput.h"
#include "mithra/readdata.h"
#include "mithra/solver.h"
#include "mithra/fdtdSC.h"

extern Inform *gmsg;

// Class Undulator
// ------------------------------------------------------------------------

Undulator::Undulator():
    Component(),
    nSlices_m(1),
    k_m(0.0),
    lambda_m(0.0),
    numPeriods_m(0.0),
    fname_m(""),
    meshLength_m(3, 0.0),
    meshResolution_m(3, 0.0),
    truncationOrder_m(2),
    totalTime_m(0.0),
    dtBunch_m(0.0)
{ }


Undulator::Undulator(const Undulator &right):
    Component(right),
    nSlices_m(right.nSlices_m),
    k_m(right.k_m),
    lambda_m(right.lambda_m),
    numPeriods_m(right.numPeriods_m),
    fname_m(right.fname_m),
    meshLength_m(right.meshLength_m),
    meshResolution_m(right.meshResolution_m),
    truncationOrder_m(right.truncationOrder_m),
    totalTime_m(right.totalTime_m),
    dtBunch_m(right.dtBunch_m)
{ }


Undulator::Undulator(const std::string &name):
    Component(name),
    nSlices_m(1),
    k_m(0.0),
    lambda_m(0.0),
    numPeriods_m(0.0),
    fname_m(""),
    meshLength_m(3, 0.0),
    meshResolution_m(3, 0.0),
    truncationOrder_m(2),
    totalTime_m(0.0),
    dtBunch_m(0.0)
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

void Undulator::apply(PartBunchBase<double, 3> *itsBunch, CoordinateSystemTrafo const& refToLocalCSTrafo) {
    // Get local coordinates w.r.t. undulator.
    const unsigned int localNum = itsBunch->getLocalNum();
    for (unsigned int i = 0; i < localNum; ++i) {
        itsBunch->R[i] = refToLocalCSTrafo.transformTo(itsBunch->R[i]);
        itsBunch->P[i] = refToLocalCSTrafo.rotateTo(itsBunch->P[i]);
    }

    itsBunch->calcBeamParameters();
    *gmsg << "Bunch before undulator in local coordinate system: " << endl;
    itsBunch->print(*gmsg);    

    *gmsg << " :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::" << endl;
    *gmsg << " MITHRA-2.0: Completely Numerical Calculation of Free Electron Laser Radiation" << endl;
    *gmsg << " Version 2.0, Copyright 2019, Arya Fallahi" << endl;
    *gmsg << " Written by Arya Fallahi, IT'IS Foundation, Zurich, Switzerland" << endl;
    *gmsg << " :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::" << endl;

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
    *gmsg << "Bunch parameters have been transferred to the full-wave solver." << endl;

    MITHRA::Undulator undulator;
    undulator.k_ = getK();
    undulator.lu_ = getLambda();
    undulator.length_ = getNumPeriods();
    double fringe = 2 * undulator.lu_;  // Default fringe field length.
    undulator.dist_ = fringe - itsBunch->get_maxExtent()[2];  // Bunch-head to undulator distance.
    std::vector<MITHRA::Undulator> undulators;
    undulators.push_back(undulator);
    *gmsg << "Undulator parameters have been transferred to the full-wave solver." << endl;

    MITHRA::Mesh mesh;
    mesh.initialize();
    mesh.lengthScale_ = 1.0;  // OPAL uses metres
    mesh.timeScale_ = 1.0;  // OPAL uses seconds
    mesh.meshCenter_ = MITHRA::FieldVector<double> (0.0);
    mesh.meshLength_ = getMeshLength();
    mesh.meshResolution_ = getMeshResolution();
    mesh.totalTime_ = getTotalTime();
    mesh.truncationOrder_ = getTruncationOrder();
    mesh.spaceCharge_ = true;
    mesh.optimizePosition_ = true;
    *gmsg << "Mesh parameters have been transferred to the full-wave solver." << endl;

    MITHRA::Seed seed;
    std::vector<MITHRA::ExtField> externalFields;
    std::vector<MITHRA::FreeElectronLaser> FELs;
    
    // Get filename with desired output data.
    std::list<std::string> jobFile = MITHRA::read_file((getFilename()).c_str());
    MITHRA::cleanJobFile(jobFile);
    MITHRA::ParseDarius parser (jobFile, mesh, bunch, seed, undulators, externalFields, FELs);
    parser.setJobParameters();
    
    MITHRA::FdTdSC   fdtdsc   (mesh, bunch, seed, undulators, externalFields, FELs);
    // Transfer particles to MITHRA full-wave solver.
    MITHRA::Charge charge;
    charge.q = itsBunch->getChargePerParticle() / (-Physics::q_e);
    for (unsigned int i = 0; i < localNum; ++i) {
        for (unsigned int d = 0; d < 3; ++d) {
            charge.rnp[d] = itsBunch->R[i][d];
            charge.gbnp[d] = itsBunch->P[i][d];
        }
        fdtdsc.chargeVectorn_.push_back(charge);
    }
    *gmsg << "Particles have been transferred to the full-wave solver." << endl;

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
    
    // Run the full-wave solver
    fdtdsc.solve();
}

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
