// ------------------------------------------------------------------------
// $RCSfile: ParallelCyclotronTracker.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.1 $initialLocalNum_m
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: ParallelCyclotronTracker
//   The class for tracking particles with 3D space charge in Cyclotrons and FFAs
//
// ------------------------------------------------------------------------
//
// $Date: 2007/10/17 04:00:08 $
// $Author: adelmann, yang, winklehner $
//
// ------------------------------------------------------------------------

#include "Algorithms/ParallelCyclotronTracker.h"

#include <fstream>
#include <iostream>
#include <limits>
#include <vector>

#include "AbstractObjects/Element.h"
#include "AbstractObjects/OpalData.h"

#include "AbsBeamline/CCollimator.h"
#include "AbsBeamline/Corrector.h"
#include "AbsBeamline/Cyclotron.h"
#include "AbsBeamline/CyclotronValley.h"
#include "AbsBeamline/Degrader.h"
#include "AbsBeamline/Diagnostic.h"
#include "AbsBeamline/Drift.h"
#include "AbsBeamline/ElementBase.h"
#include "AbsBeamline/Lambertson.h"
#include "AbsBeamline/Offset.h"
#include "AbsBeamline/Marker.h"
#include "AbsBeamline/Monitor.h"
#include "AbsBeamline/Multipole.h"
#include "AbsBeamline/MultipoleT.h"
#include "AbsBeamline/MultipoleTBase.h"
#include "AbsBeamline/MultipoleTStraight.h"
#include "AbsBeamline/MultipoleTCurvedConstRadius.h"
#include "AbsBeamline/MultipoleTCurvedVarRadius.h"
#include "AbsBeamline/PluginElement.h"
#include "AbsBeamline/Probe.h"
#include "AbsBeamline/RBend.h"
#include "AbsBeamline/Ring.h"
#include "AbsBeamline/RFCavity.h"
#include "AbsBeamline/RFQuadrupole.h"
#include "AbsBeamline/SBend.h"
#include "AbsBeamline/SBend3D.h"
#include "AbsBeamline/ScalingFFAMagnet.h"
#include "AbsBeamline/Separator.h"
#include "AbsBeamline/Septum.h"
#include "AbsBeamline/Solenoid.h"
#include "AbsBeamline/Stripper.h"
#include "AbsBeamline/VariableRFCavity.h"
#include "AbsBeamline/VariableRFCavityFringeField.h"

#include "Algorithms/Ctunes.h"
#include "Algorithms/PolynomialTimeDependence.h"

#include "Beamlines/Beamline.h"
#include "Beamlines/FlaggedBeamline.h"

#include "Elements/OpalBeamline.h"

#include "Physics/Physics.h"

#include "Utilities/OpalException.h"
#include "Utilities/Options.h"

#include "BasicActions/DumpFields.h"
#include "BasicActions/DumpEMFields.h"

#include "Structure/BoundaryGeometry.h"
#include "Structure/DataSink.h"
#include "Structure/H5PartWrapperForPC.h"

//FIXME Remove headers and dynamic_cast in readOneBunchFromFile
#include "Algorithms/PartBunch.h"
#ifdef ENABLE_AMR
    #include "Algorithms/AmrPartBunch.h"
#endif

using Physics::pi;
using Physics::q_e;


const double c_mmtns = Physics::c * 1.0e-6; // m/s --> mm/ns

Vector_t const ParallelCyclotronTracker::xaxis = Vector_t(1.0, 0.0, 0.0);
Vector_t const ParallelCyclotronTracker::yaxis = Vector_t(0.0, 1.0, 0.0);
Vector_t const ParallelCyclotronTracker::zaxis = Vector_t(0.0, 0.0, 1.0);

extern Inform *gmsg;

/**
 * Constructor ParallelCyclotronTracker
 *
 * @param beamline
 * @param bunch
 * @param ds
 * @param reference
 * @param revBeam
 * @param revTrack
 * @param maxSTEPS
 * @param timeIntegrator
 */
ParallelCyclotronTracker::ParallelCyclotronTracker(const Beamline &beamline,
                                                   PartBunchBase<double, 3> *bunch,
                                                   DataSink &ds,
                                                   const PartData &reference,
                                                   bool revBeam, bool revTrack,
                                                   int maxSTEPS, int timeIntegrator,
                                                   int numBunch)
    : Tracker(beamline, bunch, reference, revBeam, revTrack)
    , itsMBDump_m(new MultiBunchDump())
    , bgf_m(nullptr)
    , maxSteps_m(maxSTEPS)
    , numBunch_m(numBunch)
    , lastDumpedStep_m(0)
    , eta_m(0.01)
    , myNode_m(Ippl::myNode())
    , initialLocalNum_m(bunch->getLocalNum())
    , initialTotalNum_m(bunch->getTotalNum())
    , onebunch_m(OpalData::getInstance()->getInputBasename() + "-onebunch.h5")
    , opalRing_m(nullptr)
    , itsStepper_mp(nullptr)
{
    itsBeamline = dynamic_cast<Beamline *>(beamline.clone());
    itsDataSink = &ds;
    multiBunchMode_m = MB_MODE::NONE;

    IntegrationTimer_m = IpplTimings::getTimer("Integration");
    TransformTimer_m   = IpplTimings::getTimer("Frametransform");
    DumpTimer_m        = IpplTimings::getTimer("Dump");
    BinRepartTimer_m   = IpplTimings::getTimer("Binaryrepart");

    // FIXME Change track command
    if ( initialTotalNum_m == 1 ) {
        mode_m = MODE::SINGLE;
    } else if ( initialTotalNum_m == 2 ) {
        mode_m = MODE::SEO;
    } else if ( initialTotalNum_m > 2 ) {
        mode_m = MODE::BUNCH;
    } else
        mode_m = MODE::UNDEFINED;

    if ( timeIntegrator == 0 ) {
        stepper_m = stepper::INTEGRATOR::RK4;
    } else if ( timeIntegrator == 1) {
        stepper_m = stepper::INTEGRATOR::LF2;
    } else if ( timeIntegrator == 2) {
        stepper_m = stepper::INTEGRATOR::MTS;
    } else
        stepper_m = stepper::INTEGRATOR::UNDEFINED;
}

/**
 * Destructor ParallelCyclotronTracker
 *
 */
ParallelCyclotronTracker::~ParallelCyclotronTracker() {
    for(Component* component : myElements) {
        delete(component);
    }
    for(auto fd : FieldDimensions) {
        delete(fd);
    }
    delete itsBeamline;
    // delete opalRing_m;
}

/// set the working sub-mode for multi-bunch mode: "FORCE" or "AUTO"
void ParallelCyclotronTracker::setMultiBunchMode(const std::string& mbmode)
{
    if ( mbmode.compare("FORCE") == 0 ) {
        *gmsg << "FORCE mode: The multi bunches will be injected consecutively" << endl
              << "            after each revolution, until get \"TURNS\" bunches." << endl;
        multiBunchMode_m = MB_MODE::FORCE;
    } else if ( mbmode.compare("AUTO") == 0 ) {
        *gmsg << "AUTO mode: The multi bunches will be injected only when the" << endl
              << "           distance between two neighboring bunches is below" << endl
              << "           the limitation. The control parameter is set to "
              << CoeffDBunches_m << endl;
        multiBunchMode_m = MB_MODE::AUTO;
    } else if ( mbmode.compare("NONE") == 0 )
        multiBunchMode_m = MB_MODE::NONE;
    else
        throw OpalException("ParallelCyclotronTracker::setMultiBunchMode()",
                            "MBMODE name \"" + mbmode + "\" unknown.");
}

void ParallelCyclotronTracker::setMultiBunchBinning(std::string binning) {

    binning = Util::toUpper(binning);

    if ( binning.compare("BUNCH") == 0 ) {
        *gmsg << "Use 'BUNCH' injection for binnning." << endl;
        binningType_m = MB_BINNING::BUNCH;
    } else if ( binning.compare("GAMMA") == 0 ) {
        *gmsg << "Use 'GAMMA' for binning." << endl;
        binningType_m = MB_BINNING::GAMMA;
    } else {
        throw OpalException("ParallelCyclotronTracker::setMultiBunchBinning()",
                            "MB_BINNING name \"" + binning + "\" unknown.");
    }
}

/**
 *
 *
 * @param fn Base file name
 */
void ParallelCyclotronTracker::bgf_main_collision_test() {
    if(!bgf_m) return;

    Inform msg("bgf_main_collision_test ");

    /**
     *Here we check if a particle is outside the domain, flag it for deletion
     */

    Vector_t intecoords = 0.0;

    // This has to match the dT in the rk4 pusher
    double dtime = itsBunch_m->getdT() * getHarmonicNumber();

    int triId = 0;
    for(size_t i = 0; i < itsBunch_m->getLocalNum(); i++) {
        int res = bgf_m->partInside(itsBunch_m->R[i], itsBunch_m->P[i], dtime, itsBunch_m->PType[i], itsBunch_m->Q[i], intecoords, triId);
        //int res = bgf_m->partInside(itsBunch_m->R[i]*1.0e-3, itsBunch_m->P[i], dtime, itsBunch_m->PType[i], itsBunch_m->Q[i], intecoords, triId);
        if(res >= 0) {
            itsBunch_m->Bin[i] = -1;
        }
    }
}


/**
 *
 *
 * @param fn Base file name
 */
void ParallelCyclotronTracker::openFiles(size_t numFiles, std::string SfileName) {

    for (unsigned int i=0; i<numFiles; i++) {
        std::string SfileName2 = SfileName;
        if (i<=2) {
            SfileName2 += std::string("-Angle" + std::to_string(int(i)) + ".dat");
        }
        else {
            // for single Particle Mode, output after each turn, to define matched initial phase ellipse.
            SfileName2 += std::string("-afterEachTurn.dat");
        }

        outfTheta_m.emplace_back(new std::ofstream(SfileName2.c_str()));
        outfTheta_m.back()->precision(8);
        outfTheta_m.back()->setf(std::ios::scientific, std::ios::floatfield);
        *outfTheta_m.back() << "# r [mm]        beta_r*gamma       "
                            << "theta [deg]     beta_theta*gamma        "
                            << "z [mm]          beta_z*gamma" << std::endl;
    }
}

/**
 * Close all files related to
 * special output in the Cyclotron
 * mode.
 */
void ParallelCyclotronTracker::closeFiles() {
    for (auto & file : outfTheta_m) {
        file->close();
    }
}

/**
 *
 * @param ring
 */
void ParallelCyclotronTracker::visitRing(const Ring &ring) {

    *gmsg << "* ----------------------------- Adding Ring ------------------------------ *" << endl;

    delete opalRing_m;

    opalRing_m = dynamic_cast<Ring*>(ring.clone());

    myElements.push_back(opalRing_m);

    opalRing_m->initialise(itsBunch_m);

    referenceR = opalRing_m->getBeamRInit();
    referencePr = opalRing_m->getBeamPRInit();
    referenceTheta = opalRing_m->getBeamPhiInit();

    if(referenceTheta <= -180.0 || referenceTheta > 180.0) {
        throw OpalException("Error in ParallelCyclotronTracker::visitRing",
                            "PHIINIT is out of [-180, 180)!");
    }

    referenceZ = 0.0;
    referencePz = 0.0;

    referencePtot = itsReference.getGamma() * itsReference.getBeta();
    referencePt = sqrt(referencePtot * referencePtot - referencePr * referencePr);

    if(referencePtot < 0.0)
        referencePt *= -1.0;

    sinRefTheta_m = sin(referenceTheta * Physics::deg2rad);
    cosRefTheta_m = cos(referenceTheta * Physics::deg2rad);

    double BcParameter[8] = {}; // zero initialise array

    buildupFieldList(BcParameter, ElementBase::RING, opalRing_m);

    // Finally print some diagnostic
    *gmsg << "* Initial beam radius = " << referenceR << " [mm] " << endl;
    *gmsg << "* Initial gamma = " << itsReference.getGamma() << endl;
    *gmsg << "* Initial beta = " << itsReference.getBeta() << endl;
    *gmsg << "* Total reference momentum   = " << referencePtot * 1000.0
          << " [MCU]" << endl;
    *gmsg << "* Reference azimuthal momentum  = " << referencePt * 1000.0
          << " [MCU]" << endl;
    *gmsg << "* Reference radial momentum     = " << referencePr * 1000.0
          << " [MCU]" << endl;
    *gmsg << "* " << opalRing_m->getSymmetry() << " fold field symmetry "
          << endl;
    *gmsg << "* Harmonic number h= " << opalRing_m->getHarmonicNumber() << " "
          << endl;
}

/**
 *
 *
 * @param cycl
 */
void ParallelCyclotronTracker::visitCyclotron(const Cyclotron &cycl) {

    *gmsg << "* -------------------------- Adding Cyclotron ---------------------------- *" << endl;

    Cyclotron *elptr = dynamic_cast<Cyclotron *>(cycl.clone());
    myElements.push_back(elptr);

    // Is this a Spiral Inflector Simulation? If yes, we'll give the user some
    // useful information
    spiral_flag = elptr->getSpiralFlag();

    if(spiral_flag) {

        *gmsg << endl << "* This is a Spiral Inflector Simulation! This means the following:" << endl;
        *gmsg         << "* 1.) It is up to the user to provide appropriate geometry, electric and magnetic fields!" << endl;
        *gmsg         << "*     (Use BANDRF type cyclotron and use RFMAPFN to load both magnetic" << endl;
        *gmsg         << "*     and electric fields, setting SUPERPOSE to an array of TRUE values.)" << endl;
        *gmsg         << "* 2.) For high currents it is strongly recommended to use the SAAMG fieldsolver," << endl;
        *gmsg         << "*     FFT does not give the correct results (boundary conditions are missing)." << endl;
        *gmsg         << "* 3.) The whole geometry will be meshed and used for the fieldsolver." << endl;
        *gmsg         << "*     There will be no transformations of the bunch into a local frame und consequently," << endl;
        *gmsg         << "*     the problem will be treated non-relativistically!" << endl;
        *gmsg         << "*     (This is not an issue for spiral inflectors as they are typically < 100 keV/amu.)" << endl;
        *gmsg << endl << "* Note: For now, multi-bunch mode (MBM) needs to be de-activated for spiral inflector" << endl;
        *gmsg         << "* and space charge needs to be solved every time-step. numBunch_m and scSolveFreq are reset." << endl;
        numBunch_m = 1;

    }

    // Fresh run (no restart):
    if(!OpalData::getInstance()->inRestartRun()) {

        // Get reference values from cyclotron element
        // For now, these are still stored in mm. should be the only ones. -DW
        referenceR     = elptr->getRinit();
        referenceTheta = elptr->getPHIinit();
        referenceZ     = elptr->getZinit();
        referencePr    = elptr->getPRinit();
        referencePz    = elptr->getPZinit();

        if(referenceTheta <= -180.0 || referenceTheta > 180.0) {
            throw OpalException("Error in ParallelCyclotronTracker::visitCyclotron",
                                "PHIINIT is out of [-180, 180)!");
        }

        referencePtot =  itsReference.getGamma() * itsReference.getBeta();

        // Calculate reference azimuthal (tangential) momentum from total-, z- and radial momentum:
        float insqrt = referencePtot * referencePtot - \
            referencePr * referencePr - referencePz * referencePz;

        if(insqrt < 0) {

            if(insqrt > -1.0e-10) {

                referencePt = 0.0;

            } else {

                throw OpalException("Error in ParallelCyclotronTracker::visitCyclotron",
                                    "Pt imaginary!");
            }

        } else {

            referencePt = sqrt(insqrt);
        }

        if(referencePtot < 0.0)
            referencePt *= -1.0;
        // End calculate referencePt

        // Restart a run:
    } else {

        // If the user wants to save the restarted run in local frame,
        // make sure the previous h5 file was local too
      if (Options::psDumpFrame != Options::GLOBAL) {
        if (!previousH5Local) {
                throw OpalException("Error in ParallelCyclotronTracker::visitCyclotron",
                                    "You are trying a local restart from a global h5 file!");
        }
            // Else, if the user wants to save the restarted run in global frame,
            // make sure the previous h5 file was global too
      } else {
        if (previousH5Local) {
                throw OpalException("Error in ParallelCyclotronTracker::visitCyclotron",
                                    "You are trying a global restart from a local h5 file!");
            }
        }

        // Adjust some of the reference variables from the h5 file
        referencePhi *= Physics::deg2rad;
        referencePsi *= Physics::deg2rad;
        referencePtot = bega;
        if(referenceTheta <= -180.0 || referenceTheta > 180.0) {
            throw OpalException("Error in ParallelCyclotronTracker::visitCyclotron",
                                "PHIINIT is out of [-180, 180)!");
        }
    }

    sinRefTheta_m = sin(referenceTheta * Physics::deg2rad);
    cosRefTheta_m = cos(referenceTheta * Physics::deg2rad);

    *gmsg << endl;
    *gmsg << "* Bunch global starting position:" << endl;
    *gmsg << "* RINIT = " << referenceR  << " [mm]" << endl;
    *gmsg << "* PHIINIT = " << referenceTheta << " [deg]" << endl;
    *gmsg << "* ZINIT = " << referenceZ << " [mm]" << endl;
    *gmsg << endl;
    *gmsg << "* Bunch global starting momenta:" << endl;
    *gmsg << "* Initial gamma = " << itsReference.getGamma() << endl;
    *gmsg << "* Initial beta = " << itsReference.getBeta() << endl;
    *gmsg << "* Reference total momentum (beta * gamma) = " << referencePtot * 1000.0 << " [MCU]" << endl;
    *gmsg << "* Reference azimuthal momentum (Pt) = " << referencePt * 1000.0 << " [MCU]" << endl;
    *gmsg << "* Reference radial momentum (Pr) = " << referencePr * 1000.0 << " [MCU]" << endl;
    *gmsg << "* Reference axial momentum (Pz) = " << referencePz * 1000.0 << " [MCU]" << endl;
    *gmsg << endl;

    double sym = elptr->getSymmetry();
    *gmsg << "* " << sym << "-fold field symmetry " << endl;

    // ckr: this just returned the default value as defined in Component.h
    // double rff = elptr->getRfFrequ();
    // *gmsg << "* Rf frequency= " << rff << " [MHz]" << endl;

    std::string fmfn = elptr->getFieldMapFN();
    *gmsg << "* Field map file name = " << fmfn << " " << endl;

    std::string type = elptr->getCyclotronType();
    *gmsg << "* Type of cyclotron = " << type << " " << endl;

    double rmin = elptr->getMinR();
    double rmax = elptr->getMaxR();
    *gmsg << "* Radial aperture = " << rmin << " ... " << rmax<<" [m] "<< endl;

    double zmin = elptr->getMinZ();
    double zmax = elptr->getMaxZ();
    *gmsg << "* Vertical aperture = " << zmin << " ... " << zmax<<" [m]"<< endl;

    double h = elptr->getCyclHarm();
    *gmsg << "* Number of trimcoils = " << elptr->getNumberOfTrimcoils() << endl;
    *gmsg << "* Harmonic number h = " << h << " " << endl;

    /**
     * To ease the initialise() function, set a integral parameter fieldflag internally.
     * Its value is  by the option "TYPE" of the element  "CYCLOTRON"
     * fieldflag = 1, readin PSI format measured field file (default)
     * fieldflag = 2, readin carbon cyclotron field file created by Jianjun Yang, TYPE=CARBONCYCL
     * fieldflag = 3, readin ANSYS format file for CYCIAE-100 created by Jianjun Yang, TYPE=CYCIAE
     * fieldflag = 4, readin AVFEQ format file for Riken cyclotrons
     * fieldflag = 5, readin FFA format file for MSU/FNAL FFA
     * fieldflag = 6, readin both median plane B field map and 3D E field map of RF cavity for compact cyclotron
     * fieldflag = 7, read in fields for Daniel's synchrocyclotron simulations
     */
    int fieldflag = elptr->getFieldFlag(type);

    // Read in cyclotron field maps (midplane + 3D fields if desired).
    elptr->initialise(itsBunch_m, fieldflag, elptr->getBScale());

    double BcParameter[8] = {};

    BcParameter[0] = 0.001 * elptr->getRmin();
    BcParameter[1] = 0.001 * elptr->getRmax();

    // Store inner radius and outer radius of cyclotron field map in the list
    buildupFieldList(BcParameter, ElementBase::CYCLOTRON, elptr);
}

/**
 * Not implemented and most probable never used
 *
 */
void ParallelCyclotronTracker::visitBeamBeam(const BeamBeam &) {
    *gmsg << "In BeamBeam tracker is missing " << endl;
}

/**
 *
 *
 * @param coll
 */
void ParallelCyclotronTracker::visitCCollimator(const CCollimator &coll) {

    *gmsg << "* --------- Collimator -----------------------------" << endl;

    CCollimator* elptr = dynamic_cast<CCollimator *>(coll.clone());
    myElements.push_back(elptr);

    double xstart = elptr->getXStart();
    *gmsg << "* Xstart= " << xstart << " [mm]" << endl;

    double xend = elptr->getXEnd();
    *gmsg << "* Xend= " << xend << " [mm]" << endl;

    double ystart = elptr->getYStart();
    *gmsg << "* Ystart= " << ystart << " [mm]" << endl;

    double yend = elptr->getYEnd();
    *gmsg << "* Yend= " << yend << " [mm]" << endl;

    double zstart = elptr->getZStart();
    *gmsg << "* Zstart= " << zstart << " [mm]" << endl;

    double zend = elptr->getZEnd();
    *gmsg << "* Zend= " << zend << " [mm]" << endl;

    double width = elptr->getWidth();
    *gmsg << "* Width= " << width << " [mm]" << endl;

    elptr->initialise(itsBunch_m);

    double BcParameter[8] = {};

    BcParameter[0] = 0.001 * xstart ;
    BcParameter[1] = 0.001 * xend;
    BcParameter[2] = 0.001 * ystart ;
    BcParameter[3] = 0.001 * yend;
    BcParameter[4] = 0.001 * width ;

    buildupFieldList(BcParameter, ElementBase::CCOLLIMATOR, elptr);
}

/**
 *
 *
 * @param corr
 */
void ParallelCyclotronTracker::visitCorrector(const Corrector &corr) {
    *gmsg << "In Corrector; L= " << corr.getElementLength() << endl;
    myElements.push_back(dynamic_cast<Corrector *>(corr.clone()));
}

/**
 *
 *
 * @param degrader
 */
void ParallelCyclotronTracker::visitDegrader(const Degrader &deg) {
    *gmsg << "In Degrader; L= " << deg.getElementLength() << endl;
    myElements.push_back(dynamic_cast<Degrader *>(deg.clone()));

}


/**
 *
 *
 * @param diag
 */
void ParallelCyclotronTracker::visitDiagnostic(const Diagnostic &diag) {
    *gmsg << "In Diagnostic; L= " << diag.getElementLength() << endl;
    myElements.push_back(dynamic_cast<Diagnostic *>(diag.clone()));
}

/**
 *
 *
 * @param drift
 */
void ParallelCyclotronTracker::visitDrift(const Drift &drift) {
    *gmsg << "In drift L= " << drift.getElementLength() << endl;
    myElements.push_back(dynamic_cast<Drift *>(drift.clone()));
}

/**
 *
 *
 *  @param
 */
void ParallelCyclotronTracker::visitFlexibleCollimator(const FlexibleCollimator &) {

}

/**
 *
 *
 * @param lamb
 */
void ParallelCyclotronTracker::visitLambertson(const Lambertson &lamb) {
    *gmsg << "In Lambertson; L= " << lamb.getElementLength() << endl;
    myElements.push_back(dynamic_cast<Lambertson *>(lamb.clone()));
}

void ParallelCyclotronTracker::visitOffset(const Offset & off) {
    if (opalRing_m == NULL)
        throw OpalException(
                            "ParallelCylcotronTracker::visitOffset",
                            "Attempt to place an offset when Ring not defined");
    Offset* offNonConst = const_cast<Offset*>(&off);
    offNonConst->updateGeometry(opalRing_m->getNextPosition(),
                                opalRing_m->getNextNormal());
    opalRing_m->appendElement(off);
}

/**
 *
 *
 * @param marker
 */
void ParallelCyclotronTracker::visitMarker(const Marker &marker) {
    //   *gmsg << "In Marker; L= " << marker.getElementLength() << endl;
    myElements.push_back(dynamic_cast<Marker *>(marker.clone()));
    // Do nothing.
}

/**
 *
 *
 * @param corr
 */
void ParallelCyclotronTracker::visitMonitor(const Monitor &corr) {
    //   *gmsg << "In Monitor; L= " << corr.getElementLength() << endl;
    myElements.push_back(dynamic_cast<Monitor *>(corr.clone()));
    //   applyDrift(flip_s * corr.getElementLength());
}

/**
 *
 *
 * @param mult
 */
void ParallelCyclotronTracker::visitMultipole(const Multipole &mult) {
    *gmsg << "In Multipole; L= " << mult.getElementLength() << " however the element is missing " << endl;
    myElements.push_back(dynamic_cast<Multipole *>(mult.clone()));
}

/**
 *
 *
 * @param multT
 */
void ParallelCyclotronTracker::visitMultipoleT(const MultipoleT &multT) {
    *gmsg << "Adding MultipoleT" << endl;
    if (opalRing_m != NULL) {
        opalRing_m->appendElement(multT);
    } else {
        throw OpalException("ParallelCyclotronTracker::visitMultipoleT",
                            "Need to define a RINGDEFINITION to use MultipoleT element");
    }
    myElements.push_back(dynamic_cast<MultipoleT *>(multT.clone()));
}

/**
 *
 *
 * @param multTstraight
 */
void ParallelCyclotronTracker::visitMultipoleTStraight(const MultipoleTStraight &multTstraight) {
    *gmsg << "Adding MultipoleTStraight" << endl;
    if (opalRing_m != NULL) {
        opalRing_m->appendElement(multTstraight);
    } else {
        throw OpalException("ParallelCyclotronTracker::visitMultipoleTStraight",
                            "Need to define a RINGDEFINITION to use MultipoleTStraight element");
    }
    myElements.push_back(dynamic_cast<MultipoleTStraight *>(multTstraight.clone()));
}

/**
 *
 *
 * @param multTccurv
 */
void ParallelCyclotronTracker::visitMultipoleTCurvedConstRadius(const MultipoleTCurvedConstRadius &multTccurv) {
    *gmsg << "Adding MultipoleTCurvedConstRadius" << endl;
    if (opalRing_m != NULL) {
        opalRing_m->appendElement(multTccurv);
    } else {
        throw OpalException("ParallelCyclotronTracker::visitMultipoleTCurvedConstRadius",
                            "Need to define a RINGDEFINITION to use MultipoleTCurvedConstRadius element");
    }
    myElements.push_back(dynamic_cast<MultipoleTCurvedConstRadius *>(multTccurv.clone()));
}

/**
 *
 *
 * @param multTvcurv
 */
void ParallelCyclotronTracker::visitMultipoleTCurvedVarRadius(const MultipoleTCurvedVarRadius &multTvcurv) {
    *gmsg << "Adding MultipoleTCurvedVarRadius" << endl;
    if (opalRing_m != NULL) {
        opalRing_m->appendElement(multTvcurv);
    } else {
        throw OpalException("ParallelCyclotronTracker::visitMultipoleTCurvedVarRadius",
                            "Need to define a RINGDEFINITION to use MultipoleTCurvedVarRadius element");
    }
    myElements.push_back(dynamic_cast<MultipoleTCurvedVarRadius *>(multTvcurv.clone()));
}

/**
 *
 *
 * @param prob
 */
void ParallelCyclotronTracker::visitProbe(const Probe &prob) {
    *gmsg << "* -----------  Probe -------------------------------" << endl;
    Probe *elptr = dynamic_cast<Probe *>(prob.clone());
    myElements.push_back(elptr);

    double xstart = elptr->getXStart();
    *gmsg << "XStart= " << xstart << " [mm]" << endl;

    double xend = elptr->getXEnd();
    *gmsg << "XEnd= " << xend << " [mm]" << endl;

    double ystart = elptr->getYStart();
    *gmsg << "YStart= " << ystart << " [mm]" << endl;

    double yend = elptr->getYEnd();
    *gmsg << "YEnd= " << yend << " [mm]" << endl;

    // initialise, do nothing
    elptr->initialise(itsBunch_m);

    double BcParameter[8] = {};

    BcParameter[0] = 0.001 * xstart ;
    BcParameter[1] = 0.001 * xend;
    BcParameter[2] = 0.001 * ystart ;
    BcParameter[3] = 0.001 * yend;
    BcParameter[4] = 0.001 ; // width

    // store probe parameters in the list
    buildupFieldList(BcParameter, ElementBase::PROBE, elptr);
}

/**
 *
 *
 * @param bend
 */
void ParallelCyclotronTracker::visitRBend(const RBend &bend) {
    *gmsg << "In RBend; L= " << bend.getElementLength() << " however the element is missing " << endl;
    myElements.push_back(dynamic_cast<RBend *>(bend.clone()));
}

void ParallelCyclotronTracker::visitSBend3D(const SBend3D &bend) {
    *gmsg << "Adding SBend3D" << endl;
    if (opalRing_m != NULL)
        opalRing_m->appendElement(bend);
    else
        throw OpalException("ParallelCyclotronTracker::visitSBend3D",
                            "Need to define a RINGDEFINITION to use SBend3D element");
}

void ParallelCyclotronTracker::visitScalingFFAMagnet(const ScalingFFAMagnet &bend) {
    *gmsg << "Adding ScalingFFAMagnet" << endl;
    if (opalRing_m != NULL) {
        opalRing_m->appendElement(bend);
    } else {
        throw OpalException("ParallelCyclotronTracker::visitScalingFFAMagnet",
                            "Need to define a RINGDEFINITION to use ScalingFFAMagnet element");
    }
}

void ParallelCyclotronTracker::visitVariableRFCavity(const VariableRFCavity &cav) {
    *gmsg << "Adding Variable RF Cavity" << endl;
    if (opalRing_m != NULL)
        opalRing_m->appendElement(cav);
    else
        throw OpalException("ParallelCyclotronTracker::visitVariableRFCavity",
                            "Need to define a RINGDEFINITION to use VariableRFCavity element");
}

void ParallelCyclotronTracker::visitVariableRFCavityFringeField
                                  (const VariableRFCavityFringeField &cav) {
    *gmsg << "Adding Variable RF Cavity with Fringe Field" << endl;
    if (opalRing_m != NULL)
        opalRing_m->appendElement(cav);
    else
        throw OpalException("ParallelCyclotronTracker::visitVariableRFCavityFringeField",
                            "Need to define a RINGDEFINITION to use VariableRFCavity element");
}

/**
 *
 *
 * @param as
 */
void ParallelCyclotronTracker::visitRFCavity(const RFCavity &as) {

    *gmsg << "* --------- RFCavity ------------------------------" << endl;

    RFCavity *elptr = dynamic_cast<RFCavity *>(as.clone());
    myElements.push_back(elptr);

    if((elptr->getComponentType() != "SINGLEGAP") && (elptr->getComponentType() != "DOUBLEGAP")) {
        *gmsg << (elptr->getComponentType()) << endl;
        throw OpalException("ParallelCyclotronTracker::visitRFCavity",
                            "The ParallelCyclotronTracker can only play with cyclotron type RF system currently ...");
    }

    double rmin = elptr->getRmin();
    *gmsg << "* Minimal radius of cavity= " << rmin << " [mm]" << endl;

    double rmax = elptr->getRmax();
    *gmsg << "* Maximal radius of cavity= " << rmax << " [mm]" << endl;

    double rff = elptr->getCycFrequency();
    *gmsg << "* RF frequency (2*pi*f)= " << rff << " [rad/s]" << endl;

    std::string fmfn = elptr->getFieldMapFN();
    *gmsg << "* RF Field map file name= " << fmfn << endl;

    double angle = elptr->getAzimuth();
    *gmsg << "* Cavity azimuth position= " << angle << " [deg] " << endl;

    double gap = elptr->getGapWidth();
    *gmsg << "* Cavity gap width= " << gap << " [mm] " << endl;

    double pdis = elptr->getPerpenDistance();
    *gmsg << "* Cavity Shift distance= " << pdis << " [mm] " << endl;

    double phi0 = elptr->getPhi0();
    *gmsg << "* Initial RF phase (t=0)= " << phi0 << " [deg] " << endl;

    /*
      Setup time dependence and in case of no
      timedependence use a polynom with  a_0 = 1 and a_k = 0, k = 1,2,3.
    */

    std::shared_ptr<AbstractTimeDependence> freq_atd = nullptr;
    std::shared_ptr<AbstractTimeDependence> ampl_atd = nullptr;
    std::shared_ptr<AbstractTimeDependence> phase_atd = nullptr;

    dvector_t  unityVec;
    unityVec.push_back(1.);
    unityVec.push_back(0.);
    unityVec.push_back(0.);
    unityVec.push_back(0.);

    if (elptr->getFrequencyModelName() != "") {
        freq_atd = AbstractTimeDependence::getTimeDependence(elptr->getFrequencyModelName());
        *gmsg << "* Variable frequency RF Model name " << elptr->getFrequencyModelName() << endl;
    }
    else
        freq_atd = std::shared_ptr<AbstractTimeDependence>(new PolynomialTimeDependence(unityVec));

    if (elptr->getAmplitudeModelName() != "") {
        ampl_atd = AbstractTimeDependence::getTimeDependence(elptr->getAmplitudeModelName());
        *gmsg << "* Variable amplitude RF Model name " << elptr->getAmplitudeModelName() << endl;
    }
    else
        ampl_atd = std::shared_ptr<AbstractTimeDependence>(new PolynomialTimeDependence(unityVec));

    if (elptr->getPhaseModelName() != "") {
        phase_atd = AbstractTimeDependence::getTimeDependence(elptr->getPhaseModelName());
        *gmsg << "* Variable phase RF Model name " << elptr->getPhaseModelName() << endl;
    }
    else
        phase_atd = std::shared_ptr<AbstractTimeDependence>(new PolynomialTimeDependence(unityVec));

    // read cavity voltage profile data from file.
    elptr->initialise(itsBunch_m, freq_atd, ampl_atd, phase_atd);

    double BcParameter[8] = {};

    BcParameter[0] = 0.001 * rmin;
    BcParameter[1] = 0.001 * rmax;
    BcParameter[2] = 0.001 * pdis;
    BcParameter[3] = angle;

    buildupFieldList(BcParameter, ElementBase::RFCAVITY, elptr);
}

/**
 *
 *
 * @param rfq
 */
void ParallelCyclotronTracker::visitRFQuadrupole(const RFQuadrupole &rfq) {
    *gmsg << "In RFQuadrupole; L = " << rfq.getElementLength() << " however the element is missing " << endl;
    myElements.push_back(dynamic_cast<RFQuadrupole *>(rfq.clone()));
}

/**
 *
 *
 * @param bend
 */
void ParallelCyclotronTracker::visitSBend(const SBend &bend) {
    *gmsg << "In SBend; L = " << bend.getElementLength() << " however the element is missing " << endl;
    myElements.push_back(dynamic_cast<SBend *>(bend.clone()));
}

/**
 *
 *
 * @param sep
 */
void ParallelCyclotronTracker::visitSeparator(const Separator &sep) {
    *gmsg << "In Seapator L= " << sep.getElementLength() << " however the element is missing " << endl;
    myElements.push_back(dynamic_cast<Separator *>(sep.clone()));
}

/**
 *
 *
 * @param sept
 */
void ParallelCyclotronTracker::visitSeptum(const Septum &sept) {

    *gmsg << endl << "* -----------------------------  Septum ------------------------------- *" << endl;

    Septum *elptr = dynamic_cast<Septum *>(sept.clone());
    myElements.push_back(elptr);

    double xstart = elptr->getXStart();
    *gmsg << "XStart = " << xstart << " [mm]" << endl;

    double xend = elptr->getXEnd();
    *gmsg << "XEnd = " << xend << " [mm]" << endl;

    double ystart = elptr->getYStart();
    *gmsg << "YStart = " << ystart << " [mm]" << endl;

    double yend = elptr->getYEnd();
    *gmsg << "YEnd = " << yend << " [mm]" << endl;

    double width = elptr->getWidth();
    *gmsg << "Width = " << width << " [mm]" << endl;


    // initialise, do nothing
    elptr->initialise(itsBunch_m);

    double BcParameter[8] = {};

    BcParameter[0] = 0.001 * xstart ;
    BcParameter[1] = 0.001 * xend;
    BcParameter[2] = 0.001 * ystart ;
    BcParameter[3] = 0.001 * yend;
    BcParameter[4] = 0.001 * width ;

    // store septum parameters in the list
    buildupFieldList(BcParameter, ElementBase::SEPTUM, elptr);
}

/**
 *
 *
 * @param solenoid
 */
void ParallelCyclotronTracker::visitSolenoid(const Solenoid &solenoid) {
    myElements.push_back(dynamic_cast<Solenoid *>(solenoid.clone()));
    Component *elptr = *(--myElements.end());
    if(!elptr->hasAttribute("ELEMEDGE")) {
        *gmsg << "Solenoid: no position of the element given!" << endl;
        return;
    }
}

/**
 *
 *
 * @param pplate
 */
void ParallelCyclotronTracker::visitParallelPlate(const ParallelPlate &pplate) {//do nothing

    //*gmsg << "ParallelPlate: not in use in ParallelCyclotronTracker!" << endl;

    //buildupFieldList(startField, endField, elptr);

}

/**
 *
 *
 * @param cv
 */
void ParallelCyclotronTracker::visitCyclotronValley(const CyclotronValley &cv) {
    // Do nothing here.
}

/**
 *
 *
 * @param stripper
 */

void ParallelCyclotronTracker::visitStripper(const Stripper &stripper) {

    *gmsg << "* ---------Stripper------------------------------" << endl;

    Stripper *elptr = dynamic_cast<Stripper *>(stripper.clone());
    myElements.push_back(elptr);

    double xstart = elptr->getXStart();
    *gmsg << "XStart= " << xstart << " [mm]" << endl;

    double xend = elptr->getXEnd();
    *gmsg << "XEnd= " << xend << " [mm]" << endl;

    double ystart = elptr->getYStart();
    *gmsg << "YStart= " << ystart << " [mm]" << endl;

    double yend = elptr->getYEnd();
    *gmsg << "YEnd= " << yend << " [mm]" << endl;

    double opcharge = elptr->getOPCharge();
    *gmsg << "Charge of outcoming particle = +e * " << opcharge << endl;

    double opmass = elptr->getOPMass();
    *gmsg << "* Mass of the outcoming particle = " << opmass << " [GeV/c^2]" << endl;

    elptr->initialise(itsBunch_m);

    double BcParameter[8] = {};

    BcParameter[0] = 0.001 * xstart ;
    BcParameter[1] = 0.001 * xend;
    BcParameter[2] = 0.001 * ystart ;
    BcParameter[3] = 0.001 * yend;
    BcParameter[4] = 0.001; // width
    BcParameter[5] = opcharge;
    BcParameter[6] = opmass;

    buildupFieldList(BcParameter, ElementBase::STRIPPER, elptr);
}

/**
 *
 *
 * @param BcParameter
 * @param ElementType
 * @param elptr
 */
void ParallelCyclotronTracker::buildupFieldList(double BcParameter[], ElementBase::ElementType elementType, Component *elptr) {
    beamline_list::iterator sindex;

    type_pair *localpair = new type_pair();
    localpair->first = elementType;

    for(int i = 0; i < 8; i++)
        *(((localpair->second).first) + i) = *(BcParameter + i);

    (localpair->second).second = elptr;

    // always put cyclotron as the first element in the list.
    if(elementType == ElementBase::RING || elementType == ElementBase::CYCLOTRON) {
        sindex = FieldDimensions.begin();
    } else {
        sindex = FieldDimensions.end();
    }
    FieldDimensions.insert(sindex, localpair);

}

/**
 *
 *
 * @param bl
 */
void ParallelCyclotronTracker::visitBeamline(const Beamline &bl) {
    const FlaggedBeamline* fbl = static_cast<const FlaggedBeamline*>(&bl);
    fbl->iterate(*this, false);//*dynamic_cast<BeamlineVisitor *>(this), false);
}

void ParallelCyclotronTracker::checkNumPart(std::string s) {
    int nlp = itsBunch_m->getLocalNum();
    int minnlp = 0;
    int maxnlp = 111111;
    reduce(nlp, minnlp, OpMinAssign());
    reduce(nlp, maxnlp, OpMaxAssign());
    *gmsg << s << " min local particle number " << minnlp << " max local particle number: " << maxnlp << endl;
}

/**
 *
 *
 */
void ParallelCyclotronTracker::execute() {

    /*
      Initialize common variables and structures
      for the integrators
    */

    step_m = 0;
    restartStep0_m = 0;

    // Record how many bunches have already been injected. ONLY FOR MPM
    BunchCount_m = itsBunch_m->getNumBunch();

    itsBeamline->accept(*this);
    if (opalRing_m != NULL)
        opalRing_m->lockRing();

    // Display the selected elements
    *gmsg << "* -------------------------------------" << endl;
    *gmsg << "* The selected Beam line elements are :" << endl;

    for(auto fd : FieldDimensions)
        *gmsg << "* -> " <<  ElementBase::getTypeString(fd->first) << endl;

    *gmsg << "* -------------------------------------" << endl;

    // Get BoundaryGeometry that is already initialized
    bgf_m = OpalData::getInstance()->getGlobalGeometry();

    // External field arrays for dumping
    for(int k = 0; k < 2; k++)
        FDext_m[k] = Vector_t(0.0, 0.0, 0.0);

    extE_m = Vector_t(0.0, 0.0, 0.0);
    extB_m = Vector_t(0.0, 0.0, 0.0);
    DumpFields::writeFields((((*FieldDimensions.begin())->second).second));
    DumpEMFields::writeFields((((*FieldDimensions.begin())->second).second));

    function_t func = std::bind(&ParallelCyclotronTracker::getFieldsAtPoint,
                                this,
                                std::placeholders::_1,
                                std::placeholders::_2,
                                std::placeholders::_3,
                                std::placeholders::_4);

    switch ( stepper_m )
    {
        case stepper::INTEGRATOR::RK4:
            *gmsg << "* 4th order Runge-Kutta integrator" << endl;
            itsStepper_mp.reset(new RK4<function_t>(func));
            break;
        case stepper::INTEGRATOR::LF2:
            *gmsg << "* 2nd order Leap-Frog integrator" << endl;
            itsStepper_mp.reset(new LF2<function_t>(func));
            break;
        case stepper::INTEGRATOR::MTS:
            *gmsg << "* Multiple time stepping (MTS) integrator" << endl;
            break;
        case stepper::INTEGRATOR::UNDEFINED:
        default:
            itsStepper_mp.reset(nullptr);
            throw OpalException("ParallelCyclotronTracker::execute",
                                "Invalid name of TIMEINTEGRATOR in Track command");
    }

    if ( stepper_m == stepper::INTEGRATOR::MTS)
        MtsTracker();
    else
        GenericTracker();


    *gmsg << "* ----------------------------------------------- *" << endl;
    *gmsg << "* Finalizing i.e. write data and close files :" << endl;
    for(auto fd : FieldDimensions) {
        ((fd->second).second)->finalise();
    }
    *gmsg << "* ----------------------------------------------- *" << endl;
}


void ParallelCyclotronTracker::MtsTracker() {
    /*
     * variable             unit        meaning
     * ------------------------------------------------
     * t                    [ns]        time
     * dt                   [ns]        time step
     * oldReferenceTheta    [rad]       azimuth angle
     * itsBunch_m->R        [m]         particle position
     *
     */

    double t = 0, dt = 0, oldReferenceTheta = 0;
    std::tie(t, dt, oldReferenceTheta) = initializeTracking_m();

    int const numSubsteps = std::max(Options::mtsSubsteps, 1);

    *gmsg << "MTS: Number of substeps per step is " << numSubsteps << endl;

    double const dt_inner = dt / double(numSubsteps);

    *gmsg << "MTS: The inner time step is therefore " << dt_inner << " [ns]" << endl;

//     int SteptoLastInj = itsBunch_m->getSteptoLastInj();

    bool flagTransition = false; // flag to determine when to transit from single-bunch to multi-bunches mode


    *gmsg << "* ---------------------------- Start tracking ----------------------------" << endl;

    if ( itsBunch_m->hasFieldSolver() )
        computeSpaceChargeFields_m();

    for(; (step_m < maxSteps_m) && (itsBunch_m->getTotalNum()>0); step_m++) {

        bool dumpEachTurn = false;

        if(step_m % Options::sptDumpFreq == 0) {
            singleParticleDump();
        }

        Ippl::Comm->barrier();

        // First half kick from space charge force
        if(itsBunch_m->hasFieldSolver()) {
            kick(0.5 * dt);
        }

        // Substeps for external field integration
        for(int n = 0; n < numSubsteps; ++n)
            borisExternalFields(dt_inner);


        // bunch injection
        // TODO: Where is correct location for this piece of code? Beginning/end of step? Before field solve?
        if(numBunch_m > 1)
            injectBunch_m(flagTransition);

        if ( itsBunch_m->hasFieldSolver() ) {
            computeSpaceChargeFields_m();
        } else {
            // if field solver is not available , only update bunch, to transfer particles between nodes if needed,
            // reset parameters such as LocalNum, initialTotalNum_m.
            // INFOMSG("No space charge Effects are included!"<<endl;);
            if((step_m % Options::repartFreq * 100) == 0) { //TODO: why * 100?
                Vector_t const meanP = calcMeanP();
                double const phi = calculateAngle(meanP(0), meanP(1)) - 0.5 * pi;
                Vector_t const meanR = calcMeanR();
                globalToLocal(itsBunch_m->R, phi, meanR);
                itsBunch_m->updateNumTotal();
                repartition();
                localToGlobal(itsBunch_m->R, phi, meanR);
            }
        }

        // Second half kick from space charge force
        if(itsBunch_m->hasFieldSolver())
            kick(0.5 * dt);

        // recalculate bingamma and reset the BinID for each particles according to its current gamma
        if((itsBunch_m->weHaveBins()) && BunchCount_m > 1) {
            if(step_m % Options::rebinFreq == 0) {
                updateParticleBins_m();
            }
        }

        // dump some data after one push in single particle tracking
        if ( mode_m == MODE::SINGLE ) {

            unsigned int i = 0;

            double temp_meanTheta = calculateAngle2(itsBunch_m->R[i](0),
                                                    itsBunch_m->R[i](1)); // [-pi, pi]

            dumpThetaEachTurn_m(t, itsBunch_m->R[i], itsBunch_m->P[i],
                                temp_meanTheta, dumpEachTurn);

            dumpAzimuthAngles_m(t, itsBunch_m->R[i], itsBunch_m->P[i],
                                oldReferenceTheta, temp_meanTheta);

            oldReferenceTheta = temp_meanTheta;
        } else if ( mode_m == MODE::BUNCH ) {
            // both for single bunch and multi-bunch
            // avoid dump at the first step
            // dumpEachTurn has not been changed in first push
            if((step_m > 10) && ((step_m + 1) % itsBunch_m->getStepsPerTurn()) == 0) {
                ++turnnumber_m;
                dumpEachTurn = true;

                *gmsg << endl;
                *gmsg << "*** Finished turn " << turnnumber_m - 1
                      << ", Total number of live particles: "
                      << itsBunch_m->getTotalNum() << endl;
            }

            // Recalculate bingamma and reset the BinID for each particles according to its current gamma
            if (itsBunch_m->weHaveBins() && BunchCount_m > 1 &&
                step_m % Options::rebinFreq == 0)
            {
                updateParticleBins_m();
            }
        }

        // printing + updating bunch parameters + updating t
        update_m(t, dt, dumpEachTurn);


    }

    // Some post-integration stuff
    *gmsg << endl;
    *gmsg << "* ---------------------------- DONE TRACKING PARTICLES -------------------------------- * " << endl;


    //FIXME
    dvector_t Ttime = dvector_t();
    dvector_t Tdeltr = dvector_t();
    dvector_t Tdeltz = dvector_t();
    ivector_t TturnNumber = ivector_t();

    finalizeTracking_m(Ttime, Tdeltr, Tdeltz, TturnNumber);
}


/**
 *
 *
 */
void ParallelCyclotronTracker::GenericTracker() {
    /*
     * variable             unit        meaning
     * ------------------------------------------------
     * t                    [ns]        time
     * dt                   [ns]        time step
     * oldReferenceTheta    [rad]       azimuth angle
     * itsBunch_m->R        [m]         particle position
     *
     */
    // Generic Tracker that has three modes defined by timeIntegrator_m:
    // 0 --- RK-4 (default)
    // 1 --- LF-2
    // (2 --- MTS ... not yet implemented in generic tracker)
    // numBunch_m determines the number of bunches in multibunch mode (MBM, 1 for OFF)
    // Total number of particles determines single particle mode (SPM, 1 particle) or
    // Static Equilibrium Orbit Mode (SEO, 2 particles)

    double t = 0, dt = 0, oldReferenceTheta = 0;
    std::tie(t, dt, oldReferenceTheta) = initializeTracking_m();

    // If initialTotalNum_m = 2, trigger SEO mode and prepare for transverse tuning calculation
    // Where is the IF here? -DW
    dvector_t Ttime, Tdeltr, Tdeltz;
    ivector_t TturnNumber;

    // Apply the plugin elements: probe, collimator, stripper, septum once before first step
    applyPluginElements(dt);

    // Destroy particles if they are marked as Bin = -1 in the plugin elements
    // or out of global aperture
    deleteParticle();

    /********************************
     *     Main integration loop    *
     ********************************/
    *gmsg << endl;
    *gmsg << "* --------------------------------- Start tracking ------------------------------------ *" << endl;

    for(; (step_m < maxSteps_m) && (itsBunch_m->getTotalNum()>0); step_m++) {

        bool dumpEachTurn = false;

        switch (mode_m)
        {
            case MODE::SEO:
            { // initialTotalNum_m == 2
                seoMode_m(t, dt, dumpEachTurn, Ttime, Tdeltr, Tdeltz, TturnNumber);
                break;
            }
            case MODE::SINGLE:
            { // initialTotalNum_m == 1
                singleMode_m(t, dt, dumpEachTurn, oldReferenceTheta);
                break;
            }
            case MODE::BUNCH:
            { // initialTotalNum_m > 2
                // Start Tracking for number of particles > 2 (i.e. not single and not SEO mode)
                bunchMode_m(t, dt, dumpEachTurn);
                break;
            }
            case MODE::UNDEFINED:
            default:
                throw OpalException("ParallelCyclotronTracker::GenericTracker()",
                                    "No such tracking mode.");
        }
        // Update bunch and some parameters and output some info
        update_m(t, dt, dumpEachTurn);

    } // end for: the integration is DONE after maxSteps_m steps or if all particles are lost!

    // Some post-integration stuff
    *gmsg << endl;
    *gmsg << "* ---------------------------- DONE TRACKING PARTICLES -------------------------------- * " << endl;

    finalizeTracking_m(Ttime, Tdeltr, Tdeltz, TturnNumber);
}

bool ParallelCyclotronTracker::getFieldsAtPoint(const double &t, const size_t &Pindex, Vector_t &Efield, Vector_t &Bfield) {

    bool outOfBound = this->computeExternalFields_m(Pindex, t, Efield, Bfield);

    // For runs without space charge effects, override this step to save time
    if(itsBunch_m->hasFieldSolver()) {

        // Don't do for reference particle
        if(itsBunch_m->ID[Pindex] != 0) {

            // add external Field and self space charge field
            Efield += itsBunch_m->Ef[Pindex];
            Bfield += itsBunch_m->Bf[Pindex];
        }
    }

    return outOfBound;
}


/**
 *
 *
 * @param Rold
 * @param Rnew
 * @param elptr
 * @param Dold
 *
 * @return
 */

bool ParallelCyclotronTracker::checkGapCross(Vector_t Rold, Vector_t Rnew,
                                             RFCavity * rfcavity, double &Dold)
{
    bool flag = false;
    double sinx = rfcavity->getSinAzimuth();
    double cosx = rfcavity->getCosAzimuth();
    // TODO: Presumably this is still in mm, so for now, change to m -DW
    double PerpenDistance = 0.001 * rfcavity->getPerpenDistance();
    double distNew = (Rnew[0] * sinx - Rnew[1] * cosx) - PerpenDistance;
    double distOld = (Rold[0] * sinx - Rold[1] * cosx) - PerpenDistance;
    if(distOld > 0.0 && distNew <= 0.0) flag = true;
    // This parameter is used correct cavity phase
    Dold = distOld;
    return flag;
}

bool ParallelCyclotronTracker::RFkick(RFCavity * rfcavity, const double t, const double dt, const int Pindex){
    // For OPAL 2.0: As long as the RFCavity is in mm, we have to change R to mm here -DW
    double radius = sqrt(pow(1000.0 * itsBunch_m->R[Pindex](0), 2.0) + pow(1000.0 * itsBunch_m->R[Pindex](1), 2.0)
                         - pow(rfcavity->getPerpenDistance() , 2.0));
    double rmin = rfcavity->getRmin();
    double rmax = rfcavity->getRmax();
    double nomalRadius = (radius - rmin) / (rmax - rmin);
    double tempP[3];
    if(nomalRadius <= 1.0 && nomalRadius >= 0.0) {

        for(int j = 0; j < 3; j++)
            tempP[j] = itsBunch_m->P[Pindex](j);  //[px,py,pz]  units: dimensionless

        // here evaluate voltage and conduct momenta kicking;
        rfcavity->getMomentaKick(nomalRadius, tempP, t, dt, itsBunch_m->ID[Pindex], itsBunch_m->getM(), itsBunch_m->getQ()); // t : ns

        for(int j = 0; j < 3; j++)
            itsBunch_m->P[Pindex](j) = tempP[j];
        return true;
    }
    return false;
}


struct adder : public std::unary_function<double, void> {
    adder() : sum(0) {}
    double sum;
    void operator()(double x) { sum += x; }
};

/**
 *
 *
 * @param t
 * @param r
 * @param z
 * @param lastTurn
 * @param nur
 * @param nuz
 *
 * @return
 */
bool ParallelCyclotronTracker::getTunes(dvector_t &t, dvector_t &r, dvector_t &z,
                                        int lastTurn, double &nur, double &nuz) {
    TUNE_class *tune;

    int Ndat = t.size();

    /*
      remove mean
    */
    double rsum =  for_each(r.begin(), r.end(), adder()).sum;

    for(int i = 0; i < Ndat; i++)
        r[i] -= rsum;

    double zsum =  for_each(z.begin(), z.end(), adder()).sum;

    for(int i = 0; i < Ndat; i++)
        z[i] -= zsum;
    double ti = *(t.begin());
    double tf = t[t.size()-1];
    double T = (tf - ti);

    t.clear();
    double dt = T / Ndat;
    double time = 0.0;
    for(int i = 0; i < Ndat; i++) {
        t.push_back(i);
        time += dt;
    }

    T = t[Ndat-1];

    *gmsg << endl;
    *gmsg << "* ************************************* nuR ******************************************* *" << endl;
    *gmsg << endl << "* ===> " << Ndat << " data points  Ti=" << ti << " Tf= " << tf << " -> T= " << T << endl;

    int nhis_lomb = 10;
    int  stat = 0;
    // book tune class
    tune = new TUNE_class();
    stat = tune->lombAnalysis(t, r, nhis_lomb, T / lastTurn);
    if(stat != 0)
        *gmsg << "* TUNE: Lomb analysis failed" << endl;
    *gmsg << "* ************************************************************************************* *" << endl;

    delete tune;
    tune = NULL;
    // FIXME: FixMe: need to come from the inputfile
    nhis_lomb = 100;

    if(zsum != 0.0) {

        *gmsg << endl;
        *gmsg << "* ************************************* nuZ ******************************************* *" << endl;
        *gmsg << endl << "* ===> " << Ndat << " data points  Ti=" << ti << " Tf= " << tf << " -> T= " << T << endl;

        // book tune class
        tune = new TUNE_class();
        stat = tune->lombAnalysis(t, z, nhis_lomb, T / lastTurn);
        if(stat != 0)
            *gmsg << "* TUNE: Lomb analysis failed" << endl;
        *gmsg << "* ************************************************************************************* *" << endl;

        delete tune;
        tune = NULL;
    }
    return true;
}

void ParallelCyclotronTracker::saveOneBunch() {

    *gmsg << endl;
    *gmsg << "* Store beam to H5 file for multibunch simulation ... ";

    Ppos_t coord, momentum;
    ParticleAttrib<double> mass, charge;
    ParticleAttrib<short> ptype;

    std::size_t localNum = itsBunch_m->getLocalNum();

    coord.create(localNum);
    coord = itsBunch_m->R;

    momentum.create(localNum);
    momentum = itsBunch_m->P;

    mass.create(localNum);
    mass = itsBunch_m->M;

    charge.create(localNum);
    charge = itsBunch_m->Q;

    ptype.create(localNum);
    ptype = itsBunch_m->PType;

    std::map<std::string, double> additionalAttributes = {
        std::make_pair("REFPR", 0.0),
        std::make_pair("REFPT", 0.0),
        std::make_pair("REFPZ", 0.0),
        std::make_pair("REFR", 0.0),
        std::make_pair("REFTHETA", 0.0),
        std::make_pair("REFZ", 0.0),
        std::make_pair("AZIMUTH", 0.0),
        std::make_pair("ELEVATION", 0.0),
        std::make_pair("B-ref_x",  0.0),
        std::make_pair("B-ref_z",  0.0),
        std::make_pair("B-ref_y",  0.0),
        std::make_pair("E-ref_x",  0.0),
        std::make_pair("E-ref_z",  0.0),
        std::make_pair("E-ref_y",  0.0),
        std::make_pair("B-head_x", 0.0),
        std::make_pair("B-head_z", 0.0),
        std::make_pair("B-head_y", 0.0),
        std::make_pair("E-head_x", 0.0),
        std::make_pair("E-head_z", 0.0),
        std::make_pair("E-head_y", 0.0),
        std::make_pair("B-tail_x", 0.0),
        std::make_pair("B-tail_z", 0.0),
        std::make_pair("B-tail_y", 0.0),
        std::make_pair("E-tail_x", 0.0),
        std::make_pair("E-tail_z", 0.0),
        std::make_pair("E-tail_y", 0.0)
    };

    H5PartWrapperForPC h5wrapper(onebunch_m, H5_O_WRONLY);
    h5wrapper.writeHeader();
    h5wrapper.writeStep(itsBunch_m, additionalAttributes);
    h5wrapper.close();

    *gmsg << "Done." << endl;
}


bool ParallelCyclotronTracker::readOneBunchFromFile(const size_t BinID) {

    *gmsg << endl;
    *gmsg << "* Read beam from H5 file for multibunch simulation ... ";

    std::size_t localNum = itsBunch_m->getLocalNum();

    /*
     * 2nd argument: 0  --> step zero is fine since the file has only this step
     * 3rd argument: "" --> onebunch_m is used
     * 4th argument: H5_O_RDONLY does not work with this constructor
     */
    H5PartWrapperForPC h5wrapper(onebunch_m, 0, "", H5_O_WRONLY);

    size_t numParticles = h5wrapper.getNumParticles();

    if ( numParticles == 0 ) {
        throw OpalException("ParallelCyclotronTracker::readOneBunchFromFile()",
                            "No particles in file " + onebunch_m + ".");
    }

    size_t numParticlesPerNode = numParticles / Ippl::getNodes();

    size_t firstParticle = numParticlesPerNode * Ippl::myNode();
    size_t lastParticle = firstParticle + numParticlesPerNode - 1;
    if (Ippl::myNode() == Ippl::getNodes() - 1)
        lastParticle = numParticles - 1;

    PAssert_LT(firstParticle, lastParticle +1);

    numParticles = lastParticle - firstParticle + 1;

    //FIXME
    std::unique_ptr<PartBunchBase<double, 3> > tmpBunch = 0;
#ifdef ENABLE_AMR
    AmrPartBunch* amrbunch_p = dynamic_cast<AmrPartBunch*>(itsBunch_m);
    if ( amrbunch_p != 0 ) {
        tmpBunch.reset(new AmrPartBunch(&itsReference,
                                        amrbunch_p->getAmrParticleBase()));
    } else
#endif
        tmpBunch.reset(new PartBunch(&itsReference));

    tmpBunch->create(numParticles);

    h5wrapper.readStep(tmpBunch.get(), firstParticle, lastParticle);
    h5wrapper.close();

    itsBunch_m->create(numParticles);

    for(size_t ii = 0; ii < numParticles; ++ ii, ++ localNum) {
        itsBunch_m->R[localNum] = tmpBunch->R[ii];
        itsBunch_m->P[localNum] = tmpBunch->P[ii];
        itsBunch_m->M[localNum] = tmpBunch->M[ii];
        itsBunch_m->Q[localNum] = tmpBunch->Q[ii];
        itsBunch_m->PType[localNum] = ParticleType::REGULAR;
        itsBunch_m->Bin[localNum] = BinID;
        itsBunch_m->bunchNum[localNum] = BunchCount_m - 1;
    }

    itsBunch_m->boundp();

    *gmsg << "Done." << endl;

    return true;
}


void ParallelCyclotronTracker::updateParticleBins_m() {
    switch ( binningType_m ) {
        case MB_BINNING::GAMMA:
            itsBunch_m->resetPartBinID2(eta_m);
            break;
        case MB_BINNING::BUNCH:
            itsBunch_m->resetPartBinBunch();
            break;
        default:
            itsBunch_m->resetPartBinID2(eta_m);
    }
}


double ParallelCyclotronTracker::getHarmonicNumber() const {
    if (opalRing_m != NULL)
        return opalRing_m->getHarmonicNumber();
    Cyclotron* elcycl = dynamic_cast<Cyclotron*>(((*FieldDimensions.begin())->second).second);
    if (elcycl != NULL)
        return elcycl->getCyclHarm();
    throw OpalException("ParallelCyclotronTracker::getHarmonicNumber()",
                        std::string("The first item in the FieldDimensions list does not ")
                        +std::string("seem to be an Ring or a Cyclotron element"));
}


Vector_t ParallelCyclotronTracker::calcMeanR() const {
    Vector_t meanR(0.0, 0.0, 0.0);

    for(unsigned int i = 0; i < itsBunch_m->getLocalNum(); ++i) {
        for(int d = 0; d < 3; ++d) {
            meanR(d) += itsBunch_m->R[i](d);
        }
    }

    reduce(meanR, meanR, OpAddAssign());
    return meanR / Vector_t(itsBunch_m->getTotalNum());
}

Vector_t ParallelCyclotronTracker::calcMeanP() const {
    Vector_t meanP(0.0, 0.0, 0.0);

    for(unsigned int i = 0; i < itsBunch_m->getLocalNum(); ++i) {
        for(int d = 0; d < 3; ++d) {
            meanP(d) += itsBunch_m->P[i](d);
        }
    }

    reduce(meanP, meanP, OpAddAssign());
    return meanP / Vector_t(itsBunch_m->getTotalNum());
}

void ParallelCyclotronTracker::repartition() {
    if((step_m % Options::repartFreq) == 0) {
        IpplTimings::startTimer(BinRepartTimer_m);
        itsBunch_m->do_binaryRepart();
        Ippl::Comm->barrier();
        IpplTimings::stopTimer(BinRepartTimer_m);
    }
}

void ParallelCyclotronTracker::globalToLocal(ParticleAttrib<Vector_t> & particleVectors,
                                             double phi, Vector_t const translationToGlobal) {
    particleVectors -= translationToGlobal;

    Tenzor<double, 3> const rotation( cos(phi), sin(phi), 0,
                                      -sin(phi), cos(phi), 0,
                                      0,        0, 1); // clockwise rotation

    for(unsigned int i = 0; i < itsBunch_m->getLocalNum(); ++i) {

        particleVectors[i] = dot(rotation, particleVectors[i]);
    }
}

void ParallelCyclotronTracker::localToGlobal(ParticleAttrib<Vector_t> & particleVectors,
                                             double phi, Vector_t const translationToGlobal) {

    Tenzor<double, 3> const rotation(cos(phi), -sin(phi), 0,
                                     sin(phi),  cos(phi), 0,
                                     0,         0, 1); // counter-clockwise rotation

    for(unsigned int i = 0; i < itsBunch_m->getLocalNum(); ++i) {

        particleVectors[i] = dot(rotation, particleVectors[i]);
    }

    particleVectors += translationToGlobal;
}


inline void ParallelCyclotronTracker::globalToLocal(ParticleAttrib<Vector_t> & particleVectors,
                                                    Quaternion_t const quaternion,
                                                    Vector_t const meanR) {

    // Translation from global to local
    particleVectors -= meanR;

    // Rotation to align P_mean with x-axis
    rotateWithQuaternion(particleVectors, quaternion);
}

inline void ParallelCyclotronTracker::localToGlobal(ParticleAttrib<Vector_t> & particleVectors,
                                                    Quaternion_t const quaternion,
                                                    Vector_t const meanR) {

    // Reverse the quaternion by multiplying the axis components (x,y,z) with -1
    Quaternion_t reverseQuaternion = quaternion * -1.0;
    reverseQuaternion(0) *= -1.0;

    // Rotation back to original P_mean direction
    rotateWithQuaternion(particleVectors, reverseQuaternion);

    // Translation from local to global
    particleVectors += meanR;
}

inline void ParallelCyclotronTracker::globalToLocal(ParticleAttrib<Vector_t> & particleVectors,
                                                    double const phi,
                                                    double const psi,
                                                    Vector_t const meanR) {

    //double const tolerance = 1.0e-4; // TODO What is a good angle threshold? -DW

    // Translation from global to local
    particleVectors -= meanR;

    // Rotation to align P_mean with x-axis
    rotateAroundZ(particleVectors, phi);

    // If theta is large enough (i.e. there is significant momentum in z direction)
    // rotate around x-axis next
    //if (fabs(psi) > tolerance)
    rotateAroundX(particleVectors, psi);
}

inline void ParallelCyclotronTracker::globalToLocal(Vector_t & myVector,
                                                    double const phi,
                                                    double const psi,
                                                    Vector_t const meanR) {

    //double const tolerance = 1.0e-4; // TODO What is a good angle threshold? -DW

    // Translation from global to local
    myVector -= meanR;

    // Rotation to align P_mean with x-axis
    rotateAroundZ(myVector, phi);

    // If theta is large enough (i.e. there is significant momentum in z direction)
    // rotate around x-axis next
    //if (fabs(psi) > tolerance)
    rotateAroundX(myVector, psi);
}

inline void ParallelCyclotronTracker::localToGlobal(ParticleAttrib<Vector_t> & particleVectors,
                                                    double const phi,
                                                    double const psi,
                                                    Vector_t const meanR) {

    //double const tolerance = 1.0e-4; // TODO What is a good angle threshold? -DW

    // If theta is large enough (i.e. there is significant momentum in z direction)
    // rotate back around x-axis next
    //if (fabs(psi) > tolerance)
    rotateAroundX(particleVectors, -psi);

    // Rotation to align P_mean with x-axis
    rotateAroundZ(particleVectors, -phi);

    // Translation from local to global
    particleVectors += meanR;
}

inline void ParallelCyclotronTracker::localToGlobal(Vector_t & myVector,
                                                    double const phi,
                                                    double const psi,
                                                    Vector_t const meanR) {

    //double const tolerance = 1.0e-4; // TODO What is a good angle threshold? -DW

    // If theta is large enough (i.e. there is significant momentum in z direction)
    // rotate back around x-axis next
    //if (fabs(psi) > tolerance)
    rotateAroundX(myVector, -psi);

    // Rotation to align P_mean with x-axis
    rotateAroundZ(myVector, -phi);

    // Translation from local to global
    myVector += meanR;
}

inline void ParallelCyclotronTracker::rotateWithQuaternion(ParticleAttrib<Vector_t> & particleVectors, Quaternion_t const quaternion) {

    Vector_t const quaternionVectorComponent = Vector_t(quaternion(1), quaternion(2), quaternion(3));
    double const quaternionScalarComponent = quaternion(0);

    for(unsigned int i = 0; i < itsBunch_m->getLocalNum(); ++i) {

        particleVectors[i] = 2.0f * dot(quaternionVectorComponent, particleVectors[i]) * quaternionVectorComponent +
            (quaternionScalarComponent * quaternionScalarComponent -
             dot(quaternionVectorComponent, quaternionVectorComponent)) * particleVectors[i] + 2.0f *
            quaternionScalarComponent * cross(quaternionVectorComponent, particleVectors[i]);
    }
}

inline void ParallelCyclotronTracker::normalizeQuaternion(Quaternion_t & quaternion){

    double tolerance = 1.0e-10;
    double length2 = dot(quaternion, quaternion);

    if (fabs(length2) > tolerance && fabs(length2 - 1.0f) > tolerance) {

        double length = sqrt(length2);
        quaternion /= length;
    }
}

inline void ParallelCyclotronTracker::normalizeVector(Vector_t & vector) {

    double tolerance = 1.0e-10;
    double length2 = dot(vector, vector);

    if (fabs(length2) > tolerance && fabs(length2 - 1.0f) > tolerance) {

        double length = sqrt(length2);
        vector /= length;
    }
}

inline void ParallelCyclotronTracker::rotateAroundZ(ParticleAttrib<Vector_t> & particleVectors, double const phi) {
    // Clockwise rotation of particles 'particleVectors' by 'phi' around Z axis

    Tenzor<double, 3> const rotation( cos(phi), sin(phi), 0,
                                      -sin(phi), cos(phi), 0,
                                      0,        0, 1);

    for(unsigned int i = 0; i < itsBunch_m->getLocalNum(); ++i) {

        particleVectors[i] = dot(rotation, particleVectors[i]);
    }
}

inline void ParallelCyclotronTracker::rotateAroundZ(Vector_t & myVector, double const phi) {
    // Clockwise rotation of single vector 'myVector' by 'phi' around Z axis

    Tenzor<double, 3> const rotation( cos(phi), sin(phi), 0,
                                      -sin(phi), cos(phi), 0,
                                      0,        0, 1);

    myVector = dot(rotation, myVector);
}

inline void ParallelCyclotronTracker::rotateAroundX(ParticleAttrib<Vector_t> & particleVectors, double const psi) {
    // Clockwise rotation of particles 'particleVectors' by 'psi' around X axis

    Tenzor<double, 3> const rotation(1,  0,          0,
                                     0,  cos(psi), sin(psi),
                                     0, -sin(psi), cos(psi));

    for(unsigned int i = 0; i < itsBunch_m->getLocalNum(); ++i) {

        particleVectors[i] = dot(rotation, particleVectors[i]);
    }
}

inline void ParallelCyclotronTracker::rotateAroundX(Vector_t & myVector, double const psi) {
    // Clockwise rotation of single vector 'myVector' by 'psi' around X axis

    Tenzor<double, 3> const rotation(1,  0,          0,
                                     0,  cos(psi), sin(psi),
                                     0, -sin(psi), cos(psi));

    myVector = dot(rotation, myVector);
}

inline void ParallelCyclotronTracker::getQuaternionTwoVectors(Vector_t u, Vector_t v, Quaternion_t & quaternion) {
    // four vector (w,x,y,z) of the quaternion of P_mean with the positive x-axis

    normalizeVector(u);
    normalizeVector(v);

    double k_cos_theta = dot(u, v);
    double k = sqrt(dot(u, u) * dot(v, v));
    double tolerance1 = 1.0e-5;
    double tolerance2 = 1.0e-8;
    Vector_t resultVectorComponent;

    if (fabs(k_cos_theta / k + 1.0) < tolerance1) {
        // u and v are almost exactly antiparallel so we need to do
        // 180 degree rotation around any vector orthogonal to u

        resultVectorComponent = cross(u, xaxis);

        // If by chance u is parallel to xaxis, use zaxis instead
        if (dot(resultVectorComponent, resultVectorComponent) < tolerance2) {

            resultVectorComponent = cross(u, zaxis);
        }

        double halfAngle = 0.5 * pi;
        double sinHalfAngle = sin(halfAngle);

        resultVectorComponent *= sinHalfAngle;

        k = 0.0;
        k_cos_theta = cos(halfAngle);

    } else {

        resultVectorComponent = cross(u, v);

    }

    quaternion(0) = k_cos_theta + k;
    quaternion(1) = resultVectorComponent(0);
    quaternion(2) = resultVectorComponent(1);
    quaternion(3) = resultVectorComponent(2);

    normalizeQuaternion(quaternion);
}


void ParallelCyclotronTracker::push(double h) {
    /* h   [ns]
     * dt1 [ns]
     * dt2 [ns]
     */
    IpplTimings::startTimer(IntegrationTimer_m);

    // h [ns] --> h [s]
    h *= 1.0e-9;

    std::list<CavityCrossData> cavCrossDatas;
    for(beamline_list::iterator sindex = ++(FieldDimensions.begin()); sindex != FieldDimensions.end(); ++sindex) {
        if(((*sindex)->first) == ElementBase::RFCAVITY) {
            RFCavity * cav = static_cast<RFCavity *>(((*sindex)->second).second);
            CavityCrossData ccd = {cav, cav->getSinAzimuth(), cav->getCosAzimuth(), cav->getPerpenDistance() * 0.001};
            cavCrossDatas.push_back(ccd);
        }
    }
    for(unsigned int i = 0; i < itsBunch_m->getLocalNum(); ++i) {
        Vector_t const oldR = itsBunch_m->R[i];
        double const gamma = sqrt(1.0 + dot(itsBunch_m->P[i], itsBunch_m->P[i]));
        double const c_gamma = Physics::c / gamma;
        Vector_t const v = itsBunch_m->P[i] * c_gamma;
        itsBunch_m->R[i] += h * v;
        for(std::list<CavityCrossData>::const_iterator ccd = cavCrossDatas.begin(); ccd != cavCrossDatas.end(); ++ccd) {
            double const distNew = (itsBunch_m->R[i][0] * ccd->sinAzimuth - itsBunch_m->R[i][1] * ccd->cosAzimuth) - ccd->perpenDistance;
            bool tagCrossing = false;
            double distOld;
            if(distNew <= 0.0) {
                distOld = (oldR[0] * ccd->sinAzimuth - oldR[1] * ccd->cosAzimuth) - ccd->perpenDistance;
                if(distOld > 0.0) tagCrossing = true;
            }
            if(tagCrossing) {
                double const dt1 = distOld / sqrt(dot(v, v));
                double const dt2 = h - dt1;

                // Retrack particle from the old postion to cavity gap point
                itsBunch_m->R[i] = oldR + dt1 * v;

                // Momentum kick
                //itsBunch_m->R[i] *= 1000.0; // RFkick uses [itsBunch_m->R] == mm
                RFkick(ccd->cavity, itsBunch_m->getT() * 1.0e9, dt1, i);
                //itsBunch_m->R[i] *= 0.001;

                itsBunch_m->R[i] += dt2 * itsBunch_m->P[i] * c_gamma;
            }
        }
    }
    itsBunch_m->setT(itsBunch_m->getT() + h);

    // Path length update
    double const dotP = dot(itsBunch_m->P[0], itsBunch_m->P[0]);
    double const gamma = sqrt(1.0 + dotP);
    PathLength_m += h * sqrt(dotP) * Physics::c / gamma;
    itsBunch_m->setLPath(PathLength_m);

    IpplTimings::stopTimer(IntegrationTimer_m);
}


void ParallelCyclotronTracker::kick(double h) {
    IpplTimings::startTimer(IntegrationTimer_m);

    BorisPusher pusher;
    double const q = itsBunch_m->Q[0] / q_e; // For now all particles have the same charge
    double const M = itsBunch_m->M[0] * 1.0e9; // For now all particles have the same rest energy

    for(unsigned int i = 0; i < itsBunch_m->getLocalNum(); ++i) {

        pusher.kick(itsBunch_m->R[i], itsBunch_m->P[i],
                    itsBunch_m->Ef[i], itsBunch_m->Bf[i],
                    h * 1.0e-9, M, q);
    }
    IpplTimings::stopTimer(IntegrationTimer_m);
}


void ParallelCyclotronTracker::borisExternalFields(double h) {
    // h in [ns]

    // push particles for first half step
    push(0.5 * h);

    // Evaluate external fields
    IpplTimings::startTimer(IntegrationTimer_m);
    for(unsigned int i = 0; i < itsBunch_m->getLocalNum(); ++i) {

        itsBunch_m->Ef[i] = Vector_t(0.0, 0.0, 0.0);
        itsBunch_m->Bf[i] = Vector_t(0.0, 0.0, 0.0);

        computeExternalFields_m(i, itsBunch_m->getT() * 1e9 /*ns*/,
                                itsBunch_m->Ef[i], itsBunch_m->Bf[i]);
    }
    IpplTimings::stopTimer(IntegrationTimer_m);

    // Kick particles for full step
    kick(h);

    // push particles for second half step
    push(0.5 * h);

    // apply the plugin elements: probe, collimator, stripper, septum
    applyPluginElements(h);
    // destroy particles if they are marked as Bin=-1 in the plugin elements or out of global apeture
    deleteParticle();
}


void ParallelCyclotronTracker::applyPluginElements(const double dt) {
    // Plugin Elements are all defined in mm, change beam to mm before applying
    itsBunch_m->R *= Vector_t(1000.0);

    for(beamline_list::iterator sindex = ++(FieldDimensions.begin()); sindex != FieldDimensions.end(); ++sindex) {
        ElementBase::ElementType type = ((*sindex)->first);
        if( type == ElementBase::CCOLLIMATOR ||
            type == ElementBase::PROBE       ||
            type == ElementBase::SEPTUM      ||
            type == ElementBase::STRIPPER) {

            PluginElement* element = static_cast<PluginElement *>(((*sindex)->second).second);
            bool flag = element->check(itsBunch_m,
                                       turnnumber_m,
                                       itsBunch_m->getT() * 1e9  /*[ns]*/,
                                       dt);
            if (type == ElementBase::STRIPPER && flag == true) {
                itsBunch_m->updateNumTotal();
                *gmsg << "* Total number of particles after stripping = " << itsBunch_m->getTotalNum() << endl;
            }
        }
    }

    itsBunch_m->R *= Vector_t(0.001);
}

bool ParallelCyclotronTracker::deleteParticle(){
    // Update immediately if any particles are lost during this step

    bool flagNeedUpdate = (min(itsBunch_m->Bin) < 0);
    allreduce(flagNeedUpdate, 1, std::logical_or<bool>());

    if(flagNeedUpdate) {
        size_t locLostParticleNum = 0;

        const int leb = itsBunch_m->getLastemittedBin();
        std::unique_ptr<size_t[]> localBinCount;

        if ( itsBunch_m->weHaveBins() ) {
            localBinCount = std::unique_ptr<size_t[]>(new size_t[leb]);
            for (int i = 0; i < leb; ++i)
                localBinCount[i] = 0;
        }

        for(unsigned int i = 0; i < itsBunch_m->getLocalNum(); i++) {
            if(itsBunch_m->Bin[i] < 0) {
                ++locLostParticleNum;
                itsBunch_m->destroy(1, i);
            } else if ( itsBunch_m->weHaveBins() ) {
                /* we need to count the local number of particles
                 * per energy bin
                 */
                ++localBinCount[itsBunch_m->Bin[i]];
            }
        }

        if ( itsBunch_m->weHaveBins() ) {
            // set the local bin count
            for (int i = 0; i < leb; ++i) {
                itsBunch_m->setLocalBinCount(localBinCount[i], i);
            }
        }

        /* We need to destroy the particles now
         * before we compute the means. We also
         * have to update the total number of particles
         * otherwise the statistics are wrong.
         */
        itsBunch_m->performDestroy(true);
        size_t totalnum = 0;
        size_t localnum = itsBunch_m->getLocalNum();
        allreduce(&localnum, &totalnum, 1, std::plus<size_t>());
        itsBunch_m->setTotalNum(totalnum);

        size_t globLostParticleNum = 0;
        reduce(locLostParticleNum, globLostParticleNum, 1, std::plus<size_t>());

        *gmsg << "At step " << step_m
              << ", lost "  << globLostParticleNum << " particles "
              << "on stripper, collimator, septum, or out of cyclotron aperture"
              << endl;

        if (totalnum == 0) {
            return flagNeedUpdate;
        }

        Vector_t const meanR = calcMeanR();
        Vector_t const meanP = calcMeanP();

        // Bunch (local) azimuth at meanR w.r.t. y-axis
        double const phi = calculateAngle(meanP(0), meanP(1)) - 0.5 * pi;

        // Bunch (local) elevation at meanR w.r.t. xy plane
        double const psi = 0.5 * pi - acos(meanP(2) / sqrt(dot(meanP, meanP)));

        // For statistics data, the bunch is transformed into a local coordinate system
        // with meanP in direction of y-axis -DW
        globalToLocal(itsBunch_m->R, phi, psi, meanR);
        globalToLocal(itsBunch_m->P, phi, psi, Vector_t(0.0)); // P should be rotated, but not shifted

        //itsBunch_m->R *= Vector_t(0.001); // mm --> m

        // Now destroy particles and update pertinent parameters in local frame
        // Note that update() will be called within boundp() -DW
        itsBunch_m->boundp();
        //itsBunch_m->update();

        itsBunch_m->calcBeamParameters();

        //itsBunch_m->R *= Vector_t(1000.0); // m --> mm

        localToGlobal(itsBunch_m->R, phi, psi, meanR);
        localToGlobal(itsBunch_m->P, phi, psi, Vector_t(0.0));

        // If particles were deleted, recalculate bingamma and reset BinID for remaining particles
        if ( itsBunch_m->weHaveBins() )
            updateParticleBins_m();
    }
    return flagNeedUpdate;
}

void ParallelCyclotronTracker::initTrackOrbitFile() {

    std::string f = OpalData::getInstance()->getInputBasename() + std::string("-trackOrbit.dat");

    outfTrackOrbit_m.setf(std::ios::scientific, std::ios::floatfield);
    outfTrackOrbit_m.precision(8);

    if(myNode_m == 0) {

        if(OpalData::getInstance()->inRestartRun()) {

            outfTrackOrbit_m.open(f.c_str(), std::ios::app);
            outfTrackOrbit_m << "# Restart at integration step " << itsBunch_m->getLocalTrackStep() << std::endl;
        } else {

            outfTrackOrbit_m.open(f.c_str());
            outfTrackOrbit_m << "# The six-dimensional phase space data in the global Cartesian coordinates" << std::endl;
            outfTrackOrbit_m << "# Part. ID    x [m]       beta_x*gamma       y [m]      beta_y*gamma        z [m]      beta_z*gamma" << std::endl;
        }
    }
}

void ParallelCyclotronTracker::initDistInGlobalFrame() {
    if(!OpalData::getInstance()->inRestartRun()) {
        // Start a new run (no restart)

        double const initialReferenceTheta = referenceTheta * Physics::deg2rad;
        PathLength_m = 0.0;

        // TODO: Replace with TracerParticle
        // Force the initial phase space values of the particle with ID = 0 to zero,
        // to set it as a reference particle.
        if(initialTotalNum_m > 2) {
            for(size_t i = 0; i < initialLocalNum_m; ++i) {
                if(itsBunch_m->ID[i] == 0) {
                    itsBunch_m->R[i] = Vector_t(0.0);
                    itsBunch_m->P[i] = Vector_t(0.0);
                }
            }
        }

        // Initialize global R
        //itsBunch_m->R *= Vector_t(1000.0); // m --> mm

        // NEW OPAL 2.0: Immediately change to m -DW
        Vector_t const initMeanR = Vector_t(0.001 * referenceR * cosRefTheta_m, // mm --> m
                                            0.001 * referenceR * sinRefTheta_m, // mm --> m
                                            0.001 * referenceZ);                // mm --> m

        localToGlobal(itsBunch_m->R, initialReferenceTheta, initMeanR);

        // Initialize global P (Cartesian, but input P_ref is in Pr, Ptheta, Pz,
        // so translation has to be done before the rotation this once)
        // Cave: In the local frame, the positive y-axis is the direction of movement -DW
        for(size_t i = 0; i < initialLocalNum_m; ++i) {
            itsBunch_m->P[i](0) += referencePr;
            itsBunch_m->P[i](1) += referencePt;
            itsBunch_m->P[i](2) += referencePz;
        }

        // Out of the three coordinates of meanR (R, Theta, Z) only the angle
        // changes the momentum vector...
        localToGlobal(itsBunch_m->P, initialReferenceTheta);

        // Initialize the bin number of the first bunch to 0
        for(size_t i = 0; i < initialLocalNum_m; ++i) {
            itsBunch_m->Bin[i] = 0;
        }

        // Backup initial distribution if multi bunch mode
        if ((initialTotalNum_m > 2) && (numBunch_m > 1) && (multiBunchMode_m == MB_MODE::FORCE))
            saveOneBunch();

        // Else: Restart from the distribution in the h5 file
    } else {

        // Do a local frame restart (we have already checked that the old h5 file was saved in local
        // frame as well).
        if((Options::psDumpFrame != Options::GLOBAL)) {

            *gmsg << "* Restart in the local frame" << endl;

            PathLength_m = 0.0;
            //itsBunch_m->R *= Vector_t(1000.0); // m --> mm

            // referenceR and referenceZ are already in mm
            // New OPAL 2.0: Init in m -DW
            Vector_t const initMeanR = Vector_t(0.001 * referenceR * cosRefTheta_m,
                                                0.001 * referenceR * sinRefTheta_m,
                                                0.001 * referenceZ);

            // Do the tranformations
            localToGlobal(itsBunch_m->R, referencePhi, referencePsi, initMeanR);
            localToGlobal(itsBunch_m->P, referencePhi, referencePsi);

            // Initialize the bin number of the first bunch to 0
            for(size_t i = 0; i < initialLocalNum_m; ++i) {
                itsBunch_m->Bin[i] = 0;
            }

            // Or do a global frame restart (no transformations necessary)
        } else {

            *gmsg << "* Restart in the global frame" << endl;

            PathLength_m = itsBunch_m->getLPath();
            //itsBunch_m->R *= Vector_t(1000.0); // m --> mm
        }
    }

    // ------------ Get some Values ---------------------------------------------------------- //
    Vector_t const meanR = calcMeanR();
    Vector_t const meanP = calcMeanP();

    // Bunch (local) azimuth at meanR w.r.t. y-axis
    double const phi = calculateAngle(meanP(0), meanP(1)) - 0.5 * pi;

    // Bunch (local) elevation at meanR w.r.t. xy plane
    double const psi = 0.5 * pi - acos(meanP(2) / sqrt(dot(meanP, meanP)));

    if(multiBunchMode_m == MB_MODE::AUTO) {

        RLastTurn_m = std::sqrt(meanR[0] * meanR[0] + meanR[1] * meanR[1]);  // [m]
        RThisTurn_m = RLastTurn_m;

        if(OpalData::getInstance()->inRestartRun()) {

            *gmsg << "Radial position at restart position = ";

        } else {

            *gmsg << "Initial radial position = ";
        }
        // New OPAL 2.0: Init in m -DW
        *gmsg << RThisTurn_m << " m" << endl;
    }

    // Do boundp and repartition in the local frame at beginning of this run
    globalToLocal(itsBunch_m->R, phi, psi, meanR);
    globalToLocal(itsBunch_m->P, phi, psi); // P should be rotated, but not shifted

    //itsBunch_m->R *= Vector_t(0.001); // mm --> m

    itsBunch_m->boundp();

    checkNumPart(std::string("* Before repartition: "));
    repartition();
    checkNumPart(std::string("* After repartition:  "));

    itsBunch_m->calcBeamParameters();

    *gmsg << endl << "* *********************** Bunch information in local frame: ************************";
    *gmsg << *itsBunch_m << endl;

    //itsBunch_m->R *= Vector_t(1000.0); // m --> mm

    localToGlobal(itsBunch_m->R, phi, psi, meanR);
    localToGlobal(itsBunch_m->P, phi, psi);

    // Save initial distribution if not a restart
    if(!OpalData::getInstance()->inRestartRun()) {

        step_m -= 1;

        bunchDumpPhaseSpaceData();
        bunchDumpStatData();

        step_m += 1;
    }

    // Print out the Bunch information at beginning of the run. Because the bunch information
    // displays in units of m we have to change back and forth one more time.
    //itsBunch_m->R *= Vector_t(0.001); // mm --> m

    itsBunch_m->calcBeamParameters();

    *gmsg << endl << "* *********************** Bunch information in global frame: ***********************";
    *gmsg << *itsBunch_m << endl;

    //itsBunch_m->R *= Vector_t(1000.0); // m --> mm
}

//TODO: This can be completely rewritten with TracerParticle -DW
void ParallelCyclotronTracker::singleParticleDump() {
    IpplTimings::startTimer(DumpTimer_m);

    if(Ippl::getNodes() > 1 ) {

        double x;
        int id;
        dvector_t tmpr;
        ivector_t tmpi;

        int tag = Ippl::Comm->next_tag(IPPL_APP_TAG4, IPPL_APP_CYCLE);

        // for all nodes, find the location of particle with ID = 0 & 1 in bunch containers
        int found[2] = {-1, -1};
        int counter = 0;

        for(size_t i = 0; i < itsBunch_m->getLocalNum(); ++i) {
            if(itsBunch_m->ID[i] == 0) {
                found[counter] = i;
                counter++;
            }
            if(itsBunch_m->ID[i] == 1) {
                found[counter] = i;
                counter++;
            }
        }

        if(myNode_m == 0) {
            int notReceived = Ippl::getNodes() - 1;
            int numberOfPart = 0;
            // receive other nodes
            while(notReceived > 0) {

                int node = COMM_ANY_NODE;
                Message *rmsg =  Ippl::Comm->receive_block(node, tag);

                if(rmsg == nullptr)
                    ERRORMSG("Could not receive from client nodes in main." << endl);

                --notReceived;

                rmsg->get(&numberOfPart);

                for(int i = 0; i < numberOfPart; ++i) {
                    rmsg->get(&id);
                    tmpi.push_back(id);
                    for(int ii = 0; ii < 6; ii++) {
                        rmsg->get(&x);
                        tmpr.push_back(x);
                    }
                }
                delete rmsg;
            }
            // own node
            for(int i = 0; i < counter; ++i) {

                tmpi.push_back(itsBunch_m->ID[found[i]]);

                for(int j = 0; j < 3; ++j) {

                    tmpr.push_back(itsBunch_m->R[found[i]](j));
                    tmpr.push_back(itsBunch_m->P[found[i]](j));
                }
            }
            // store
            dvector_t::iterator itParameter = tmpr.begin();

            for(auto tmpid : tmpi) {

                outfTrackOrbit_m << "ID" << tmpid;

                if (tmpid == 0) { // for stat file
                    itsBunch_m->RefPartR_m[0] = *itParameter;
                    itsBunch_m->RefPartR_m[1] = *(itParameter + 2);
                    itsBunch_m->RefPartR_m[2] = *(itParameter + 4);
                    itsBunch_m->RefPartP_m[0] = *(itParameter + 1);
                    itsBunch_m->RefPartP_m[1] = *(itParameter + 3);
                    itsBunch_m->RefPartP_m[2] = *(itParameter + 5);
                }
                for(int ii = 0; ii < 6; ii++) {

                    outfTrackOrbit_m << " " << *itParameter;
                    ++itParameter;
                }

                outfTrackOrbit_m << std::endl;
            }
        } else {
            // for other nodes
            Message *smsg = new Message();
            smsg->put(counter);

            for(int i = 0; i < counter; i++) {

                smsg->put(itsBunch_m->ID[found[i]]);

                for(int j = 0; j < 3; j++) {

                    smsg->put(itsBunch_m->R[found[i]](j));
                    smsg->put(itsBunch_m->P[found[i]](j));
                }
            }

            if(!Ippl::Comm->send(smsg, 0, tag))
                ERRORMSG("Ippl::Comm->send(smsg, 0, tag) failed " << endl);
        }

        Ippl::Comm->barrier();

    } else {

        for(size_t i = 0; i < itsBunch_m->getLocalNum(); i++) {

            if(itsBunch_m->ID[i] == 0 || itsBunch_m->ID[i] == 1) {

                outfTrackOrbit_m << "ID" << itsBunch_m->ID[i] << " ";
                outfTrackOrbit_m << itsBunch_m->R[i](0) << " " << itsBunch_m->P[i](0) << " ";
                outfTrackOrbit_m << itsBunch_m->R[i](1) << " " << itsBunch_m->P[i](1) << " ";
                outfTrackOrbit_m << itsBunch_m->R[i](2) << " " << itsBunch_m->P[i](2) << std::endl;

                if (itsBunch_m->ID[i] == 0) { // for stat file
                    itsBunch_m->RefPartR_m = itsBunch_m->R[i];
                    itsBunch_m->RefPartP_m = itsBunch_m->P[i];
                }
            }
        }
    }

    IpplTimings::stopTimer(DumpTimer_m);
}

void ParallelCyclotronTracker::bunchDumpStatData(){

    // dump stat file per bin in case of multi-bunch mode
    if ( multiBunchMode_m != MB_MODE::NONE ) {
        double phi = 0.0, psi = 0.0;
        Vector_t meanR = calcMeanR();

        // Bunch (global) angle w.r.t. x-axis (cylinder coordinates)
        double theta = calculateAngle(meanR(0), meanR(1)) * Physics::rad2deg;

        if(Options::psDumpFrame != Options::GLOBAL) {
            Vector_t meanP = calcMeanP();

            // Bunch (local) azimuth at meanR w.r.t. y-axis
            phi = calculateAngle(meanP(0), meanP(1)) - 0.5 * pi;

            // Bunch (local) elevation at meanR w.r.t. xy plane
            psi = 0.5 * pi - acos(meanP(2) / sqrt(dot(meanP, meanP)));

            // Rotate so Pmean is in positive y direction. No shift, so that normalized emittance and
            // unnormalized emittance as well as centroids are calculated correctly in
            // PartBunch::calcBeamParameters()
            globalToLocal(itsBunch_m->R, phi, psi, meanR);
            globalToLocal(itsBunch_m->P, phi, psi);
        }

        bunchDumpStatDataPerBin(theta);

        if(Options::psDumpFrame != Options::GLOBAL) {
            localToGlobal(itsBunch_m->R, phi, psi, meanR);
            localToGlobal(itsBunch_m->P, phi, psi);
        }

        return;
    }

    IpplTimings::startTimer(DumpTimer_m);

    /*
      in case before a bunchDumpPhasespaceData has
      happen, the calcBeamParameters() is not needed.
      maybe we have to introduce a dirty-stat-data flag.
      this flag would be set if the bunch is dirty and
      reset by calcBeamParameters
    */
    itsBunch_m->R *= Vector_t(0.001); // mm --> m
    itsBunch_m->calcBeamParameters();
    itsBunch_m->R *= Vector_t(1000.0); // m --> mm

    // --------------------------------- Get some Values ---------------------------------------- //
    double const E = itsBunch_m->get_meanKineticEnergy();
    double const temp_t = itsBunch_m->getT() * 1e9; // s -> ns
    Vector_t meanR;
    Vector_t meanP;
    if (Options::psDumpFrame == Options::BUNCH_MEAN) {
        meanR = calcMeanR();
        meanP = calcMeanP();
    } else if (itsBunch_m->getLocalNum() > 0) {
        meanR = itsBunch_m->R[0];
        meanP = itsBunch_m->P[0];
    }
    double phi = 0;
    double psi = 0;

    // Bunch (global) angle w.r.t. x-axis (cylinder coordinates)
    double azimuth = calculateAngle(meanR(0), meanR(1)) * Physics::rad2deg;

    // --------------  Calculate the external fields at the center of the bunch ----------------- //
    beamline_list::iterator DumpSindex = FieldDimensions.begin();

    extE_m = Vector_t(0.0, 0.0, 0.0);
    extB_m = Vector_t(0.0, 0.0, 0.0);

    (((*DumpSindex)->second).second)->apply(meanR, meanP, temp_t, extE_m, extB_m);

    // If we are saving in local frame, bunch and fields at the bunch center have to be rotated
    // TODO: Make decision if we maybe want to always save statistics data in local frame? -DW
    if(Options::psDumpFrame != Options::GLOBAL) {
        // -------------------- ----------- Do Transformations ---------------------------------- //
        // Bunch (local) azimuth at meanR w.r.t. y-axis
        phi = calculateAngle(meanP(0), meanP(1)) - 0.5 * pi;

        // Bunch (local) elevation at meanR w.r.t. xy plane
        psi = 0.5 * pi - acos(meanP(2) / sqrt(dot(meanP, meanP)));

        // Rotate so Pmean is in positive y direction. No shift, so that normalized emittance and
        // unnormalized emittance as well as centroids are calculated correctly in
        // PartBunch::calcBeamParameters()
        globalToLocal(extB_m, phi, psi);
        globalToLocal(extE_m, phi, psi);
        globalToLocal(itsBunch_m->R, phi, psi);
        globalToLocal(itsBunch_m->P, phi, psi);
    }

    //itsBunch_m->R *= Vector_t(0.001); // mm -> m

    FDext_m[0] = extB_m / 10.0; // kgauss --> T
    FDext_m[1] = extE_m;        // kV/mm? -DW

    // Save the stat file
    itsDataSink->writeStatData(itsBunch_m, FDext_m, E, azimuth);

    //itsBunch_m->R *= Vector_t(1000.0); // m -> mm

    // If we are in local mode, transform back after saving
    if(Options::psDumpFrame != Options::GLOBAL) {
        localToGlobal(itsBunch_m->R, phi, psi);
        localToGlobal(itsBunch_m->P, phi, psi);
    }

    IpplTimings::stopTimer(DumpTimer_m);
}


void ParallelCyclotronTracker::bunchDumpStatDataPerBin(const double& azimuth) {
    IpplTimings::startTimer(DumpTimer_m);

    for (short b = 0; b < BunchCount_m; ++b) {

        MultiBunchDump::beaminfo_t binfo;

        if ( itsBunch_m->calcBunchBeamParameters(binfo, b) ) {
            binfo.azimuth = azimuth;
            itsMBDump_m->writeData(binfo, b);
        }
    }

    IpplTimings::stopTimer(DumpTimer_m);
}


void ParallelCyclotronTracker::bunchDumpPhaseSpaceData() {
    // --------------------------------- Particle dumping --------------------------------------- //
    // Note: Don't dump when
    // 1. after one turn
    // in order to sychronize the dump step for multi-bunch and single bunch for comparison
    // with each other during post-process phase.
    // ------------------------------------------------------------------------------------------ //
    IpplTimings::startTimer(DumpTimer_m);

    // --------------------------------- Get some Values ---------------------------------------- //
    double const temp_t = itsBunch_m->getT() * 1.0e9; // s -> ns

    Vector_t meanR;
    Vector_t meanP;

    // in case of multi-bunch mode take always bunch mean (although it takes all bunches)
    if (Options::psDumpFrame == Options::BUNCH_MEAN || multiBunchMode_m != MB_MODE::NONE ) {
        meanR = calcMeanR();
        meanP = calcMeanP();
    } else if (itsBunch_m->getLocalNum() > 0) {
        meanR = itsBunch_m->R[0];
        meanP = itsBunch_m->P[0];
    }

    double const betagamma_temp = sqrt(dot(meanP, meanP));
    double const E = itsBunch_m->get_meanKineticEnergy();

    // Bunch (global) angle w.r.t. x-axis (cylinder coordinates)
    double const theta = atan2(meanR(1), meanR(0));

    // Bunch (local) azimuth at meanR w.r.t. y-axis
    double const phi = calculateAngle(meanP(0), meanP(1)) - 0.5 * pi;

    // Bunch (local) elevation at meanR w.r.t. xy plane
    double const psi = 0.5 * pi - acos(meanP(2) / sqrt(dot(meanP, meanP)));

    // ---------------- Re-calculate reference values in format of input values ----------------- //
    // Position:
    // New OPAL 2.0: Init in m (back to mm just for output) -DW
    referenceR = 1000.0 * sqrt(meanR(0) * meanR(0) + meanR(1) * meanR(1));
    referenceTheta = theta / Physics::deg2rad;
    referenceZ = 1000.0 * meanR(2);

    // Momentum in Theta-hat, R-hat, Z-hat coordinates at position meanR:
    referencePtot = betagamma_temp;
    referencePz = meanP(2);
    referencePr = meanP(0) * cos(theta) + meanP(1) * sin(theta);
    referencePt = sqrt(referencePtot * referencePtot - \
                       referencePz * referencePz - referencePr * referencePr);

    // -----  Calculate the external fields at the center of the bunch (Cave: Global Frame) ----- //
    beamline_list::iterator DumpSindex = FieldDimensions.begin();

    extE_m = Vector_t(0.0, 0.0, 0.0);
    extB_m = Vector_t(0.0, 0.0, 0.0);

    (((*DumpSindex)->second).second)->apply(meanR, meanP, temp_t, extE_m, extB_m);
    FDext_m[0] = extB_m * 0.1; // kgauss --> T
    FDext_m[1] = extE_m;

    // -------------- If flag psDumpFrame is global, dump bunch in global frame ------------- //
    if (Options::psDumpFreq < std::numeric_limits<int>::max() ){
        bool dumpLocalFrame    = true;
        std::string dumpString = "local";
        if (Options::psDumpFrame == Options::GLOBAL) {
            dumpLocalFrame = false;
            dumpString     = "global";
        }

        if (dumpLocalFrame == true) {
            // ---------------- If flag psDumpFrame is local, dump bunch in local frame ---------------- //

            // The bunch is transformed into a local coordinate system with meanP in direction of y-axis
            globalToLocal(itsBunch_m->R, phi, psi, meanR);
            //globalToLocal(itsBunch_m->R, phi, psi, meanR * Vector_t(0.001));
            globalToLocal(itsBunch_m->P, phi, psi); // P should only be rotated

            globalToLocal(extB_m, phi, psi);
            globalToLocal(extE_m, phi, psi);
        }

        FDext_m[0] = extB_m * 0.1; // kgauss --> T
        FDext_m[1] = extE_m;

        lastDumpedStep_m = itsDataSink->writePhaseSpace_cycl(itsBunch_m, // Local and in m
                                                             FDext_m, E,
                                                             referencePr,
                                                             referencePt,
                                                             referencePz,
                                                             referenceR,
                                                             referenceTheta,
                                                             referenceZ,
                                                             phi / Physics::deg2rad, // P_mean azimuth
                                                             // at ref. R/Th/Z
                                                             psi / Physics::deg2rad, // P_mean elevation
                                                             // at ref. R/Th/Z
                                                             dumpLocalFrame);        // Flag localFrame

        if (dumpLocalFrame == true) {
            // Return to global frame
            localToGlobal(itsBunch_m->R, phi, psi, meanR);
            //localToGlobal(itsBunch_m->R, phi, psi, meanR * Vector_t(0.001));
            localToGlobal(itsBunch_m->P, phi, psi);
        }


        // Tell user in which mode we are dumping
        // New: no longer dumping for num part < 3, omit phase space dump number info
        if (lastDumpedStep_m == -1){
          *gmsg << endl << "* Integration step " << step_m + 1
                << " (no phase space dump for <= 2 particles)" << endl;
        } else {
          *gmsg << endl << "* Phase space dump " << lastDumpedStep_m
                << " (" << dumpString << " frame) at integration step " << step_m + 1 << endl;
        }
    }

    // Print dump information on screen
    *gmsg << "* T = " << temp_t << " ns"
          << ", Live Particles: " << itsBunch_m->getTotalNum() << endl;
    *gmsg << "* E = " << E << " MeV"
          << ", beta * gamma = " << betagamma_temp << endl;
    *gmsg << "* Bunch position: R =  " << referenceR << " mm"
          << ", Theta = " << referenceTheta << " Deg"
          << ", Z = " << referenceZ << " mm" << endl;
    *gmsg << "* Local Azimuth = " << phi / Physics::deg2rad << " Deg"
          << ", Local Elevation = " << psi / Physics::deg2rad << " Deg" << endl;

    IpplTimings::stopTimer(DumpTimer_m);
}


void ParallelCyclotronTracker::update_m(double& t, const double& dt,
                                        const bool& dumpEachTurn)
{
    // Reference parameters
    t += dt;
    itsBunch_m->setT(t * 1.0e-9);
    itsBunch_m->setLocalTrackStep((step_m + 1));
    if (!(step_m + 1 % 1000))
        *gmsg << "Step " << step_m + 1 << endl;

    if (itsBunch_m->getLocalNum() > 0) {
        double tempP2 = dot(itsBunch_m->P[0], itsBunch_m->P[0]);
        double tempGamma = sqrt(1.0 + tempP2);
        double tempBeta = sqrt(tempP2) / tempGamma;

        PathLength_m += c_mmtns * dt / 1000.0 * tempBeta; // unit: m
        itsBunch_m->setLPath(PathLength_m);
    }

    // Here is global frame, don't do itsBunch_m->boundp();

    if (itsBunch_m->getTotalNum()>0) {
        // Only dump last step if we have particles left.
        // Check separately for phase space (ps) and statistics (stat) data dump frequency
        if ( mode_m != MODE::SEO && ( ((step_m + 1) % Options::psDumpFreq == 0) ||
                                      (Options::psDumpEachTurn && dumpEachTurn)))
        {
            // Write phase space data to h5 file
            bunchDumpPhaseSpaceData();
        }

        if ( mode_m != MODE::SEO && ( ((step_m + 1) % Options::statDumpFreq == 0) ||
                                      (Options::psDumpEachTurn && dumpEachTurn)))
        {
            // Write statistics data to stat file
            bunchDumpStatData();
        }
    }
}


std::tuple<double, double, double> ParallelCyclotronTracker::initializeTracking_m() {
    // Read in some control parameters
    setup_m.scSolveFreq         = (spiral_flag) ? 1 : Options::scSolveFreq;
    setup_m.stepsPerTurn        = itsBunch_m->getStepsPerTurn();

    // Define 3 special azimuthal angles where dump particle's six parameters
    // at each turn into 3 ASCII files. only important in single particle tracking
    azimuth_angle_m.resize(3);
    azimuth_angle_m[0] = 0.0;
    azimuth_angle_m[1] = 22.5 * Physics::deg2rad;
    azimuth_angle_m[2] = 45.0 * Physics::deg2rad;

    double harm       = getHarmonicNumber();
    double dt         = itsBunch_m->getdT() * 1.0e9 * harm; // time step size (s --> ns)
    double t          = itsBunch_m->getT() * 1.0e9;               // current time   (s --> ns)


    double oldReferenceTheta      = referenceTheta * Physics::deg2rad; // init here, reset each step
    setup_m.deltaTheta            = pi / (setup_m.stepsPerTurn);    // half of the average angle per step

    //int stepToLastInj = itsBunch_m->getSteptoLastInj(); // TODO: Do we need this? -DW

    // Record how many bunches have already been injected. ONLY FOR MBM
    BunchCount_m = itsBunch_m->getNumBunch();

    initTrackOrbitFile();

    // Get data from h5 file for restart run and reset current step
    // to last step of previous simulation
    if(OpalData::getInstance()->inRestartRun()) {

        restartStep0_m = itsBunch_m->getLocalTrackStep();
        step_m = restartStep0_m;

        *gmsg << "* Restart at integration step " << restartStep0_m << endl;
    }

    setup_m.stepsNextCheck = step_m + setup_m.stepsPerTurn; // Steps to next check for transition

    initDistInGlobalFrame();

    if ( multiBunchMode_m != MB_MODE::NONE ) {
        updateParticleBins_m();
    }

    turnnumber_m = 1;

    // --- Output to user --- //
    *gmsg << "* Beginning of this run is at t = " << t << " [ns]" << endl;
    *gmsg << "* The time step is set to dt = " << dt << " [ns]" << endl;

    if(numBunch_m > 1) {

        *gmsg << "* MBM: Time interval between neighbour bunches is set to "
              << setup_m.stepsPerTurn * dt << "[ns]" << endl;
        *gmsg << "* MBM: The particles energy bin reset frequency is set to "
              << Options::rebinFreq << endl;
    }

    *gmsg << "* Single particle trajectory dump frequency is set to " << Options::sptDumpFreq << endl;
    *gmsg << "* The frequency to solve space charge fields is set to " << setup_m.scSolveFreq << endl;
    *gmsg << "* The repartition frequency is set to " << Options::repartFreq << endl;


    switch ( mode_m )
    {
        case MODE::SEO:
            *gmsg << endl;
            *gmsg << "* ------------------------- STATIC EQUILIBRIUM ORBIT MODE ----------------------------- *" << endl
                  << "* Instruction: When the total particle number is equal to 2, SEO mode is triggered      *" << endl
                  << "* automatically. This mode does NOT include any RF cavities. The initial distribution   *" << endl
                  << "* file must be specified. In the file the first line is for reference particle and the  *" << endl
                  << "* second line is for off-center particle. The tune is calculated by FFT routines based  *" << endl
                  << "* on these two particles.                                                               *" << endl
                  << "* ---------------- NOTE: SEO MODE ONLY WORKS SERIALLY ON SINGLE NODE ------------------ *" << endl;

            if(Ippl::getNodes() != 1)
                throw OpalException("Error in ParallelCyclotronTracker::initializeTracking_m",
                                    "SEO MODE ONLY WORKS SERIALLY ON SINGLE NODE!");
            break;
        case MODE::SINGLE:
            *gmsg << endl;
            *gmsg << "* ------------------------------ SINGLE PARTICLE MODE --------------------------------- *" << endl
                  << "* Instruction: When the total particle number is equal to 1, single particle mode is    *" << endl
                  << "* triggered automatically. The initial distribution file must be specified which should *" << endl
                  << "* contain only one line for the single particle                                         *" << endl
                  << "* ---------NOTE: SINGLE PARTICLE MODE ONLY WORKS SERIALLY ON A SINGLE NODE ------------ *" << endl;

            if(Ippl::getNodes() != 1)
                throw OpalException("Error in ParallelCyclotronTracker::initializeTracking_m",
                                    "SINGLE PARTICLE MODE ONLY WORKS SERIALLY ON A SINGLE NODE!");

            // For single particle mode open output files
            openFiles(azimuth_angle_m.size() + 1, OpalData::getInstance()->getInputBasename());
            break;
        case MODE::BUNCH:
            break;
        case MODE::UNDEFINED:
        default:
            throw OpalException("ParallelCyclotronTracker::initializeTracking_m()",
                                "No such tracking mode.");
    }

    return std::make_tuple(t, dt, oldReferenceTheta);
}


void ParallelCyclotronTracker::finalizeTracking_m(dvector_t& Ttime,
                                                  dvector_t& Tdeltr,
                                                  dvector_t& Tdeltz, ivector_t& TturnNumber)
{
    for(size_t ii = 0; ii < (itsBunch_m->getLocalNum()); ii++) {
        if(itsBunch_m->ID[ii] == 0) {
            double FinalMomentum2 = pow(itsBunch_m->P[ii](0), 2.0) + pow(itsBunch_m->P[ii](1), 2.0) + pow(itsBunch_m->P[ii](2), 2.0);
            double FinalEnergy = (sqrt(1.0 + FinalMomentum2) - 1.0) * itsBunch_m->getM() * 1.0e-6;
            *gmsg << "* Final energy of reference particle = " << FinalEnergy << " [MeV]" << endl;
            *gmsg << "* Total phase space dump number(includes the initial distribution) = " << lastDumpedStep_m + 1 << endl;
            *gmsg << "* One can restart simulation from the last dump step (--restart " << lastDumpedStep_m << ")" << endl;
        }
    }

    Ippl::Comm->barrier();

    switch ( mode_m )
    {
        case MODE::SEO:
        {
            // Calculate tunes after tracking.
            *gmsg << endl;
            *gmsg << "* **************** The result for tune calulation (NO space charge) ******************* *" << endl
                  << "* Number of tracked turns: " << TturnNumber.back() << endl;
            double nur, nuz;
            getTunes(Ttime, Tdeltr, Tdeltz, TturnNumber.back(), nur, nuz);
            break;
        }
        case MODE::SINGLE:
            closeFiles();
            // fall through
        case MODE::BUNCH:  // we do nothing
        case MODE::UNDEFINED:
        default:
        {
            // not for multibunch
            if (!(itsBunch_m->weHaveBins())) {
                *gmsg << "*" << endl;
                *gmsg << "* Finished during turn " << turnnumber_m << " (" << turnnumber_m - 1 << " turns completed)" << endl;
                *gmsg << "* Cave: Turn number is not correct for restart mode"<< endl;
            }
        }
    }

    Ippl::Comm->barrier();

    if (myNode_m == 0)
        outfTrackOrbit_m.close();

    *gmsg << endl << "* *********************** Bunch information in global frame: ***********************";

    if (itsBunch_m->getTotalNum() > 0){
        // Print out the Bunch information at end of the run.
        itsBunch_m->calcBeamParameters();
        *gmsg << *itsBunch_m << endl;

    } else {

        *gmsg << endl << "* No Particles left in bunch!" << endl;
        *gmsg << "* **********************************************************************************" << endl;

    }
}


void ParallelCyclotronTracker::seoMode_m(double& t, const double dt, bool& dumpEachTurn,
                                         dvector_t& Ttime, dvector_t& Tdeltr,
                                         dvector_t& Tdeltz, ivector_t& TturnNumber)
{

    // 2 particles: Trigger SEO mode
    // (Switch off cavity and calculate betatron oscillation tuning)
    double r_tuning[2], z_tuning[2] ;

    IpplTimings::startTimer(IntegrationTimer_m);
    for(size_t i = 0; i < (itsBunch_m->getLocalNum()); i++) {

        if((step_m % Options::sptDumpFreq == 0)) {

            outfTrackOrbit_m << "ID" << (itsBunch_m->ID[i]);
            outfTrackOrbit_m << " " << itsBunch_m->R[i](0)
                             << " " << itsBunch_m->P[i](0)
                             << " " << itsBunch_m->R[i](1)
                             << " " << itsBunch_m->P[i](1)
                             << " " << itsBunch_m->R[i](2)
                             << " " << itsBunch_m->P[i](2)
                             << std::endl;
        }

        double OldTheta = calculateAngle(itsBunch_m->R[i](0), itsBunch_m->R[i](1));
        r_tuning[i] = itsBunch_m->R[i](0) * cos(OldTheta) +
                      itsBunch_m->R[i](1) * sin(OldTheta);

        z_tuning[i] = itsBunch_m->R[i](2);


        // Integrate for one step in the lab Cartesian frame (absolute value).
        itsStepper_mp->advance(itsBunch_m, i, t, dt);

        if( (i == 0) && (step_m > 10) && ((step_m%setup_m.stepsPerTurn) == 0))
            turnnumber_m++;

    } // end for: finish one step tracking for all particles for initialTotalNum_m == 2 mode
    IpplTimings::stopTimer(IntegrationTimer_m);

    // store dx and dz for future tune calculation if higher precision needed, reduce freqSample.
    if(step_m % Options::sptDumpFreq == 0) {

        Ttime.push_back(t * 1.0e-9);
        Tdeltz.push_back(z_tuning[1]);
        Tdeltr.push_back(r_tuning[1] - r_tuning[0]);
        TturnNumber.push_back(turnnumber_m);
    }
}


void ParallelCyclotronTracker::singleMode_m(double& t, const double dt,
                                            bool& dumpEachTurn, double& oldReferenceTheta) {
    // 1 particle: Trigger single particle mode

    // ********************************************************************************** //
    // * This was moved here b/c collision should be tested before the actual           * //
    // * timestep (bgf_main_collision_test() predicts the next step automatically)      * //

    // apply the plugin elements: probe, collimator, stripper, septum
    applyPluginElements(dt);

    // check if we lose particles at the boundary
    bgf_main_collision_test();
    // ********************************************************************************** //

    if (itsBunch_m->getLocalNum() == 0) return; // might happen if particle is in collimator

    IpplTimings::startTimer(IntegrationTimer_m);

    unsigned int i = 0; // we only have a single particle

    if ( step_m % Options::sptDumpFreq == 0 ) {

        outfTrackOrbit_m << "ID" <<itsBunch_m->ID[i]
                         << " " << itsBunch_m->R[i](0)
                         << " " << itsBunch_m->P[i](0)
                         << " " << itsBunch_m->R[i](1)
                         << " " << itsBunch_m->P[i](1)
                         << " " << itsBunch_m->R[i](2)
                         << " " << itsBunch_m->P[i](2)
                         << std::endl;
    }

    double temp_meanTheta = calculateAngle2(itsBunch_m->R[i](0),
                                            itsBunch_m->R[i](1)); // [-pi, pi]


    dumpThetaEachTurn_m(t, itsBunch_m->R[i], itsBunch_m->P[i],
                        temp_meanTheta, dumpEachTurn);

    dumpAzimuthAngles_m(t, itsBunch_m->R[i], itsBunch_m->P[i],
                        oldReferenceTheta, temp_meanTheta);

    oldReferenceTheta = temp_meanTheta;

    // used for gap crossing checking
    Vector_t Rold = itsBunch_m->R[i]; // [x,y,z]    (mm)
    Vector_t Pold = itsBunch_m->P[i]; // [px,py,pz] (beta*gamma)

    // integrate for one step in the lab Cartesian frame (absolute value).
    itsStepper_mp->advance(itsBunch_m, i, t, dt);

    // If gap crossing happens, do momenta kicking (not if gap crossing just happened)
    if (itsBunch_m->cavityGapCrossed[i] == true)
        itsBunch_m->cavityGapCrossed[i] = false;
    else
        gapCrossKick_m(i, t, dt, Rold, Pold);
    IpplTimings::stopTimer(IntegrationTimer_m);

    // Destroy particles if they are marked as Bin = -1 in the plugin elements
    // or out of the global aperture
    deleteParticle();
}


void ParallelCyclotronTracker::bunchMode_m(double& t, const double dt, bool& dumpEachTurn) {

    // Flag for transition from single-bunch to multi-bunches mode
    static bool flagTransition = false;

    // single particle dumping
    if(step_m % Options::sptDumpFreq == 0)
        singleParticleDump();

    // Find out if we need to do bunch injection
    if (numBunch_m > 1)
        injectBunch_m(flagTransition);

//     oldReferenceTheta = calculateAngle2(meanR(0), meanR(1));

    // Calculate SC field before each time step and keep constant during integration.
    // Space Charge effects are included only when total macropaticles number is NOT LESS THAN 1000.
    if (itsBunch_m->hasFieldSolver()) {

        if (step_m % setup_m.scSolveFreq == 0) {

            IpplTimings::startTimer(TransformTimer_m);

            computeSpaceChargeFields_m();

            IpplTimings::stopTimer(TransformTimer_m);

        } else {
            // If we are not solving for the space charge fields at this time step
            // we will apply the fields from the previous step and have to rotate them
            // accordingly. For this we find the quaternion between the previous mean momentum (PreviousMeanP)
            // and the current mean momentum (meanP) and rotate the fields with this quaternion.

            Vector_t const meanP = calcMeanP();

            Quaternion_t quaternionToNewMeanP;

            getQuaternionTwoVectors(PreviousMeanP, meanP, quaternionToNewMeanP);

            // Reset PreviousMeanP. Cave: This HAS to be after the quaternion is calculated!
            PreviousMeanP = calcMeanP();

            // Rotate the fields
            globalToLocal(itsBunch_m->Ef, quaternionToNewMeanP);
            globalToLocal(itsBunch_m->Bf, quaternionToNewMeanP);
        }
    }

    // *** This was moved here b/c collision should be tested before the **********************
    // *** actual timestep (bgf_main_collision_test() predicts the next step automatically) -DW
    // Apply the plugin elements: probe, collimator, stripper, septum
    applyPluginElements(dt);

    // check if we lose particles at the boundary
    bgf_main_collision_test();

    IpplTimings::startTimer(IntegrationTimer_m);

    for(size_t i = 0; i < itsBunch_m->getLocalNum(); i++) {

        // used for gap crossing checking
        Vector_t Rold = itsBunch_m->R[i]; // [x,y,z]    (mm)
        Vector_t Pold = itsBunch_m->P[i]; // [px,py,pz] (beta*gamma)

        // Integrate for one step in the lab Cartesian frame (absolute value).
        itsStepper_mp->advance(itsBunch_m, i, t, dt);

        // If gap crossing happens, do momenta kicking (not if gap crossing just happened)
        if (itsBunch_m->cavityGapCrossed[i] == true)
            itsBunch_m->cavityGapCrossed[i] = false;
        else
            gapCrossKick_m(i, t, dt, Rold, Pold);
    }

    IpplTimings::stopTimer(IntegrationTimer_m);

    // Destroy particles if they are marked as Bin = -1 in the plugin elements
    // or out of global aperture
    deleteParticle();

    // Recalculate bingamma and reset the BinID for each particles according to its current gamma
    if (itsBunch_m->weHaveBins() && BunchCount_m > 1 && step_m % Options::rebinFreq == 0)
    {
        updateParticleBins_m();
    }

    // Some status output for user after each turn
    if ( (step_m > 10) && ((step_m + 1) % setup_m.stepsPerTurn) == 0) {
        turnnumber_m++;
        dumpEachTurn = true;

        *gmsg << endl;
        *gmsg << "*** Finished turn " << turnnumber_m - 1
              << ", Total number of live particles: "
              << itsBunch_m->getTotalNum() << endl;
    }

    Ippl::Comm->barrier();
}


void ParallelCyclotronTracker::gapCrossKick_m(size_t i, double t,
                                              double dt,
                                              const Vector_t& Rold,
                                              const Vector_t& Pold)
{
    for (beamline_list::iterator sindex = ++(FieldDimensions.begin());
        sindex != FieldDimensions.end(); ++sindex)
    {
        bool tag_crossing = false;
        double DistOld = 0.0; //mm
        RFCavity * rfcav;

        if (((*sindex)->first) == ElementBase::RFCAVITY) {
            // here check gap cross in the list, if do , set tag_crossing to TRUE
            rfcav = static_cast<RFCavity *>(((*sindex)->second).second);
            tag_crossing = checkGapCross(Rold, itsBunch_m->R[i], rfcav, DistOld);
        }

        if ( tag_crossing ) {
            itsBunch_m->cavityGapCrossed[i] = true;

            double oldMomentum2  = dot(Pold, Pold);
            double oldBetgam = sqrt(oldMomentum2);
            double oldGamma = sqrt(1.0 + oldMomentum2);
            double oldBeta = oldBetgam / oldGamma;
            double dt1 = DistOld / (Physics::c * oldBeta * 1.0e-6); // ns
            double dt2 = dt - dt1;

            // retrack particle from the old postion to cavity gap point
            // restore the old coordinates and momenta
            itsBunch_m->R[i] = Rold;
            itsBunch_m->P[i] = Pold;

            if (dt / dt1 < 1.0e9)
                itsStepper_mp->advance(itsBunch_m, i, t, dt1);

            // Momentum kick
            RFkick(rfcav, t, dt1, i);

            /* Retrack particle from cavity gap point for
             * the left time to finish the entire timestep
             */
            if (dt / dt2 < 1.0e9)
                itsStepper_mp->advance(itsBunch_m, i, t, dt2);
        }
    }
}


void ParallelCyclotronTracker::dumpAzimuthAngles_m(const double& t,
                                                   const Vector_t& R,
                                                   const Vector_t& P,
                                                   const double& oldReferenceTheta,
                                                   const double& temp_meanTheta)
{
    for (unsigned int i=0; i<=2; i++) {
        if ((oldReferenceTheta < azimuth_angle_m[i] - setup_m.deltaTheta) &&
            (  temp_meanTheta >= azimuth_angle_m[i] - setup_m.deltaTheta)) {
            *(outfTheta_m[i]) << "#Turn number = " << turnnumber_m
                              << ", Time = " << t << " [ns]" << std::endl
                              << " " << std::hypot(R(0), R(1))
                              << " " << P(0) * std::cos(temp_meanTheta) + P(1) * std::sin(temp_meanTheta)
                              << " " << temp_meanTheta * Physics::rad2deg
                              << " " << - P(0) * std::sin(temp_meanTheta) + P(1) * std::cos(temp_meanTheta)
                              << " " << R(2)
                              << " " << P(2) << std::endl;
            }
        }
}

void ParallelCyclotronTracker::dumpThetaEachTurn_m(const double& t,
                                                   const Vector_t& R,
                                                   const Vector_t& P,
                                                   const double& temp_meanTheta,
                                                   bool& dumpEachTurn)
{
    if ((step_m > 10) && ((step_m + 1) % setup_m.stepsPerTurn) == 0) {

        ++turnnumber_m;

        dumpEachTurn = true;

        *gmsg << "* SPT: Finished turn " << turnnumber_m - 1 << endl;

        *(outfTheta_m[3]) << "#Turn number = " << turnnumber_m
                          << ", Time = " << t << " [ns]" << std::endl
                          << " " << std::sqrt(R(0) * R(0) + R(1) * R(1))
                          << " " << P(0) * std::cos(temp_meanTheta) +
                                    P(1) * std::sin(temp_meanTheta)
                          << " " << temp_meanTheta / Physics::deg2rad
                          << " " << - P(0) * std::sin(temp_meanTheta) +
                                      P(1) * std::cos(temp_meanTheta)
                          << " " << R(2)
                          << " " << P(2) << std::endl;
    }
}


void ParallelCyclotronTracker::computeSpaceChargeFields_m() {
    // Firstly reset E and B to zero before fill new space charge field data for each track step
    itsBunch_m->Bf = Vector_t(0.0);
    itsBunch_m->Ef = Vector_t(0.0);

    if (spiral_flag && itsBunch_m->getFieldSolverType() == "SAAMG") {
        // --- Single bunch mode with spiral inflector --- //

        // If we are doing a spiral inflector simulation and are using the SAAMG solver
        // we don't rotate or translate the bunch and gamma is 1.0 (non-relativistic).

        //itsBunch_m->R *= Vector_t(0.001); // mm --> m

        IpplTimings::stopTimer(TransformTimer_m);

        itsBunch_m->setGlobalMeanR(Vector_t(0.0, 0.0, 0.0));
        itsBunch_m->setGlobalToLocalQuaternion(Quaternion_t(1.0, 0.0, 0.0, 0.0));

        itsBunch_m->computeSelfFields_cycl(1.0);

        IpplTimings::startTimer(TransformTimer_m);

        //itsBunch_m->R *= Vector_t(1000.0); // m --> mm

    } else {

        Vector_t const meanR = calcMeanR(); // (m)

        PreviousMeanP = calcMeanP();

        // Since calcMeanP takes into account all particles of all bins (TODO: Check this! -DW)
        // Using the quaternion method with PreviousMeanP and yaxis should give the correct result

        Quaternion_t quaternionToYAxis;

        getQuaternionTwoVectors(PreviousMeanP, yaxis, quaternionToYAxis);

#if COORD_PRINT
        if ( BunchCount_m > 1 ) {
            std::stringstream ss;
            ss << "global-coordinates-step-" << step_m;
            std::ofstream out(ss.str());
            for (uint i = 0; i < itsBunch_m->getLocalNum(); ++i) {
                out << itsBunch_m->R[i](0) << " "
                    << itsBunch_m->R[i](1) << " "
                    << itsBunch_m->R[i](2) << std::endl;
            }
            out.close();
        }
#endif

        globalToLocal(itsBunch_m->R, quaternionToYAxis, meanR);

#if COORD_PRINT
        if ( BunchCount_m > 1 ) {
            std::stringstream ss;
            ss << "local-coordinates-step-" << step_m;
            std::ofstream out(ss.str());
            for (uint i = 0; i < itsBunch_m->getLocalNum(); ++i) {
                out << itsBunch_m->R[i](0) << " "
                    << itsBunch_m->R[i](1) << " "
                    << itsBunch_m->R[i](2) << std::endl;
            }
            out.close();
        }
#endif

        //itsBunch_m->R *= Vector_t(0.001); // mm --> m

        if ((step_m + 1) % Options::boundpDestroyFreq == 0)
            itsBunch_m->boundp_destroy();
        else
            itsBunch_m->boundp();

        IpplTimings::stopTimer(TransformTimer_m);

        if ((itsBunch_m->weHaveBins()) && BunchCount_m > 1) {
            // --- Multibunche mode --- //

            // Calculate gamma for each energy bin
            itsBunch_m->calcGammas_cycl();

            repartition();

            // Calculate space charge field for each energy bin
            for(int b = 0; b < itsBunch_m->getLastemittedBin(); b++) {

                itsBunch_m->setBinCharge(b, itsBunch_m->getChargePerParticle());
                //itsBunch_m->setGlobalMeanR(0.001 * meanR);
                itsBunch_m->setGlobalMeanR(meanR);
                itsBunch_m->setGlobalToLocalQuaternion(quaternionToYAxis);
                itsBunch_m->computeSelfFields_cycl(b);
            }

            itsBunch_m->Q = itsBunch_m->getChargePerParticle();

        } else {
            // --- Single bunch mode --- //

            double temp_meangamma = std::sqrt(1.0 + dot(PreviousMeanP, PreviousMeanP));

            repartition();

            //itsBunch_m->setGlobalMeanR(0.001 * meanR);
            itsBunch_m->setGlobalMeanR(meanR);
            itsBunch_m->setGlobalToLocalQuaternion(quaternionToYAxis);

            itsBunch_m->computeSelfFields_cycl(temp_meangamma);

        }


        IpplTimings::startTimer(TransformTimer_m);

        //scale coordinates back
        //itsBunch_m->R *= Vector_t(1000.0); // m --> mm

        // Transform coordinates back to global
        localToGlobal(itsBunch_m->R, quaternionToYAxis, meanR);

        // Transform self field back to global frame (rotate only)
        localToGlobal(itsBunch_m->Ef, quaternionToYAxis);
        localToGlobal(itsBunch_m->Bf, quaternionToYAxis);
    }
}


bool ParallelCyclotronTracker::computeExternalFields_m(const size_t& i, const double& t,
                                                       Vector_t& Efield, Vector_t& Bfield)
{
    beamline_list::iterator sindex = FieldDimensions.begin();

    // Flag whether a particle is out of field
    bool outOfBound = (((*sindex)->second).second)->apply(i, t, Efield, Bfield);

    Bfield *= 0.10;  // kGauss --> T
    Efield *= 1.0e6; // kV/mm  --> V/m

    return outOfBound;
}


void ParallelCyclotronTracker::injectBunch_m(bool& flagTransition) {
    if ((BunchCount_m == 1) && (multiBunchMode_m == MB_MODE::AUTO) && (!flagTransition)) {

        if (step_m == setup_m.stepsNextCheck) {
            // If all of the following conditions are met, this code will be executed
            // to check the distance between two neighboring bunches:
            // 1. Only one bunch exists (BunchCount_m == 1)
            // 2. We are in multi-bunch mode, AUTO sub-mode (multiBunchMode_m == 2)
            // 3. It has been a full revolution since the last check (stepsNextCheck)

            *gmsg << "* MBM: Checking for automatically injecting new bunch ..." << endl;

            //itsBunch_m->R *= Vector_t(0.001); // mm --> m
            itsBunch_m->calcBeamParameters();
            //itsBunch_m->R *= Vector_t(1000.0); // m --> mm

            Vector_t Rmean = itsBunch_m->get_centroid(); // m

            RThisTurn_m = std::hypot(Rmean[0],Rmean[1]);

            Vector_t Rrms = itsBunch_m->get_rrms(); // m

            double XYrms = std::hypot(Rrms[0], Rrms[1]);

            // If the distance between two neighboring bunches is less than 5 times of its 2D rms size
            // start multi-bunch simulation, fill current phase space to initialR and initialP arrays
            if ((RThisTurn_m - RLastTurn_m) < CoeffDBunches_m * XYrms) {
                // since next turn, start multi-bunches
                saveOneBunch();
                flagTransition = true;
                *gmsg << "* MBM: Saving beam distribution at turn " << turnnumber_m << endl;
                *gmsg << "* MBM: After one revolution, Multi-Bunch Mode will be invoked" << endl;
            }

            setup_m.stepsNextCheck += setup_m.stepsPerTurn;

            *gmsg << "* MBM: RLastTurn = " << RLastTurn_m << " [m]" << endl;
            *gmsg << "* MBM: RThisTurn = " << RThisTurn_m << " [m]" << endl;
            *gmsg << "* MBM: XYrms = " << XYrms    << " [m]" << endl;

            RLastTurn_m = RThisTurn_m;
        }
    }

    else if ((BunchCount_m < numBunch_m) && (step_m == setup_m.stepsNextCheck)) {
        // Matthias: SteptoLastInj was used in MtsTracker, removed by DW in GenericTracker

        // If all of the following conditions are met, this code will be executed
        // to read new bunch from hdf5 format file:
        // 1. We are in multi-bunch mode (numBunch_m > 1)
        // 2. It has been a full revolution since the last check
        // 3. Number of existing bunches is less than the desired number of bunches
        // 4. FORCE mode, or AUTO mode with flagTransition = true
        // Note: restart from 1 < BunchCount < numBunch_m must be avoided.
        *gmsg << "* MBM: Step " << step_m << ", injecting a new bunch ..." << endl;

        BunchCount_m++;

        itsBunch_m->setNumBunch(BunchCount_m);

        // read initial distribution from h5 file
        switch ( multiBunchMode_m ) {
            case MB_MODE::FORCE:
            case MB_MODE::AUTO:
                readOneBunchFromFile(BunchCount_m - 1);
                updateParticleBins_m();
                break;
            case MB_MODE::NONE:
                // do nothing
            default:
                throw OpalException("ParallelCyclotronTracker::injectBunch_m()",
                                    "We shouldn't be here in single bunch mode.");
        }

        setup_m.stepsNextCheck += setup_m.stepsPerTurn;

        Ippl::Comm->barrier();

        *gmsg << "* MBM: Bunch " << BunchCount_m
              << " injected, total particle number = "
              << itsBunch_m->getTotalNum() << endl;

    } else if (BunchCount_m == numBunch_m) {
        // After this, numBunch_m is wrong but not needed anymore...
        numBunch_m--;
    }
}