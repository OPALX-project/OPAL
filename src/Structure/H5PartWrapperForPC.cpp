//
//  Copyright & License: See Copyright.readme in src directory
//

#include "Structure/H5PartWrapperForPC.h"

#include "OPALconfig.h"
#include "Algorithms/PartBunchBase.h"
#include "AbstractObjects/OpalData.h"
#include "Utilities/Options.h"
#include "Utilities/Util.h"
#include "Physics/Physics.h"

#include <boost/filesystem.hpp>

#include <sstream>
#include <set>

H5PartWrapperForPC::H5PartWrapperForPC(const std::string &fileName, h5_int32_t flags):
    H5PartWrapper(fileName, flags),
    previousH5Local_m(false),
    referenceMomentum_m(0.0),
    referenceLocation_m(0.0),
    meanE_m(0.0),
    meanMomentum_m(0.0),
    azimuth_m(0.0),
    elevation_m(0.0)
{ }

H5PartWrapperForPC::H5PartWrapperForPC(const std::string &fileName,
                                       int restartStep,
                                       std::string sourceFile,
                                       h5_int32_t flags):
    H5PartWrapper(fileName, restartStep, sourceFile, flags),
    previousH5Local_m(false),
    referenceMomentum_m(0.0),
    referenceLocation_m(0.0),
    meanE_m(0.0),
    meanMomentum_m(0.0),
    azimuth_m(0.0),
    elevation_m(0.0)
{ }

H5PartWrapperForPC::~H5PartWrapperForPC()
{ }

void H5PartWrapperForPC::readHeader() {
    h5_int64_t numFileAttributes = H5GetNumFileAttribs(file_m);

    const h5_size_t lengthAttributeName = 256;
    char attributeName[lengthAttributeName];
    h5_int64_t attributeType;
    h5_size_t numAttributeElements;
    std::set<std::string> attributeNames;

    for (h5_int64_t i = 0; i < numFileAttributes; ++ i) {
        REPORTONERROR(H5GetFileAttribInfo(file_m,
                                          i,
                                          attributeName,
                                          lengthAttributeName,
                                          &attributeType,
                                          &numAttributeElements));

        attributeNames.insert(attributeName);
    }

    if (attributeNames.find("dump frequency") != attributeNames.end()) {
        h5_int64_t dumpfreq;
        READFILEATTRIB(Int64, file_m, "dump frequency", &dumpfreq);
        OpalData::getInstance()->setRestartDumpFreq(dumpfreq);
    }
}

void H5PartWrapperForPC::readStep(PartBunchBase<double, 3>* bunch,
                                  h5_ssize_t firstParticle,
                                  h5_ssize_t lastParticle)
{
    h5_ssize_t numStepsInSource = H5GetNumSteps(file_m);
    h5_ssize_t readStep = numStepsInSource - 1;
    REPORTONERROR(H5SetStep(file_m, readStep));

    readStepHeader(bunch);
    readStepData(bunch, firstParticle, lastParticle);
}

void H5PartWrapperForPC::readStepHeader(PartBunchBase<double, 3>* bunch) {
    h5_float64_t pathLength;
    READSTEPATTRIB(Float64, file_m, "SPOS", &pathLength);
    bunch->setLPath(pathLength);

    h5_int64_t ltstep;
    READSTEPATTRIB(Int64, file_m, "LocalTrackStep", &ltstep);
    bunch->setLocalTrackStep((long long)ltstep);

    h5_int64_t gtstep;
    READSTEPATTRIB(Int64, file_m, "GlobalTrackStep", &gtstep);
    bunch->setGlobalTrackStep((long long)gtstep);

    READSTEPATTRIB(Float64, file_m, "ENERGY", &meanE_m);

    double actualT;
    READSTEPATTRIB(Float64, file_m, "TIME", &actualT);
    bunch->setT(actualT);

    h5_int64_t SteptoLastInj;
    READSTEPATTRIB(Int64, file_m, "SteptoLastInj", &SteptoLastInj);
    bunch->setSteptoLastInj((int)SteptoLastInj);

    h5_int64_t numBunch;
    READSTEPATTRIB(Int64, file_m, "NumBunch", &numBunch);
    bunch->setNumBunch((short)numBunch);

    if (predecessorOPALFlavour_m == "opal-cycl") {
        READSTEPATTRIB(Float64, file_m, "REFPR", (h5_float64_t*) &referenceMomentum_m[0]);
        READSTEPATTRIB(Float64, file_m, "REFPT", (h5_float64_t*) &referenceMomentum_m[1]);
        READSTEPATTRIB(Float64, file_m, "REFPZ", (h5_float64_t*) &referenceMomentum_m[2]);

        READSTEPATTRIB(Float64, file_m, "REFR", (h5_float64_t*) &referenceLocation_m[0]);
        READSTEPATTRIB(Float64, file_m, "REFTHETA", (h5_float64_t*) &referenceLocation_m[1]);
        READSTEPATTRIB(Float64, file_m, "REFZ", (h5_float64_t*) &referenceLocation_m[2]);

        READSTEPATTRIB(Float64, file_m, "AZIMUTH", &azimuth_m);
        READSTEPATTRIB(Float64, file_m, "ELEVATION", &elevation_m);

        h5_int64_t localDump = 0;
        h5_int64_t rc = H5ReadStepAttribInt64(file_m, "LOCAL", &localDump);
        if(rc != H5_SUCCESS) {

            ERRORMSG("H5 rc= " << rc << " in " << __FILE__ << " @ line " << __LINE__ << endl);

            std::string errormsg =
                "You are trying to restart from a legacy file that doesn't contain"
                "information on local/global frame. We are working on legacy support,"
                "but for now you have to use OPAL 1.3.0!";

            throw OpalException("H5PartWrapperForPC::readStepHeader", errormsg);
        } else {
            if (localDump == 1) previousH5Local_m = true;
        }
    } else {
        bunch->setT(0.0);
        bunch->setLocalTrackStep((long long) 0 );
    }

    double mass = bunch->getM() * 1e-6;
    meanMomentum_m = sqrt(std::pow(meanE_m,2.0) + 2 * meanE_m * mass) / mass;
}

void H5PartWrapperForPC::readStepData(PartBunchBase<double, 3>* bunch,
                                      h5_ssize_t firstParticle,
                                      h5_ssize_t lastParticle)
{
    h5_ssize_t numParticles = getNumParticles();
    if (lastParticle >= numParticles || firstParticle > lastParticle) {
        throw OpalException("H5PartWrapperForPC::readStepData",
                            "the provided particle numbers don't match the number of particles in the file");
    }
    const unsigned int yIndex = (predecessorOPALFlavour_m == "opal-t"? 2: 1);
    const unsigned int zIndex = (predecessorOPALFlavour_m == "opal-t"? 1: 2);

    REPORTONERROR(H5PartSetView(file_m, firstParticle, lastParticle));

    numParticles = lastParticle - firstParticle + 1;

    std::vector<char> buffer(numParticles * sizeof(h5_float64_t));
    h5_float64_t *f64buffer = reinterpret_cast<h5_float64_t*>(&buffer[0]);

    READDATA(Float64, file_m, "x", f64buffer);
    for(long int n = 0; n < numParticles; ++ n) {
        bunch->R[n](0) = f64buffer[n];
        bunch->Bin[n] = 0;
    }

    READDATA(Float64, file_m, "y", f64buffer);
    for(long int n = 0; n < numParticles; ++ n) {
        bunch->R[n](yIndex) = f64buffer[n];
    }

    READDATA(Float64, file_m, "z", f64buffer);
    for(long int n = 0; n < numParticles; ++ n) {
        bunch->R[n](zIndex) = f64buffer[n];
    }

    READDATA(Float64, file_m, "px", f64buffer);
    for(long int n = 0; n < numParticles; ++ n) {
        bunch->P[n](0) = f64buffer[n];
    }

    READDATA(Float64, file_m, "py", f64buffer);
    for(long int n = 0; n < numParticles; ++ n) {
        bunch->P[n](yIndex) = f64buffer[n];
    }

    READDATA(Float64, file_m, "pz", f64buffer);
    for(long int n = 0; n < numParticles; ++ n) {
        bunch->P[n](zIndex) = f64buffer[n];
    }

    READDATA(Float64, file_m, "q", f64buffer);
    for(long int n = 0; n < numParticles; ++ n) {
        bunch->Q[n] = f64buffer[n];
    }

    READDATA(Float64, file_m, "mass", f64buffer);
    for(long int n = 0; n < numParticles; ++ n) {
        bunch->M[n] = f64buffer[n];
    }

    if ( bunch->getNumBunch() > 1 ) {
        std::vector<char> ibuffer(numParticles * sizeof(h5_int64_t));
        h5_int64_t *i64buffer = reinterpret_cast<h5_int64_t*>(&ibuffer[0]);
        READDATA(Int64, file_m, "bin", i64buffer);
        for(long int n = 0; n < numParticles; ++ n) {
            bunch->Bin[n] = i64buffer[n];
        }
    }

    std::vector<char> ibuffer(numParticles * sizeof(h5_int64_t));
    h5_int64_t *i64buffer = reinterpret_cast<h5_int64_t*>(&ibuffer[0]);
    READDATA(Int64, file_m, "bunchNumber", i64buffer);
    for(long int n = 0; n < numParticles; ++ n) {
        bunch->bunchNum[n] = i64buffer[n];
    }

    REPORTONERROR(H5PartSetView(file_m, -1, -1));
}

void H5PartWrapperForPC::writeHeader() {
    std::stringstream OPAL_version;
    OPAL_version << OPAL_PROJECT_NAME << " " << OPAL_PROJECT_VERSION << " # git rev. " << Util::getGitRevision();
    WRITESTRINGFILEATTRIB(file_m, "OPAL_version", OPAL_version.str().c_str());

    // Units for the step data
    WRITESTRINGFILEATTRIB(file_m, "xUnit", "m");
    WRITESTRINGFILEATTRIB(file_m, "yUnit", "m");
    WRITESTRINGFILEATTRIB(file_m, "zUnit", "m");
    WRITESTRINGFILEATTRIB(file_m, "pxUnit", "#beta#gamma");
    WRITESTRINGFILEATTRIB(file_m, "pyUnit", "#beta#gamma");
    WRITESTRINGFILEATTRIB(file_m, "pzUnit", "#beta#gamma");
    WRITESTRINGFILEATTRIB(file_m, "qUnit", "C");
    WRITESTRINGFILEATTRIB(file_m, "massUnit", "GeV");
    WRITESTRINGFILEATTRIB(file_m, "idUnit", "1");
    WRITESTRINGFILEATTRIB(file_m, "binUnit", "1");
    WRITESTRINGFILEATTRIB(file_m, "bunchNumberUnit", "1");
    if (Options::ebDump) {
      WRITESTRINGFILEATTRIB(file_m, "Ex", "MV/m");
      WRITESTRINGFILEATTRIB(file_m, "Ey", "MV/m");
      WRITESTRINGFILEATTRIB(file_m, "Ez", "MV/m");
      WRITESTRINGFILEATTRIB(file_m, "Bx", "T");
      WRITESTRINGFILEATTRIB(file_m, "By", "T");
      WRITESTRINGFILEATTRIB(file_m, "Bz", "T");
    }

    // Units for the step attributes
    WRITESTRINGFILEATTRIB(file_m, "TIMEUnit",     "s");
    WRITESTRINGFILEATTRIB(file_m, "SPOSUnit",     "m");

    WRITESTRINGFILEATTRIB(file_m, "RefPartRUnit", "m");
    WRITESTRINGFILEATTRIB(file_m, "centroidUnit", "m");
    WRITESTRINGFILEATTRIB(file_m, "RMSXUnit",     "m");
    WRITESTRINGFILEATTRIB(file_m, "minXUnit",     "m");
    WRITESTRINGFILEATTRIB(file_m, "maxXUnit",     "m");

    WRITESTRINGFILEATTRIB(file_m, "RefPartPUnit", "#beta#gamma");
    WRITESTRINGFILEATTRIB(file_m, "RMSPUnit",     "#beta#gamma");
    WRITESTRINGFILEATTRIB(file_m, "minPUnit",     "#beta#gamma");
    WRITESTRINGFILEATTRIB(file_m, "maxPUnit",     "#beta#gamma");

    WRITESTRINGFILEATTRIB(file_m, "#varepsilonUnit",      "m rad");
    WRITESTRINGFILEATTRIB(file_m, "#varepsilon-geomUnit", "m rad");

    WRITESTRINGFILEATTRIB(file_m, "StepUnit", "1");
    WRITESTRINGFILEATTRIB(file_m, "LocalTrackStepUnit", "1");
    WRITESTRINGFILEATTRIB(file_m, "GlobalTrackStepUnit", "1");

    WRITESTRINGFILEATTRIB(file_m, "#sigmaUnit", "1");

    WRITESTRINGFILEATTRIB(file_m, "ENERGYUnit", "MeV");
    WRITESTRINGFILEATTRIB(file_m, "dEUnit",     "MeV");

    WRITESTRINGFILEATTRIB(file_m, "MASSUnit",   "GeV");
    WRITESTRINGFILEATTRIB(file_m, "CHARGEUnit", "C");

    WRITESTRINGFILEATTRIB(file_m, "NumBunchUnit",      "1");
    WRITESTRINGFILEATTRIB(file_m, "SteptoLastInjUnit", "1");
    WRITESTRINGFILEATTRIB(file_m, "LOCALUnit",         "1");

    // additional step attributes
    WRITESTRINGFILEATTRIB(file_m, "REFPRUnit",        "#beta#gamma");
    WRITESTRINGFILEATTRIB(file_m, "REFPTUnit",        "#beta#gamma");
    WRITESTRINGFILEATTRIB(file_m, "REFPZUnit",        "#beta#gamma");
    WRITESTRINGFILEATTRIB(file_m, "REFRUnit",         "mm");
    WRITESTRINGFILEATTRIB(file_m, "REFTHETAUnit",     "deg");
    WRITESTRINGFILEATTRIB(file_m, "REFZUnit",         "mm");
    WRITESTRINGFILEATTRIB(file_m, "REFAZIMUTHUnit",   "deg");
    WRITESTRINGFILEATTRIB(file_m, "REFELEVATIONUnit", "deg");

    WRITESTRINGFILEATTRIB(file_m, "spos-headUnit", "m");
    WRITESTRINGFILEATTRIB(file_m, "spos-refUnit",  "m");
    WRITESTRINGFILEATTRIB(file_m, "spos-tailUnit", "m");

    WRITESTRINGFILEATTRIB(file_m, "B-refUnit",  "T");
    WRITESTRINGFILEATTRIB(file_m, "E-refUnit",  "MV/m");
    WRITESTRINGFILEATTRIB(file_m, "B-headUnit", "T");
    WRITESTRINGFILEATTRIB(file_m, "E-headUnit", "MV/m");
    WRITESTRINGFILEATTRIB(file_m, "B-tailUnit", "T");
    WRITESTRINGFILEATTRIB(file_m, "E-tailUnit", "MV/m");

    /// Write file dump frequency.
    h5_int64_t dumpfreq = Options::psDumpFreq;
    WRITEFILEATTRIB(Int64, file_m, "dump frequency", &dumpfreq, 1);
}

void H5PartWrapperForPC::writeStep(PartBunchBase<double, 3>* bunch,
                                   const std::map<std::string, double> &additionalStepAttributes)
{
    if (bunch->getTotalNum() == 0) return;

    open(H5_O_APPENDONLY);
    writeStepHeader(bunch, additionalStepAttributes);
    writeStepData(bunch);
    close();

    ++ numSteps_m;
}

void H5PartWrapperForPC::writeStepHeader(PartBunchBase<double, 3>* bunch,
                                         const std::map<std::string, double> &additionalStepAttributes)
{
    bunch->calcBeamParameters();

    double   t          = bunch->getT();
    double   pathLength = bunch->getLPath();
    Vector_t rmin       = bunch->get_origin();
    Vector_t rmax       = bunch->get_maxExtent();
    Vector_t centroid   = bunch->get_centroid();

    Vector_t meanR = bunch->get_rmean();
    Vector_t meanP = bunch->get_pmean();
    Vector_t xsigma = bunch->get_rrms();
    Vector_t psigma = bunch->get_prms();
    Vector_t vareps = bunch->get_norm_emit();
    Vector_t geomvareps = bunch->get_emit();

    Vector_t RefPartR = bunch->RefPartR_m;
    Vector_t RefPartP = bunch->RefPartP_m;

    double meanEnergy = bunch->get_meanKineticEnergy();
    double energySpread = bunch->getdE();
    double I_0 = 4.0 * Physics::pi * Physics::epsilon_0 * Physics::c * bunch->getM() / bunch->getQ();
    double sigma = ((xsigma[0] * xsigma[0]) + (xsigma[1] * xsigma[1])) /
        (2.0 * bunch->get_gamma() * I_0 * (geomvareps[0] * geomvareps[0] + geomvareps[1] * geomvareps[1]));

    h5_int64_t localTrackStep = (h5_int64_t)bunch->getLocalTrackStep();
    h5_int64_t globalTrackStep = (h5_int64_t)bunch->getGlobalTrackStep();

    h5_int64_t numBunch = (h5_int64_t)bunch->getNumBunch();
    h5_int64_t SteptoLastInj = (h5_int64_t)bunch->getSteptoLastInj();

    double mass = 1.0e-9 * bunch->getM();
    double charge = bunch->getCharge();

    h5_int64_t localFrame = ( Options::psDumpFrame != Options::GLOBAL );

    double sposHead = 0.0;
    double sposRef = 0.0;
    double sposTail = 0.0;

    Vector_t maxP(0.0);
    Vector_t minP(0.0);
    bunch->get_PBounds(minP, maxP);

    REPORTONERROR(H5SetStep(file_m, numSteps_m));

    char const *OPALFlavour = "opal-cycl";
    WRITESTRINGSTEPATTRIB(file_m, "OPAL_flavour", OPALFlavour);
    WRITESTEPATTRIB(Float64, file_m, "TIME", &t, 1);
    WRITESTEPATTRIB(Float64, file_m, "SPOS", &pathLength, 1);

    WRITESTEPATTRIB(Float64, file_m, "RefPartR", (h5_float64_t *)&RefPartR, 3);
    WRITESTEPATTRIB(Float64, file_m, "centroid", (h5_float64_t *)&centroid, 3);
    WRITESTEPATTRIB(Float64, file_m, "RMSX",     (h5_float64_t *)&xsigma,   3);
    WRITESTEPATTRIB(Float64, file_m, "minX",     (h5_float64_t *)&rmin,     3);
    WRITESTEPATTRIB(Float64, file_m, "maxX",     (h5_float64_t *)&rmax,     3);

    WRITESTEPATTRIB(Float64, file_m, "RefPartP", (h5_float64_t *)&RefPartP, 3);
    WRITESTEPATTRIB(Float64, file_m, "RMSP",     (h5_float64_t *)&psigma,   3);
    WRITESTEPATTRIB(Float64, file_m, "minP",     (h5_float64_t *)&minP,     3);
    WRITESTEPATTRIB(Float64, file_m, "maxP",     (h5_float64_t *)&maxP,     3);

    WRITESTEPATTRIB(Float64, file_m, "#varepsilon",      (h5_float64_t *)&vareps,     3);
    WRITESTEPATTRIB(Float64, file_m, "#varepsilon-geom", (h5_float64_t *)&geomvareps, 3);

    WRITESTEPATTRIB(Int64, file_m, "Step",            &numSteps_m,      1);
    WRITESTEPATTRIB(Int64, file_m, "LocalTrackStep",  &localTrackStep,  1);
    WRITESTEPATTRIB(Int64, file_m, "GlobalTrackStep", &globalTrackStep, 1);

    WRITESTEPATTRIB(Float64, file_m, "#sigma", &sigma, 1);

    WRITESTEPATTRIB(Float64, file_m, "ENERGY", &meanEnergy, 1);
    WRITESTEPATTRIB(Float64, file_m, "dE", &energySpread, 1);

    /// Write particle mass and charge per particle. (Consider making these file attributes.)
    WRITESTEPATTRIB(Float64, file_m, "MASS", &mass, 1);
    WRITESTEPATTRIB(Float64, file_m, "CHARGE", &charge, 1);

    WRITESTEPATTRIB(Int64, file_m, "NumBunch", &numBunch, 1);
    WRITESTEPATTRIB(Int64, file_m, "SteptoLastInj", &SteptoLastInj, 1);
    WRITESTEPATTRIB(Int64, file_m, "LOCAL", &localFrame, 1);

    try {
        h5_float64_t refpr     = additionalStepAttributes.at("REFPR");
        h5_float64_t refpt     = additionalStepAttributes.at("REFPT");
        h5_float64_t refpz     = additionalStepAttributes.at("REFPZ");
        h5_float64_t refr      = additionalStepAttributes.at("REFR");
        h5_float64_t reft      = additionalStepAttributes.at("REFTHETA");
        h5_float64_t refz      = additionalStepAttributes.at("REFZ");
        h5_float64_t azimuth   = additionalStepAttributes.at("AZIMUTH");
        h5_float64_t elevation = additionalStepAttributes.at("ELEVATION");

        Vector_t referenceB(additionalStepAttributes.at("B-ref_x"),
                            additionalStepAttributes.at("B-ref_z"),
                            additionalStepAttributes.at("B-ref_y"));
        Vector_t referenceE(additionalStepAttributes.at("E-ref_x"),
                            additionalStepAttributes.at("E-ref_z"),
                            additionalStepAttributes.at("E-ref_y"));
        Vector_t headB(additionalStepAttributes.at("B-head_x"),
                      additionalStepAttributes.at("B-head_z"),
                      additionalStepAttributes.at("B-head_y"));
        Vector_t headE(additionalStepAttributes.at("E-head_x"),
                       additionalStepAttributes.at("E-head_z"),
                       additionalStepAttributes.at("E-head_y"));
        Vector_t tailB(additionalStepAttributes.at("B-tail_x"),
                      additionalStepAttributes.at("B-tail_z"),
                      additionalStepAttributes.at("B-tail_y"));
        Vector_t tailE(additionalStepAttributes.at("E-tail_x"),
                       additionalStepAttributes.at("E-tail_z"),
                       additionalStepAttributes.at("E-tail_y"));

        WRITESTEPATTRIB(Float64, file_m, "REFPR",     &refpr, 1);
        WRITESTEPATTRIB(Float64, file_m, "REFPT",     &refpt, 1);
        WRITESTEPATTRIB(Float64, file_m, "REFPZ",     &refpz, 1);
        WRITESTEPATTRIB(Float64, file_m, "REFR",      &refr, 1);
        WRITESTEPATTRIB(Float64, file_m, "REFTHETA",  &reft, 1);
        WRITESTEPATTRIB(Float64, file_m, "REFZ",      &refz, 1);
        WRITESTEPATTRIB(Float64, file_m, "AZIMUTH",   &azimuth, 1);
        WRITESTEPATTRIB(Float64, file_m, "ELEVATION", &elevation, 1);

        WRITESTEPATTRIB(Float64, file_m, "spos-head", &sposHead, 1);
        WRITESTEPATTRIB(Float64, file_m, "spos-ref",  &sposRef,  1);
        WRITESTEPATTRIB(Float64, file_m, "spos-tail", &sposTail, 1);

        WRITESTEPATTRIB(Float64, file_m, "B-ref",  (h5_float64_t *)&referenceB, 3);
        WRITESTEPATTRIB(Float64, file_m, "E-ref",  (h5_float64_t *)&referenceE, 3);
        WRITESTEPATTRIB(Float64, file_m, "B-head", (h5_float64_t *)&headB, 3);
        WRITESTEPATTRIB(Float64, file_m, "E-head", (h5_float64_t *)&headE, 3);
        WRITESTEPATTRIB(Float64, file_m, "B-tail", (h5_float64_t *)&tailB, 3);
        WRITESTEPATTRIB(Float64, file_m, "E-tail", (h5_float64_t *)&tailE, 3);
    } catch (std::out_of_range & m) {
        ERRORMSG(m.what() << endl);

        throw OpalException("H5PartWrapperForPC::writeStepHeader",
                            "some additional step attribute not found");
    }
}

void H5PartWrapperForPC::writeStepData(PartBunchBase<double, 3>* bunch) {
    /*
      find particle with ID==0
      and save index in zID

      FIXME After issue 287 is resolved this shouldn't be
            necessary anymore!
     */

    size_t IDZero = bunch->getLocalNum();
    bool found = false;
    for(size_t k = 0; k < IDZero; ++ k) {
        if (bunch->ID[k] == 0) {
            found = true;
            IDZero = k;
        }
    }

    const size_t numLocalParticles = (found? bunch->getLocalNum() - 1: bunch->getLocalNum());
    const size_t skipID = IDZero;

    std::vector<char> buffer(numLocalParticles * sizeof(h5_float64_t));
    h5_float64_t *f64buffer = reinterpret_cast<h5_float64_t*>(&buffer[0]);
    h5_int64_t   *i64buffer = reinterpret_cast<h5_int64_t*>  (&buffer[0]);


    REPORTONERROR(H5PartSetNumParticles(file_m, numLocalParticles));

    for(size_t i = 0; i < skipID; ++ i)
        f64buffer[i] =  bunch->R[i](0);
    for (size_t i = skipID; i < numLocalParticles; ++ i)
        f64buffer[i] = bunch->R[i + 1](0);

    WRITEDATA(Float64, file_m, "x", f64buffer);

    for(size_t i = 0; i < skipID; ++ i)
        f64buffer[i] =  bunch->R[i](1);
    for (size_t i = skipID; i < numLocalParticles; ++ i)
        f64buffer[i] = bunch->R[i + 1](1);

    WRITEDATA(Float64, file_m, "y", f64buffer);

    for(size_t i = 0; i < skipID; ++ i)
        f64buffer[i] =  bunch->R[i](2);
    for (size_t i = skipID; i < numLocalParticles; ++ i)
        f64buffer[i] = bunch->R[i + 1](2);

    WRITEDATA(Float64, file_m, "z", f64buffer);

    for(size_t i = 0; i < skipID; ++ i)
        f64buffer[i] =  bunch->P[i](0);
    for (size_t i = skipID; i < numLocalParticles; ++ i)
        f64buffer[i] = bunch->P[i + 1](0);

    WRITEDATA(Float64, file_m, "px", f64buffer);

    for(size_t i = 0; i < skipID; ++ i)
        f64buffer[i] =  bunch->P[i](1);
    for (size_t i = skipID; i < numLocalParticles; ++ i)
        f64buffer[i] = bunch->P[i + 1](1);

    WRITEDATA(Float64, file_m, "py", f64buffer);

    for(size_t i = 0; i < skipID; ++ i)
        f64buffer[i] =  bunch->P[i](2);
    for (size_t i = skipID; i < numLocalParticles; ++ i)
        f64buffer[i] = bunch->P[i + 1](2);

    WRITEDATA(Float64, file_m, "pz", f64buffer);

    for(size_t i = 0; i < skipID; ++ i)
        f64buffer[i] =  bunch->Q[i];
    for (size_t i = skipID; i < numLocalParticles; ++ i)
        f64buffer[i] = bunch->Q[i + 1];

    WRITEDATA(Float64, file_m, "q", f64buffer);

    for(size_t i = 0; i < skipID; ++ i)
        f64buffer[i] =  bunch->M[i];
    for (size_t i = skipID; i < numLocalParticles; ++ i)
        f64buffer[i] = bunch->M[i + 1];

    WRITEDATA(Float64, file_m, "mass", f64buffer);

    for(size_t i = 0; i < skipID; ++ i)
        i64buffer[i] =  bunch->ID[i];
    for (size_t i = skipID; i < numLocalParticles; ++ i)
        i64buffer[i] = bunch->ID[i + 1];

    WRITEDATA(Int64, file_m, "id", i64buffer);

    if ( bunch->hasBinning() ) {
        for(size_t i = 0; i < skipID; ++ i)
            i64buffer[i] =  bunch->Bin[i];
        for (size_t i = skipID; i < numLocalParticles; ++ i)
            i64buffer[i] = bunch->Bin[i + 1];

        WRITEDATA(Int64, file_m, "bin", i64buffer);
    }

    for  (size_t i = 0; i < skipID; ++ i)
        i64buffer[i] = bunch->bunchNum[i];
    for (size_t i = skipID; i < numLocalParticles; ++ i)
            i64buffer[i] = bunch->bunchNum[i + 1];
    WRITEDATA(Int64, file_m, "bunchNumber", i64buffer);

    if (Options::ebDump) {
        for(size_t i = 0; i < skipID; ++ i)
            f64buffer[i] =  bunch->Ef[i](0);
        for (size_t i = skipID; i < numLocalParticles; ++ i)
            f64buffer[i] = bunch->Ef[i + 1](0);

        WRITEDATA(Float64, file_m, "Ex", f64buffer);

        for(size_t i = 0; i < skipID; ++ i)
            f64buffer[i] =  bunch->Ef[i](1);
        for (size_t i = skipID; i < numLocalParticles; ++ i)
            f64buffer[i] = bunch->Ef[i + 1](1);

        WRITEDATA(Float64, file_m, "Ey", f64buffer);

        for(size_t i = 0; i < skipID; ++ i)
            f64buffer[i] =  bunch->Ef[i](2);
        for (size_t i = skipID; i < numLocalParticles; ++ i)
            f64buffer[i] = bunch->Ef[i + 1](2);

        WRITEDATA(Float64, file_m, "Ez", f64buffer);

        for(size_t i = 0; i < skipID; ++ i)
            f64buffer[i] =  bunch->Bf[i](0);
        for (size_t i = skipID; i < numLocalParticles; ++ i)
            f64buffer[i] = bunch->Bf[i + 1](0);

        WRITEDATA(Float64, file_m, "Bx", f64buffer);

        for(size_t i = 0; i < skipID; ++ i)
            f64buffer[i] =  bunch->Bf[i](1);
        for (size_t i = skipID; i < numLocalParticles; ++ i)
            f64buffer[i] = bunch->Bf[i + 1](1);

        WRITEDATA(Float64, file_m, "By", f64buffer);

        for(size_t i = 0; i < skipID; ++ i)
            f64buffer[i] =  bunch->Bf[i](2);
        for (size_t i = skipID; i < numLocalParticles; ++ i)
            f64buffer[i] = bunch->Bf[i + 1](2);

        WRITEDATA(Float64, file_m, "Bz", f64buffer);

    }

    /// Write space charge field map if asked for.
    if(Options::rhoDump) {
        NDIndex<3> idx = bunch->getFieldLayout().getLocalNDIndex();
        NDIndex<3> elem;
        REPORTONERROR(H5Block3dSetView(file_m,
                                       idx[0].min(), idx[0].max(),
                                       idx[1].min(), idx[1].max(),
                                       idx[2].min(), idx[2].max()));

        std::unique_ptr<h5_float64_t[]> data(new h5_float64_t[(idx[0].max() + 1)  * (idx[1].max() + 1) * (idx[2].max() + 1)]);

        int ii = 0;
        // h5block uses the fortran convention of storing data:
        // INTEGER, DIMENSION(2,3) :: a
        // => {a(1,1), a(2,1), a(1,2), a(2,2), a(1,3), a(2,3)}
        for(int i = idx[2].min(); i <= idx[2].max(); ++ i) {
            for(int j = idx[1].min(); j <= idx[1].max(); ++ j) {
                for(int k = idx[0].min(); k <= idx[0].max(); ++ k) {
                    data[ii] = bunch->getRho(k, j, i);
                    ++ ii;
                }
            }
        }
        REPORTONERROR(H5Block3dWriteScalarFieldFloat64(file_m, "rho", data.get()));

        /// Need this to align particles and fields when writing space charge map.
        REPORTONERROR(H5Block3dSetFieldOrigin(file_m, "rho",
                                              (h5_float64_t)bunch->get_origin()(0),
                                              (h5_float64_t)bunch->get_origin()(1),
                                              (h5_float64_t)bunch->get_origin()(2)));

        REPORTONERROR(H5Block3dSetFieldSpacing(file_m, "rho",
                                               (h5_float64_t)bunch->get_hr()(0),
                                               (h5_float64_t)bunch->get_hr()(1),
                                               (h5_float64_t)bunch->get_hr()(2)));
    }
}
