//
//  Copyright & License: See Copyright.readme in src directory
//

#include "Structure/H5PartWrapperForPS.h"

#include "OPALconfig.h"
#include "Algorithms/bet/EnvelopeBunch.h"
#include "AbstractObjects/OpalData.h"
#include "Utilities/Options.h"
#include "Utilities/Util.h"
#include "Physics/Physics.h"

#include <boost/filesystem.hpp>

#include <sstream>
#include <set>

H5PartWrapperForPS::H5PartWrapperForPS(const std::string &fileName, h5_int32_t flags):
    H5PartWrapper(fileName, flags)
{ }

H5PartWrapperForPS::H5PartWrapperForPS(const std::string &fileName, int restartStep, std::string sourceFile, h5_int32_t flags):
    H5PartWrapper(fileName, restartStep, sourceFile, flags)
{
    if (restartStep == -1) {
        restartStep = H5GetNumSteps(file_m) - 1 ;
        OpalData::getInstance()->setRestartStep(restartStep);
    }
}

H5PartWrapperForPS::~H5PartWrapperForPS()
{ }

void H5PartWrapperForPS::readHeader() {
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

    if (attributeNames.find("dPhiGlobal") != attributeNames.end()) {
        h5_float64_t dphi;
        READFILEATTRIB(Float64, file_m, "dPhiGlobal", &dphi);
        OpalData::getInstance()->setGlobalPhaseShift(dphi);
    }

    if (attributeNames.find("nAutoPhaseCavities") != attributeNames.end()) {
        auto opal = OpalData::getInstance();
        h5_float64_t phase;
        h5_int64_t numAutoPhaseCavities;
        if (!H5HasFileAttrib(file_m, "nAutoPhaseCavities") ||
            H5ReadFileAttribInt64(file_m, "nAutoPhaseCavities", &numAutoPhaseCavities) != H5_SUCCESS) {
            numAutoPhaseCavities = 0;
        } else {
            for(long i = 0; i < numAutoPhaseCavities; ++ i) {
                std::string elementName  = "Cav-" + std::to_string(i + 1) + "-name";
                std::string elementPhase = "Cav-" + std::to_string(i + 1) + "-value";

                char name[128];
                READFILEATTRIB(String, file_m, elementName.c_str(), name);
                READFILEATTRIB(Float64, file_m, elementPhase.c_str(), &phase);

                opal->setMaxPhase(name, phase);
            }
        }
    }
}

void H5PartWrapperForPS::readStep(PartBunchBase<double, 3>* bunch, h5_ssize_t firstParticle, h5_ssize_t lastParticle) {
    h5_ssize_t numStepsInSource = H5GetNumSteps(file_m);
    h5_ssize_t readStep = numStepsInSource - 1;
    REPORTONERROR(H5SetStep(file_m, readStep));

    readStepHeader(bunch);
    readStepData(bunch, firstParticle, lastParticle);
}

void H5PartWrapperForPS::readStepHeader(PartBunchBase<double, 3>* bunch) {
    double actualT;
    READSTEPATTRIB(Float64, file_m, "TIME", &actualT);
    bunch->setT(actualT);
}

void H5PartWrapperForPS::readStepData(PartBunchBase<double, 3>* bunch, h5_ssize_t firstParticle, h5_ssize_t lastParticle) {
    h5_ssize_t numParticles = getNumParticles();
    if (lastParticle >= numParticles || firstParticle > lastParticle) {
        throw OpalException("H5PartWrapperForPS::readStepData",
                            "the provided particle numbers don't match the number of particles in the file");
    }

    EnvelopeBunch *ebunch = static_cast<EnvelopeBunch*>(bunch);

    REPORTONERROR(H5PartSetView(file_m, firstParticle, lastParticle));

    numParticles = lastParticle - firstParticle + 1;

    std::vector<char> buffer(numParticles * sizeof(h5_float64_t));
    h5_float64_t *f64buffer = reinterpret_cast<h5_float64_t*>(&buffer[0]);

    READDATA(Float64, file_m, "x", f64buffer);
    for(long int n = 0; n < numParticles; ++ n) {
        ebunch->setX(n, f64buffer[n]);
    }

    READDATA(Float64, file_m, "y", f64buffer);
    for(long int n = 0; n < numParticles; ++ n) {
        ebunch->setY(n, f64buffer[n]);
    }

    READDATA(Float64, file_m, "z", f64buffer);
    for(long int n = 0; n < numParticles; ++ n) {
        ebunch->setZ(n, f64buffer[n]);
    }

    READDATA(Float64, file_m, "px", f64buffer);
    for(long int n = 0; n < numParticles; ++ n) {
        ebunch->setPx(n, f64buffer[n]);
    }

    READDATA(Float64, file_m, "py", f64buffer);
    for(long int n = 0; n < numParticles; ++ n) {
        ebunch->setPy(n, f64buffer[n]);
    }

    READDATA(Float64, file_m, "beta", f64buffer);
    for(long int n = 0; n < numParticles; ++ n) {
        ebunch->setBeta(n, f64buffer[n]);
    }

    READDATA(Float64, file_m, "X0", f64buffer);
    for(long int n = 0; n < numParticles; ++ n) {
        ebunch->setX0(n, f64buffer[n]);
    }

    READDATA(Float64, file_m, "Y0", f64buffer);
    for(long int n = 0; n < numParticles; ++ n) {
        ebunch->setY0(n, f64buffer[n]);
    }

    READDATA(Float64, file_m, "pX0", f64buffer);
    for(long int n = 0; n < numParticles; ++ n) {
        ebunch->setPx0(n, f64buffer[n]);
    }

    READDATA(Float64, file_m, "pY0", f64buffer);
    for(long int n = 0; n < numParticles; ++ n) {
        ebunch->setPy0(n, f64buffer[n]);
    }

    REPORTONERROR(H5PartSetView(file_m, -1, -1));
}

void H5PartWrapperForPS::writeHeader() {
    h5_int64_t dumpfreq = Options::psDumpFreq;
    h5_float64_t dphi = OpalData::getInstance()->getGlobalPhaseShift();
    std::stringstream OPAL_version;

    OPAL_version << OPAL_PROJECT_NAME << " " << OPAL_PROJECT_VERSION << " # git rev. " << Util::getGitRevision();
    WRITESTRINGFILEATTRIB(file_m, "OPAL_version", OPAL_version.str().c_str());

    WRITESTRINGFILEATTRIB(file_m, "xUnit", "m");
    WRITESTRINGFILEATTRIB(file_m, "yUnit", "m");
    WRITESTRINGFILEATTRIB(file_m, "zUnit", "m");
    WRITESTRINGFILEATTRIB(file_m, "pxUnit", "#beta#gamma");
    WRITESTRINGFILEATTRIB(file_m, "pyUnit", "#beta#gamma");
    WRITESTRINGFILEATTRIB(file_m, "pzUnit", "#beta#gamma");
    WRITESTRINGFILEATTRIB(file_m, "qUnit", "Cb");

    WRITESTRINGFILEATTRIB(file_m, "idUnit", "1");

    WRITESTRINGFILEATTRIB(file_m, "SPOSUnit", "m");
    WRITESTRINGFILEATTRIB(file_m, "TIMEUnit", "s");
    WRITESTRINGFILEATTRIB(file_m, "ENERGYUnit", "MeV");
    WRITESTRINGFILEATTRIB(file_m, "#varepsilonUnit", "m rad");
    WRITESTRINGFILEATTRIB(file_m, "#varepsilon-geomUnit", "m rad");

    WRITESTRINGFILEATTRIB(file_m, "#sigmaUnit", "1");
    WRITESTRINGFILEATTRIB(file_m, "RMSXUnit", "m");
    WRITESTRINGFILEATTRIB(file_m, "RMSRUnit", "m");
    WRITESTRINGFILEATTRIB(file_m, "RMSPUnit", "#beta#gamma");

    WRITESTRINGFILEATTRIB(file_m, "MASSUnit", "GeV");
    WRITESTRINGFILEATTRIB(file_m, "CHARGEUnit", "C");

    WRITESTRINGFILEATTRIB(file_m, "spos-headUnit", "m");
    WRITESTRINGFILEATTRIB(file_m, "E-headUnit", "MV/m");
    WRITESTRINGFILEATTRIB(file_m, "B-headUnit", "T");
    WRITESTRINGFILEATTRIB(file_m, "spos-refUnit", "m");
    WRITESTRINGFILEATTRIB(file_m, "E-refUnit", "MV/m");
    WRITESTRINGFILEATTRIB(file_m, "B-refUnit", "T");
    WRITESTRINGFILEATTRIB(file_m, "spos-tailUnit", "m");
    WRITESTRINGFILEATTRIB(file_m, "E-tailUnit", "MV/m");
    WRITESTRINGFILEATTRIB(file_m, "B-tailUnit", "T");

    WRITESTRINGFILEATTRIB(file_m, "NumBunchUnit", "1");
    WRITESTRINGFILEATTRIB(file_m, "NumPartUnit", "1");

    WRITESTRINGFILEATTRIB(file_m, "RefPartRUnit", "m");
    WRITESTRINGFILEATTRIB(file_m, "RefPartPUnit", "#beta#gamma");

    /// Write file dump frequency.
    WRITEFILEATTRIB(Int64, file_m, "dump frequency", &dumpfreq, 1);

    /// Write global phase change
    WRITEFILEATTRIB(Float64, file_m, "dPhiGlobal", &dphi, 1);
}

void H5PartWrapperForPS::writeStep(PartBunchBase<double, 3>* bunch, const std::map<std::string, double> &additionalStepAttributes) {
    if (static_cast<EnvelopeBunch*>(bunch)->getTotalNum() == 0) return;

    open(H5_O_APPENDONLY);
    writeStepHeader(bunch, additionalStepAttributes);
    writeStepData(bunch);
    close();
}

void H5PartWrapperForPS::writeStepHeader(PartBunchBase<double, 3>* bunch,
                                         const std::map<std::string, double> &additionalStepAttributes) {
    EnvelopeBunch *ebunch = static_cast<EnvelopeBunch*>(bunch);
    ebunch->calcBeamParameters();

    double   actPos   = ebunch->get_sPos();
    double   t        = ebunch->getT();
    Vector_t rmin     = ebunch->get_origin();
    Vector_t rmax     = ebunch->get_maxExtent();

    Vector_t xsigma = ebunch->get_rrms();
    Vector_t psigma = ebunch->get_prms();
    Vector_t vareps = ebunch->get_norm_emit();

    double meanEnergy = ebunch->get_meanKineticEnergy();

    double mass = 1.0e-9 * ebunch->getM();
    double charge = ebunch->getCharge();

    Vector_t minP = ebunch->minP();
    Vector_t maxP = ebunch->maxP();

    // TODO:
    Vector_t RefPartR(0.0);
    Vector_t RefPartP(0.0);
    Vector_t centroid(0.0);
    Vector_t geomvareps(0.0);
    double I_0 = 4.0 * Physics::pi * Physics::epsilon_0 * Physics::c * ebunch->getM() / ebunch->getQ();
    double sigma = ((xsigma[0] * xsigma[0]) + (xsigma[1] * xsigma[1])) /
        (2.0 * ebunch->get_gamma() * I_0 * (geomvareps[0] * geomvareps[0] + geomvareps[1] * geomvareps[1]));

    /* ------------------------------------------------------------------------ */

    REPORTONERROR(H5SetStep(file_m, numSteps_m));

    char const *OPALFlavour = "opal-env";
    WRITESTRINGSTEPATTRIB(file_m, "OPAL_flavour", OPALFlavour);
    WRITESTEPATTRIB(Float64, file_m, "TIME", &t, 1);
    WRITESTEPATTRIB(Float64, file_m, "SPOS", &actPos, 1);
    WRITESTEPATTRIB(Float64, file_m, "RMSX", (h5_float64_t *)&xsigma, 3);
    WRITESTEPATTRIB(Float64, file_m, "RMSP", (h5_float64_t *)&psigma, 3);
    WRITESTEPATTRIB(Float64, file_m, "#varepsilon", (h5_float64_t *)&vareps, 3);

    WRITESTEPATTRIB(Float64, file_m, "ENERGY", &meanEnergy, 1);

    WRITESTEPATTRIB(Float64, file_m, "RefPartR", (h5_float64_t *)&RefPartR, 3);
    WRITESTEPATTRIB(Float64, file_m, "RefPartP", (h5_float64_t *)&RefPartP, 3);
    WRITESTEPATTRIB(Float64, file_m, "centroid", (h5_float64_t *)&centroid, 3);
    WRITESTEPATTRIB(Float64, file_m, "#varepsilon-geom", (h5_float64_t *)&geomvareps, 3);

    WRITESTEPATTRIB(Float64, file_m, "minX", (h5_float64_t *)&rmin, 3);
    WRITESTEPATTRIB(Float64, file_m, "maxX", (h5_float64_t *)&rmax, 3);

    WRITESTEPATTRIB(Float64, file_m, "minP", (h5_float64_t *)&minP, 3);
    WRITESTEPATTRIB(Float64, file_m, "maxP", (h5_float64_t *)&maxP, 3);

    WRITESTEPATTRIB(Float64, file_m, "#sigma", &sigma, 1);

    /// Write particle mass and charge per particle. (Consider making these file attributes.)
    WRITESTEPATTRIB(Float64, file_m, "MASS", &mass, 1);
    WRITESTEPATTRIB(Float64, file_m, "CHARGE", &charge, 1);

    try {
        double sposHead = additionalStepAttributes.at("sposHead");
        double sposRef = additionalStepAttributes.at("sposRef");
        double sposTail = additionalStepAttributes.at("sposTail");
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

        WRITESTEPATTRIB(Float64, file_m, "spos-head", &sposHead, 1);
        WRITESTEPATTRIB(Float64, file_m, "spos-ref", &sposRef, 1);
        WRITESTEPATTRIB(Float64, file_m, "spos-tail", &sposTail, 1);

        WRITESTEPATTRIB(Float64, file_m, "B-ref", (h5_float64_t *)&referenceB, 3);
        WRITESTEPATTRIB(Float64, file_m, "E-ref", (h5_float64_t *)&referenceE, 3);
        WRITESTEPATTRIB(Float64, file_m, "B-head", (h5_float64_t *)&headB, 3);
        WRITESTEPATTRIB(Float64, file_m, "E-head", (h5_float64_t *)&headE, 3);
        WRITESTEPATTRIB(Float64, file_m, "B-tail", (h5_float64_t *)&tailB, 3);
        WRITESTEPATTRIB(Float64, file_m, "E-tail", (h5_float64_t *)&tailE, 3);
    } catch (std::out_of_range & m) {
        ERRORMSG(m.what() << endl);

        throw OpalException("H5PartWrapperForPC::writeStepHeader",
                            "some additional step attribute not found");
    }
    ++ numSteps_m;
}

void H5PartWrapperForPS::writeStepData(PartBunchBase<double, 3>* bunch) {
    EnvelopeBunch *ebunch = static_cast<EnvelopeBunch*>(bunch);

    size_t numLocalParticles = ebunch->getLocalNum();

    std::vector<char> buffer(numLocalParticles * sizeof(h5_float64_t));
    h5_float64_t *f64buffer = reinterpret_cast<h5_float64_t*>(&buffer[0]);

    REPORTONERROR(H5PartSetNumParticles(file_m, numLocalParticles));

    for(size_t i = 0; i < numLocalParticles; ++ i)
        f64buffer[i] =  ebunch->getX(i);
    WRITEDATA(Float64, file_m, "x", f64buffer);

    for(size_t i = 0; i < numLocalParticles; ++ i)
        f64buffer[i] =  ebunch->getY(i);
    WRITEDATA(Float64, file_m, "y", f64buffer);

    for(size_t i = 0; i < numLocalParticles; ++ i)
        f64buffer[i] =  ebunch->getZ(i);
    WRITEDATA(Float64, file_m, "z", f64buffer);

    for(size_t i = 0; i < numLocalParticles; ++ i)
        f64buffer[i] =  ebunch->getPx(i);
    WRITEDATA(Float64, file_m, "px", f64buffer);

    for(size_t i = 0; i < numLocalParticles; ++ i)
        f64buffer[i] =  ebunch->getPy(i);
    WRITEDATA(Float64, file_m, "py", f64buffer);

    for(size_t i = 0; i < numLocalParticles; ++ i)
        f64buffer[i] =  ebunch->getPz(i);
    WRITEDATA(Float64, file_m, "pz", f64buffer);

    for(size_t i = 0; i < numLocalParticles; ++ i)
        f64buffer[i] =  ebunch->getBeta(i);
    WRITEDATA(Float64, file_m, "beta", f64buffer);

    for(size_t i = 0; i < numLocalParticles; ++ i)
        f64buffer[i] =  ebunch->getX0(i);
    WRITEDATA(Float64, file_m, "X0", f64buffer);

    for(size_t i = 0; i < numLocalParticles; ++ i)
        f64buffer[i] =  ebunch->getY0(i);
    WRITEDATA(Float64, file_m, "Y0", f64buffer);

    for(size_t i = 0; i < numLocalParticles; ++ i)
        f64buffer[i] =  ebunch->getPx0(i);
    WRITEDATA(Float64, file_m, "pX0", f64buffer);

    for(size_t i = 0; i < numLocalParticles; ++ i)
        f64buffer[i] =  ebunch->getPy0(i);
    WRITEDATA(Float64, file_m, "pY0", f64buffer);
}

void H5PartWrapperForPS::stashPhaseSpaceEnvelope(EnvelopeBunch &bunch,
                                                 Vector_t FDext[],
                                                 double sposHead,
                                                 double sposRef,
                                                 double sposTail) {
    EnvelopeBunch *ebunch = static_cast<EnvelopeBunch*>(&bunch);

    /// Calculate bunch statistical parameters etc. Put them in the right format
    /// for H5 file write.
    ebunch->calcBeamParameters();

    size_t numLocalParticles = ebunch->getLocalNum();
    //TODO:
    stash_RefPartR.push_back(Vector_t(0.0)); //ebunch->RefPart_R;
    stash_RefPartP.push_back(Vector_t(0.0)); //ebunch->RefPart_P;
    stash_centroid.push_back(Vector_t(0.0));

    stash_actPos.push_back(ebunch->get_sPos());
    stash_t.push_back(ebunch->getT());

    stash_nTot.push_back(ebunch->getTotalNum());
    stash_nLoc.push_back(numLocalParticles);

    stash_xsigma.push_back(ebunch->sigmax());
    stash_psigma.push_back(ebunch->sigmap());
    stash_vareps.push_back(ebunch->get_norm_emit());
    stash_geomvareps.push_back(ebunch->emtn());
    stash_rmin.push_back(ebunch->minX());
    stash_rmax.push_back(ebunch->maxX());
    stash_maxP.push_back(ebunch->maxP());
    stash_minP.push_back(ebunch->minP());

    //in MeV
    stash_meanEnergy.push_back(ebunch->get_meanKineticEnergy() * 1e-6);

    stash_sposHead.push_back(sposHead);
    stash_sposRef.push_back(sposRef);
    stash_sposTail.push_back(sposTail);
    stash_Bhead.push_back(FDext[0]);
    stash_Ehead.push_back(FDext[1]);
    stash_Bref.push_back(FDext[2]);
    stash_Eref.push_back(FDext[3]);
    stash_Btail.push_back(FDext[4]);
    stash_Etail.push_back(FDext[5]);

    /// Write particle mass and charge per particle. (Consider making these
    /// file attributes.)
    stash_mass.push_back(1.0e-9 * ebunch->getM());
    stash_charge.push_back(ebunch->getChargePerParticle());

    /// Write bunch phase space.
    //for (size_t i=0; i<numLocalParticles;i++)
    //farray[i] =  bunch->getX(i);

    //for (size_t i=0; i<numLocalParticles;i++)
    //farray[i] =  bunch->getY(i);

    //for (size_t i=0; i<numLocalParticles;i++)
    //farray[i] =  bunch->getZ(i);

    //for (size_t i=0; i<numLocalParticles;i++)
    //farray[i] =  bunch->getPx(i);

    //for (size_t i=0; i<numLocalParticles;i++)
    //farray[i] =  bunch->getPy(i);

    //for (size_t i=0; i<numLocalParticles;i++)
    //farray[i] =  bunch->getPz(i);

     ++ numSteps_m;
}

void H5PartWrapperForPS::dumpStashedPhaseSpaceEnvelope() {

    size_t numLocalParticles = stash_nLoc[0];
    std::vector<char> buffer(numLocalParticles * sizeof(h5_float64_t));
    h5_float64_t *f64buffer = reinterpret_cast<h5_float64_t*>(&buffer[0]);

    //FIXME: restart step
    for(int step = 0; step < numSteps_m; ++ step) {
        numLocalParticles = stash_nLoc[step];

        /// Set current record/time step.
        REPORTONERROR(H5SetStep(file_m, step));
        REPORTONERROR(H5PartSetNumParticles(file_m, numLocalParticles));

        /// Write statistical data.
        WRITESTEPATTRIB(Float64, file_m, "SPOS", &stash_actPos[step], 1);
        WRITESTEPATTRIB(Float64, file_m, "RMSX", &stash_xsigma[step](0), 3);
        WRITESTEPATTRIB(Float64, file_m, "RMSP", &stash_psigma[step](0), 3);
        WRITESTEPATTRIB(Float64, file_m, "maxX", &stash_rmax[step](0), 3);
        WRITESTEPATTRIB(Float64, file_m, "minX", &stash_rmin[step](0), 3);
        WRITESTEPATTRIB(Float64, file_m, "maxP", &stash_maxP[step](0), 3);
        WRITESTEPATTRIB(Float64, file_m, "minP", &stash_minP[step](0), 3);
        WRITESTEPATTRIB(Float64, file_m, "centroid", &stash_centroid[step](0), 3);
        WRITESTEPATTRIB(Float64, file_m, "TIME", &stash_t[step], 1);
        WRITESTEPATTRIB(Float64, file_m, "ENERGY", &stash_meanEnergy[step], 1);
        WRITESTEPATTRIB(Float64, file_m, "RefPartR", &stash_RefPartR[step](0), 3);
        WRITESTEPATTRIB(Float64, file_m, "RefPartP", &stash_RefPartP[step](0), 3);
        WRITESTEPATTRIB(Float64, file_m, "B-head", &stash_Bhead[step](0), 3);
        WRITESTEPATTRIB(Float64, file_m, "E-head", &stash_Ehead[step](0), 3);
        WRITESTEPATTRIB(Float64, file_m, "B-ref", &stash_Bref[step](0), 3);
        WRITESTEPATTRIB(Float64, file_m, "E-ref", &stash_Eref[step](0), 3);
        WRITESTEPATTRIB(Float64, file_m, "B-tail", &stash_Btail[step](0), 3);
        WRITESTEPATTRIB(Float64, file_m, "E-tail", &stash_Etail[step](0), 3);
        WRITESTEPATTRIB(Float64, file_m, "spos-head", &stash_sposHead[step], 1);
        WRITESTEPATTRIB(Float64, file_m, "spos-ref", &stash_sposRef[step], 1);
        WRITESTEPATTRIB(Float64, file_m, "spos-tail", &stash_sposTail[step], 1);

        /// Write particle mass and charge per particle. (Consider making these
        /// file attributes.)
        WRITESTEPATTRIB(Float64, file_m, "MASS", &stash_mass[step], 1);
        WRITESTEPATTRIB(Float64, file_m, "CHARGE", &stash_charge[step], 1);
        /// Write normalized emittance.
        WRITESTEPATTRIB(Float64, file_m, "#varepsilon", &stash_vareps[step](0), 3);
        /// Write geometric emittance.
        WRITESTEPATTRIB(Float64, file_m, "#varepsilon-geom", &stash_geomvareps[step](0), 3);

        /// Write bunch phase space.
        for(size_t i = 0; i < numLocalParticles; i++)
            f64buffer[i] =  0.0;
        WRITEDATA(Float64, file_m, "x", f64buffer);
        for(size_t i = 0; i < numLocalParticles; i++)
            f64buffer[i] =  0.0; //bunch->getY(i);
        WRITEDATA(Float64, file_m, "y", f64buffer);
        for(size_t i = 0; i < numLocalParticles; i++)
            f64buffer[i] =  0.0; //bunch->getZ(i);
        WRITEDATA(Float64, file_m, "z", f64buffer);
        for(size_t i = 0; i < numLocalParticles; i++)
            f64buffer[i] =  0.0; //bunch->getPx(i);
        WRITEDATA(Float64, file_m, "px", f64buffer);
        for(size_t i = 0; i < numLocalParticles; i++)
            f64buffer[i] =  0.0; //bunch->getPy(i);
        WRITEDATA(Float64, file_m, "py", f64buffer);
        for(size_t i = 0; i < numLocalParticles; i++)
            f64buffer[i] =  0.0; //bunch->getPz(i);
        WRITEDATA(Float64, file_m, "pz", f64buffer);
        //rc = H5Fflush, file_m->file, rc = H5F_SCOPE_GLOBAL);
    }
}
