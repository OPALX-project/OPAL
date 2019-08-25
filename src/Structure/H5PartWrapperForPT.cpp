//
//  Copyright & License: See Copyright.readme in src directory
//

#include "Structure/H5PartWrapperForPT.h"

#include "OPALconfig.h"
#include "AbstractObjects/OpalData.h"
#include "Utilities/Options.h"
#include "Utilities/Util.h"
#include "Physics/Physics.h"

#include "h5core/h5_types.h"
#include <boost/filesystem.hpp>

#include <sstream>
#include <set>

H5PartWrapperForPT::H5PartWrapperForPT(const std::string &fileName, h5_int32_t flags):
    H5PartWrapper(fileName, flags)
{ }

H5PartWrapperForPT::H5PartWrapperForPT(const std::string &fileName, int restartStep, std::string sourceFile, h5_int32_t flags):
    H5PartWrapper(fileName, restartStep, sourceFile, flags)
{
    if (restartStep == -1) {
        restartStep = H5GetNumSteps(file_m) - 1 ;
        OpalData::getInstance()->setRestartStep(restartStep);
    }
}

H5PartWrapperForPT::~H5PartWrapperForPT()
{ }

void H5PartWrapperForPT::readHeader() {
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
        h5_int64_t numAutoPhaseCavities = 0;
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

void H5PartWrapperForPT::readStep(PartBunchBase<double, 3>* bunch, h5_ssize_t firstParticle, h5_ssize_t lastParticle) {
    h5_ssize_t numStepsInSource = H5GetNumSteps(file_m);
    h5_ssize_t readStep = numStepsInSource - 1;
    REPORTONERROR(H5SetStep(file_m, readStep));

    readStepHeader(bunch);
    readStepData(bunch, firstParticle, lastParticle);
}

void H5PartWrapperForPT::readStepHeader(PartBunchBase<double, 3>* bunch) {
    double actualT;
    READSTEPATTRIB(Float64, file_m, "TIME", &actualT);
    bunch->setT(actualT);

    double spos;
    READSTEPATTRIB(Float64, file_m, "SPOS", &spos);
    bunch->set_sPos(spos);

    h5_int64_t ltstep;
    READSTEPATTRIB(Int64, file_m, "LocalTrackStep", &ltstep);
    bunch->setLocalTrackStep((long long)(ltstep + 1));

    h5_int64_t gtstep;
    READSTEPATTRIB(Int64, file_m, "GlobalTrackStep", &gtstep);
    bunch->setGlobalTrackStep((long long)(gtstep + 1));

    Vector_t RefPartR;
    READSTEPATTRIB(Float64, file_m, "RefPartR", (h5_float64_t *)&RefPartR);
    bunch->RefPartR_m = RefPartR;

    Vector_t RefPartP;
    READSTEPATTRIB(Float64, file_m, "RefPartP", (h5_float64_t *)&RefPartP);
    bunch->RefPartP_m = RefPartP;

    Vector_t TaitBryant;
    READSTEPATTRIB(Float64, file_m, "TaitBryantAngles", (h5_float64_t *)&TaitBryant);
    Quaternion rotTheta(cos(0.5 * TaitBryant[0]), 0, sin(0.5 * TaitBryant[0]), 0);
    Quaternion rotPhi(cos(0.5 * TaitBryant[1]), sin(0.5 * TaitBryant[1]), 0, 0);
    Quaternion rotPsi(cos(0.5 * TaitBryant[2]), 0, 0, sin(0.5 * TaitBryant[2]));
    Quaternion rotation = rotTheta * (rotPhi * rotPsi);
    bunch->toLabTrafo_m = CoordinateSystemTrafo(-rotation.conjugate().rotate(RefPartR), rotation);
}

void H5PartWrapperForPT::readStepData(PartBunchBase<double, 3>* bunch, h5_ssize_t firstParticle, h5_ssize_t lastParticle) {
    h5_ssize_t numParticles = getNumParticles();
    if (lastParticle >= numParticles || firstParticle > lastParticle) {
        throw OpalException("H5PartWrapperForPT::readStepData",
                            "the provided particle numbers don't match the number of particles in the file");
    }

    REPORTONERROR(H5PartSetView(file_m, firstParticle, lastParticle));

    numParticles = lastParticle - firstParticle + 1;

    std::vector<char> buffer(numParticles * sizeof(h5_float64_t));
    h5_float64_t *f64buffer = reinterpret_cast<h5_float64_t*>(&buffer[0]);
    h5_int32_t *i32buffer = reinterpret_cast<h5_int32_t*>(&buffer[0]);

    READDATA(Float64, file_m, "x", f64buffer);
    for(long int n = 0; n < numParticles; ++ n) {
        bunch->R[n](0) = f64buffer[n];
        bunch->Bin[n] = 0;
    }

    READDATA(Float64, file_m, "y", f64buffer);
    for(long int n = 0; n < numParticles; ++ n) {
        bunch->R[n](1) = f64buffer[n];
    }

    READDATA(Float64, file_m, "z", f64buffer);
    for(long int n = 0; n < numParticles; ++ n) {
        bunch->R[n](2) = f64buffer[n];
    }

    READDATA(Float64, file_m, "px", f64buffer);
    for(long int n = 0; n < numParticles; ++ n) {
        bunch->P[n](0) = f64buffer[n];
    }

    READDATA(Float64, file_m, "py", f64buffer);
    for(long int n = 0; n < numParticles; ++ n) {
        bunch->P[n](1) = f64buffer[n];
    }

    READDATA(Float64, file_m, "pz", f64buffer);
    for(long int n = 0; n < numParticles; ++ n) {
        bunch->P[n](2) = f64buffer[n];
    }

    READDATA(Float64, file_m, "q", f64buffer);
    for(long int n = 0; n < numParticles; ++ n) {
        bunch->Q[n] = f64buffer[n];
    }

    READDATA(Int32, file_m, "id", i32buffer);
    for(long int n = 0; n < numParticles; ++ n) {
        bunch->ID[n] = i32buffer[n];
    }

    REPORTONERROR(H5PartSetView(file_m, -1, -1));
}

void H5PartWrapperForPT::writeHeader() {
    std::stringstream OPAL_version;
    OPAL_version << OPAL_PROJECT_NAME << " " << OPAL_PROJECT_VERSION << " # git rev. " << Util::getGitRevision();
    WRITESTRINGFILEATTRIB(file_m, "OPAL_version", OPAL_version.str().c_str());


    WRITESTRINGFILEATTRIB(file_m, "idUnit", "1");
    WRITESTRINGFILEATTRIB(file_m, "xUnit", "m");
    WRITESTRINGFILEATTRIB(file_m, "yUnit", "m");
    WRITESTRINGFILEATTRIB(file_m, "zUnit", "m");

    WRITESTRINGFILEATTRIB(file_m, "pxUnit", "#beta#gamma");
    WRITESTRINGFILEATTRIB(file_m, "pyUnit", "#beta#gamma");
    WRITESTRINGFILEATTRIB(file_m, "pzUnit", "#beta#gamma");

    WRITESTRINGFILEATTRIB(file_m, "ptypeUnit", "1");
    WRITESTRINGFILEATTRIB(file_m, "qUnit", "C");

    if (Options::ebDump) {
        WRITESTRINGFILEATTRIB(file_m, "ExUnit", "MV/m");
        WRITESTRINGFILEATTRIB(file_m, "EyUnit", "MV/m");
        WRITESTRINGFILEATTRIB(file_m, "EzUnit", "MV/m");

        WRITESTRINGFILEATTRIB(file_m, "BxUnit", "T");
        WRITESTRINGFILEATTRIB(file_m, "ByUnit", "T");
        WRITESTRINGFILEATTRIB(file_m, "BzUnit", "T");
    }

    if (Options::rhoDump) {
        WRITESTRINGFILEATTRIB(file_m, "rhoUnit", "C/m^3");
    }

    WRITESTRINGFILEATTRIB(file_m, "TaitBryantAnglesUnit", "rad");

    WRITESTRINGFILEATTRIB(file_m, "SPOSUnit", "m");
    WRITESTRINGFILEATTRIB(file_m, "TIMEUnit", "s");
    WRITESTRINGFILEATTRIB(file_m, "#gammaUnit", "1");
    WRITESTRINGFILEATTRIB(file_m, "ENERGYUnit", "MeV");
    WRITESTRINGFILEATTRIB(file_m, "#varepsilonUnit", "m rad");
    WRITESTRINGFILEATTRIB(file_m, "#varepsilon-geomUnit", "m rad");

    WRITESTRINGFILEATTRIB(file_m, "#sigmaUnit", "1");
    WRITESTRINGFILEATTRIB(file_m, "RMSXUnit", "m");
    WRITESTRINGFILEATTRIB(file_m, "centroidUnit", "m");
    WRITESTRINGFILEATTRIB(file_m, "RMSPUnit", "#beta#gamma");
    WRITESTRINGFILEATTRIB(file_m, "MEANPUnit", "#beta#gamma");

    WRITESTRINGFILEATTRIB(file_m, "minXUnit", "m");
    WRITESTRINGFILEATTRIB(file_m, "maxXUnit", "m");
    WRITESTRINGFILEATTRIB(file_m, "minPUnit", "#beta#gamma");
    WRITESTRINGFILEATTRIB(file_m, "maxPUnit", "#beta#gamma");

    WRITESTRINGFILEATTRIB(file_m, "dEUnit", "MeV");

    WRITESTRINGFILEATTRIB(file_m, "MASSUnit", "GeV");
    WRITESTRINGFILEATTRIB(file_m, "CHARGEUnit", "C");

    WRITESTRINGFILEATTRIB(file_m, "E-refUnit", "MV/m");
    WRITESTRINGFILEATTRIB(file_m, "B-refUnit", "T");

    WRITESTRINGFILEATTRIB(file_m, "StepUnit", "1");
    WRITESTRINGFILEATTRIB(file_m, "LocalTrackStepUnit", "1");
    WRITESTRINGFILEATTRIB(file_m, "GlobalTrackStepUnit", "1");
    WRITESTRINGFILEATTRIB(file_m, "NumBunchUnit", "1");
    WRITESTRINGFILEATTRIB(file_m, "NumPartUnit", "1");

    WRITESTRINGFILEATTRIB(file_m, "RefPartRUnit", "m");
    WRITESTRINGFILEATTRIB(file_m, "RefPartPUnit", "#beta#gamma");
    WRITESTRINGFILEATTRIB(file_m, "SteptoLastInjUnit", "1");

    WRITESTRINGFILEATTRIB(file_m, "dump frequencyUnit", "1")
    WRITESTRINGFILEATTRIB(file_m, "dPhiGlobalUnit", "rad")

    /// Write file dump frequency.
    h5_int64_t dumpfreq = Options::psDumpFreq;
    WRITEFILEATTRIB(Int64, file_m, "dump frequency", &dumpfreq, 1);


    /// Write global phase change
    h5_float64_t dphi = OpalData::getInstance()->getGlobalPhaseShift();
    WRITEFILEATTRIB(Float64, file_m, "dPhiGlobal", &dphi, 1);
}

void H5PartWrapperForPT::writeStep(PartBunchBase<double, 3>* bunch, const std::map<std::string, double> &additionalStepAttributes) {
    if (bunch->getTotalNum() == 0) return;

    open(H5_O_APPENDONLY);
    writeStepHeader(bunch, additionalStepAttributes);
    writeStepData(bunch);
    close();
}

void H5PartWrapperForPT::writeStepHeader(PartBunchBase<double, 3>* bunch, const std::map<std::string, double> &additionalStepAttributes) {
    bunch->calcBeamParameters();

    double   actPos   = bunch->get_sPos();
    double   t        = bunch->getT();
    Vector_t rmin     = bunch->get_origin();
    Vector_t rmax     = bunch->get_maxExtent();
    Vector_t centroid = bunch->get_centroid();

    Vector_t maxP(0.0);
    Vector_t minP(0.0);

    Vector_t xsigma = bunch->get_rrms();
    Vector_t psigma = bunch->get_prms();
    Vector_t vareps = bunch->get_norm_emit();
    Vector_t geomvareps = bunch->get_emit();
    Vector_t RefPartR = bunch->RefPartR_m;
    Vector_t RefPartP = bunch->RefPartP_m;
    Vector_t TaitBryant = Util::getTaitBryantAngles(bunch->toLabTrafo_m.getRotation());
    Vector_t pmean = bunch->get_pmean();

    double meanEnergy = bunch->get_meanKineticEnergy();
    double energySpread = bunch->getdE();
    double I_0 = 4.0 * Physics::pi * Physics::epsilon_0 * Physics::c * bunch->getM() / bunch->getQ();
    double sigma = ((xsigma[0] * xsigma[0]) + (xsigma[1] * xsigma[1])) /
        (2.0 * bunch->get_gamma() * I_0 * (geomvareps[0] * geomvareps[0] + geomvareps[1] * geomvareps[1]));

    h5_int64_t localTrackStep = (h5_int64_t)bunch->getLocalTrackStep();
    h5_int64_t globalTrackStep = (h5_int64_t)bunch->getGlobalTrackStep();

    double mass = 1.0e-9 * bunch->getM();
    double charge = bunch->getCharge();

    h5_int64_t numBunch = 1;
    h5_int64_t SteptoLastInj = 0;

    bunch->get_PBounds(minP, maxP);

    /* ------------------------------------------------------------------------ */

    REPORTONERROR(H5SetStep(file_m, numSteps_m));

    char const *OPALFlavour = "opal-t";
    WRITESTRINGSTEPATTRIB(file_m, "OPAL_flavour", OPALFlavour);
    WRITESTEPATTRIB(Float64, file_m, "SPOS", &actPos, 1);
    WRITESTEPATTRIB(Float64, file_m, "RefPartR", (h5_float64_t *)&RefPartR, 3);
    WRITESTEPATTRIB(Float64, file_m, "centroid", (h5_float64_t *)&centroid, 3);
    WRITESTEPATTRIB(Float64, file_m, "RMSX", (h5_float64_t *)&xsigma, 3);

    WRITESTEPATTRIB(Float64, file_m, "RefPartP", (h5_float64_t *)&RefPartP, 3);
    WRITESTEPATTRIB(Float64, file_m, "MEANP", (h5_float64_t *)&pmean, 3);
    WRITESTEPATTRIB(Float64, file_m, "RMSP", (h5_float64_t *)&psigma, 3);
    WRITESTEPATTRIB(Float64, file_m, "TaitBryantAngles", (h5_float64_t *)&TaitBryant, 3);

    WRITESTEPATTRIB(Float64, file_m, "#varepsilon", (h5_float64_t *)&vareps, 3);
    WRITESTEPATTRIB(Float64, file_m, "#varepsilon-geom", (h5_float64_t *)&geomvareps, 3);

    WRITESTEPATTRIB(Float64, file_m, "minX", (h5_float64_t *)&rmin, 3);
    WRITESTEPATTRIB(Float64, file_m, "maxX", (h5_float64_t *)&rmax, 3);

    WRITESTEPATTRIB(Float64, file_m, "minP", (h5_float64_t *)&minP, 3);
    WRITESTEPATTRIB(Float64, file_m, "maxP", (h5_float64_t *)&maxP, 3);

    WRITESTEPATTRIB(Int64, file_m, "Step", &numSteps_m, 1);
    WRITESTEPATTRIB(Int64, file_m, "LocalTrackStep", &localTrackStep, 1);
    WRITESTEPATTRIB(Int64, file_m, "GlobalTrackStep", &globalTrackStep, 1);

    WRITESTEPATTRIB(Float64, file_m, "#sigma", &sigma, 1);

    WRITESTEPATTRIB(Float64, file_m, "TIME", &t, 1);

    WRITESTEPATTRIB(Float64, file_m, "ENERGY", &meanEnergy, 1);
    WRITESTEPATTRIB(Float64, file_m, "dE", &energySpread, 1);

    /// Write particle mass and charge per particle. (Consider making these file attributes.)
    WRITESTEPATTRIB(Float64, file_m, "MASS", &mass, 1);

    WRITESTEPATTRIB(Float64, file_m, "CHARGE", &charge, 1);

    WRITESTEPATTRIB(Int64, file_m, "NumBunch", &numBunch, 1);

    WRITESTEPATTRIB(Int64, file_m, "SteptoLastInj", &SteptoLastInj, 1);

    try {
        Vector_t referenceB(additionalStepAttributes.at("B-ref_x"),
                            additionalStepAttributes.at("B-ref_z"),
                            additionalStepAttributes.at("B-ref_y"));
        Vector_t referenceE(additionalStepAttributes.at("E-ref_x"),
                            additionalStepAttributes.at("E-ref_z"),
                            additionalStepAttributes.at("E-ref_y"));

        WRITESTEPATTRIB(Float64, file_m, "B-ref", (h5_float64_t *)&referenceB, 3);
        WRITESTEPATTRIB(Float64, file_m, "E-ref", (h5_float64_t *)&referenceE, 3);
    } catch (std::out_of_range & m) {
        ERRORMSG(m.what() << endl);

        throw OpalException("H5PartWrapperForPC::writeStepHeader",
                            "some additional step attribute not found");
    }

    ++ numSteps_m;
}

void H5PartWrapperForPT::writeStepData(PartBunchBase<double, 3>* bunch) {
    size_t numLocalParticles = bunch->getLocalNum();

    REPORTONERROR(H5PartSetNumParticles(file_m, numLocalParticles));

    std::vector<char> buffer(numLocalParticles * sizeof(h5_float64_t));
    h5_float64_t *f64buffer = reinterpret_cast<h5_float64_t*>(&buffer[0]);
    h5_int64_t *i64buffer = reinterpret_cast<h5_int64_t*>(&buffer[0]);
    h5_int32_t *i32buffer = reinterpret_cast<h5_int32_t*>(&buffer[0]);

    for(size_t i = 0; i < numLocalParticles; ++ i)
        f64buffer[i] =  bunch->R[i](0);
    WRITEDATA(Float64, file_m, "x", f64buffer);

    for(size_t i = 0; i < numLocalParticles; ++ i)
        f64buffer[i] =  bunch->R[i](1);
    WRITEDATA(Float64, file_m, "y", f64buffer);

    for(size_t i = 0; i < numLocalParticles; ++ i)
        f64buffer[i] =  bunch->R[i](2);
    WRITEDATA(Float64, file_m, "z", f64buffer);

    for(size_t i = 0; i < numLocalParticles; ++ i)
        f64buffer[i] =  bunch->P[i](0);
    WRITEDATA(Float64, file_m, "px", f64buffer);

    for(size_t i = 0; i < numLocalParticles; ++ i)
        f64buffer[i] =  bunch->P[i](1);
    WRITEDATA(Float64, file_m, "py", f64buffer);

    for(size_t i = 0; i < numLocalParticles; ++ i)
        f64buffer[i] =  bunch->P[i](2);
    WRITEDATA(Float64, file_m, "pz", f64buffer);

    for(size_t i = 0; i < numLocalParticles; ++ i)
        f64buffer[i] =  bunch->Q[i];
    WRITEDATA(Float64, file_m, "q", f64buffer);

    for(size_t i = 0; i < numLocalParticles; ++ i)
        i64buffer[i] =  bunch->ID[i];
    WRITEDATA(Int64, file_m, "id", i64buffer);

    for(size_t i = 0; i < numLocalParticles; ++ i)
        i32buffer[i] = (h5_int32_t) bunch->PType[i];
    WRITEDATA(Int32, file_m, "ptype", i32buffer);

    if(Options::ebDump) {
        for(size_t i = 0; i < numLocalParticles; ++ i)
            f64buffer[i] =  bunch->Ef[i](0);
        WRITEDATA(Float64, file_m, "Ex", f64buffer);

        for(size_t i = 0; i < numLocalParticles; ++ i)
            f64buffer[i] =  bunch->Ef[i](1);
        WRITEDATA(Float64, file_m, "Ey", f64buffer);

        for(size_t i = 0; i < numLocalParticles; ++ i)
            f64buffer[i] =  bunch->Ef[i](2);
        WRITEDATA(Float64, file_m, "Ez", f64buffer);

        for(size_t i = 0; i < numLocalParticles; ++ i)
            f64buffer[i] =  bunch->Bf[i](0);
        WRITEDATA(Float64, file_m, "Bx", f64buffer);

        for(size_t i = 0; i < numLocalParticles; ++ i)
            f64buffer[i] =  bunch->Bf[i](1);
        WRITEDATA(Float64, file_m, "By", f64buffer);

        for(size_t i = 0; i < numLocalParticles; ++ i)
            f64buffer[i] =  bunch->Bf[i](2);
        WRITEDATA(Float64, file_m, "Bz", f64buffer);

    }

    /// Write space charge field map if asked for.
    if(Options::rhoDump) {
        NDIndex<3> idx = bunch->getFieldLayout().getLocalNDIndex();
        NDIndex<3> elem;
        h5_err_t herr = H5Block3dSetView(
                                         file_m,
                                         idx[0].min(), idx[0].max(),
                                         idx[1].min(), idx[1].max(),
                                         idx[2].min(), idx[2].max());
        reportOnError(herr, __FILE__, __LINE__);

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
        herr = H5Block3dWriteScalarFieldFloat64(file_m, "rho", data.get());
        reportOnError(herr, __FILE__, __LINE__);

        /// Need this to align particles and fields when writing space charge map.
        herr = H5Block3dSetFieldOrigin(file_m, "rho",
                                       (h5_float64_t)bunch->get_origin()(0),
                                       (h5_float64_t)bunch->get_origin()(1),
                                       (h5_float64_t)bunch->get_origin()(2));
        reportOnError(herr, __FILE__, __LINE__);

        herr = H5Block3dSetFieldSpacing(file_m, "rho",
                                        (h5_float64_t)bunch->get_hr()(0),
                                        (h5_float64_t)bunch->get_hr()(1),
                                        (h5_float64_t)bunch->get_hr()(2));
        reportOnError(herr, __FILE__, __LINE__);

    }
}
