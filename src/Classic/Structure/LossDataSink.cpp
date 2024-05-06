//
// Class LossDataSink
//   This class writes file attributes to describe phase space of loss files
//
// Copyright (c) 200x - 2020, Paul Scherrer Institut, Villigen PSI, Switzerland
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
#include "Structure/LossDataSink.h"

#include "AbstractObjects/OpalData.h"
#include "Algorithms/DistributionMoments.h"
#include "Message/GlobalComm.h"
#include "OPALconfig.h"
#include "Utilities/GeneralClassicException.h"
#include "Utilities/Options.h"
#include "Utilities/Util.h"

#include "Utility/IpplInfo.h"

#include <cmath>
#include <filesystem>

extern Inform* gmsg;

#define WRITE_FILEATTRIB_STRING( attribute, value ) {             \
    h5_int64_t h5err = H5WriteFileAttribString (H5file_m, attribute, value); \
    if (h5err <= H5_ERR) {                                              \
        std::stringstream ss;                                               \
        ss << "failed to write string attribute " << attribute << " to file " << fileName_m; \
        throw GeneralClassicException(std::string(__func__), ss.str()); \
    }\
}
#define WRITE_STEPATTRIB_FLOAT64( attribute, value, size ) { \
        h5_int64_t h5err = H5WriteStepAttribFloat64 (H5file_m, attribute, value, size); \
    if (h5err <= H5_ERR) { \
        std::stringstream ss; \
        ss << "failed to write float64 attribute " << attribute << " to file " << fileName_m; \
        throw GeneralClassicException(std::string(__func__), ss.str()); \
    }\
}
#define WRITE_STEPATTRIB_INT64( attribute, value, size ) { \
        h5_int64_t h5err = H5WriteStepAttribInt64 (H5file_m, attribute, value, size); \
    if (h5err <= H5_ERR) { \
        std::stringstream ss; \
        ss << "failed to write int64 attribute " << attribute << " to file " << fileName_m; \
        throw GeneralClassicException(std::string(__func__), ss.str()); \
    }\
}
#define WRITE_DATA_FLOAT64( name, value ) {                             \
        h5_int64_t h5err = H5PartWriteDataFloat64 (H5file_m, name, value); \
    if (h5err <= H5_ERR) { \
        std::stringstream ss; \
        ss << "failed to write float64 data " << name << " to file " << fileName_m; \
        throw GeneralClassicException(std::string(__func__), ss.str()); \
    }\
}
#define WRITE_DATA_INT64( name, value ) { \
        h5_int64_t h5err = H5PartWriteDataInt64 (H5file_m, name, value); \
    if (h5err <= H5_ERR) { \
        std::stringstream ss; \
        ss << "failed to write int64 data " << name << " to file " << fileName_m; \
        throw GeneralClassicException(std::string(__func__), ss.str()); \
    }\
}
#define SET_STEP() {               \
    h5_int64_t h5err = H5SetStep (H5file_m, H5call_m); \
    if (h5err <= H5_ERR) {                                              \
        std::stringstream ss;                                               \
        ss << "failed to set step " << H5call_m << " in file " << fileName_m; \
        throw GeneralClassicException(std::string(__func__), ss.str()); \
    }\
}
#define GET_NUM_STEPS() {                             \
    H5call_m = H5GetNumSteps( H5file_m );                        \
    if (H5call_m <= H5_ERR) {                                             \
        std::stringstream ss;                                               \
        ss << "failed to get number of steps of file " << fileName_m; \
        throw GeneralClassicException(std::string(__func__), ss.str()); \
    }\
}
#define SET_NUM_PARTICLES( num ) {             \
    h5_int64_t h5err = H5PartSetNumParticles (H5file_m, num); \
    if (h5err <= H5_ERR) {                                              \
        std::stringstream ss;                                               \
        ss << "failed to set number of particles to " << num << " in step " << H5call_m << " in file " << fileName_m; \
        throw GeneralClassicException(std::string(__func__), ss.str()); \
    }\
}

#define OPEN_FILE( fname, mode, props ) {        \
    H5file_m = H5OpenFile (fname, mode, props); \
    if (H5file_m == (h5_file_t)H5_ERR) { \
        std::stringstream ss;                                               \
        ss << "failed to open file " << fileName_m; \
        throw GeneralClassicException(std::string(__func__), ss.str()); \
    }\
}
#define CLOSE_FILE() {               \
    h5_int64_t h5err = H5CloseFile (H5file_m); \
    if (h5err <= H5_ERR) {                                              \
        std::stringstream ss;                                               \
        ss << "failed to close file " << fileName_m; \
        throw GeneralClassicException(std::string(__func__), ss.str()); \
    }\
}

namespace {
    void f64transform(const std::vector<OpalParticle>& particles,
                   unsigned int startIdx,
                   unsigned int numParticles,
                   h5_float64_t *buffer,
                   std::function<h5_float64_t(const OpalParticle&)> select) {
        std::transform(particles.begin() + startIdx,
                       particles.begin() + startIdx + numParticles,
                       buffer, select);
    }
    void i64transform(const std::vector<OpalParticle>& particles,
                   unsigned int startIdx,
                   unsigned int numParticles,
                   h5_int64_t *buffer,
                   std::function<h5_int64_t(const OpalParticle&)> select) {
        std::transform(particles.begin() + startIdx,
                       particles.begin() + startIdx + numParticles,
                       buffer, select);
    }

    void cminmax(double& min, double& max, double val) {
        if (-val > min) {
            min = -val;
        } else if (val > max) {
            max = val;
        }
    }
}

SetStatistics::SetStatistics():
    outputName_m(""),
    spos_m(0.0),
    refTime_m(0.0),
    tmean_m(0.0),
    trms_m(0.0),
    nTotal_m(0),
    RefPartR_m(0.0),
    RefPartP_m(0.0),
    rmin_m(0.0),
    rmax_m(0.0),
    rmean_m(0.0),
    pmean_m(0.0),
    rrms_m(0.0),
    prms_m(0.0),
    rprms_m(0.0),
    normEmit_m(0.0),
    rsqsum_m(0.0),
    psqsum_m(0.0),
    rpsum_m(0.0),
    eps2_m(0.0),
    eps_norm_m(0.0),
    fac_m(0.0)
    { }

LossDataSink::LossDataSink(const std::string& outfn, bool hdf5Save, CollectionType collectionType):
    h5hut_mode_m(hdf5Save),
    H5file_m(0),
    outputName_m(outfn),
    H5call_m(0),
    collectionType_m(collectionType)
{
    particles_m.clear();
    turnNumber_m.clear();
    bunchNumber_m.clear();

    if (h5hut_mode_m && !Options::enableHDF5) {
        throw GeneralClassicException(
                "LossDataSink::LossDataSink",
                "You must select an OPTION to save Loss data files\n"
                "Please, choose 'ENABLEHDF5=TRUE' or 'ASCIIDUMP=TRUE'");
    }

    OpalData::getInstance()->checkAndAddOutputFileName(outputName_m);

    if (OpalData::getInstance()->hasBunchAllocated()) {
        OpalData::getInstance()->setOpenMode(OpalData::OpenMode::APPEND);
    }
}

LossDataSink::LossDataSink(const LossDataSink& rhs):
    h5hut_mode_m(rhs.h5hut_mode_m),
    H5file_m(rhs.H5file_m),
    outputName_m(rhs.outputName_m),
    H5call_m(rhs.H5call_m),
    RefPartR_m(rhs.RefPartR_m),
    RefPartP_m(rhs.RefPartP_m),
    globalTrackStep_m(rhs.globalTrackStep_m),
    refTime_m(rhs.refTime_m),
    spos_m(rhs.spos_m),
    collectionType_m(rhs.collectionType_m)
{
    particles_m.clear();
    turnNumber_m.clear();
    bunchNumber_m.clear();
}

LossDataSink::~LossDataSink() noexcept(false) {
    if (H5file_m) {
        CLOSE_FILE ();
        H5file_m = 0;
    }
}

void LossDataSink::openH5(h5_int32_t mode) {
    h5_prop_t props = H5CreateFileProp ();
    MPI_Comm comm = Ippl::getComm();
    H5SetPropFileMPIOCollective (props, &comm);
    OPEN_FILE (fileName_m.c_str(), mode, props);
    H5CloseProp (props);
}


void LossDataSink::writeHeaderH5() {
    // Write file attributes to describe phase space to H5 file.
    std::stringstream OPAL_version;
    OPAL_version << OPAL_PROJECT_NAME << " " << OPAL_PROJECT_VERSION << " # git rev. " << Util::getGitRevision();
    WRITE_FILEATTRIB_STRING("OPAL_version", OPAL_version.str().c_str());

    WRITE_FILEATTRIB_STRING("SPOSUnit", "m");
    WRITE_FILEATTRIB_STRING("TIMEUnit", "s");
    WRITE_FILEATTRIB_STRING("RefPartRUnit", "m");
    WRITE_FILEATTRIB_STRING("RefPartPUnit", "#beta#gamma");
    WRITE_FILEATTRIB_STRING("GlobalTrackStepUnit", "1");

    WRITE_FILEATTRIB_STRING("centroidUnit", "m");
    WRITE_FILEATTRIB_STRING("RMSXUnit", "m");
    WRITE_FILEATTRIB_STRING("MEANPUnit", "#beta#gamma");
    WRITE_FILEATTRIB_STRING("RMSPUnit", "#beta#gamma");
    WRITE_FILEATTRIB_STRING("#varepsilonUnit", "m rad");
    WRITE_FILEATTRIB_STRING("#varepsilon-geomUnit", "m rad");
    WRITE_FILEATTRIB_STRING("ENERGYUnit", "MeV");
    WRITE_FILEATTRIB_STRING("dEUnit", "MeV");
    WRITE_FILEATTRIB_STRING("TotalChargeUnit", "C");
    WRITE_FILEATTRIB_STRING("TotalMassUnit", "MeV");

    WRITE_FILEATTRIB_STRING("idUnit", "1");
    WRITE_FILEATTRIB_STRING("xUnit", "m");
    WRITE_FILEATTRIB_STRING("yUnit", "m");
    WRITE_FILEATTRIB_STRING("zUnit", "m");
    WRITE_FILEATTRIB_STRING("pxUnit", "#beta#gamma");
    WRITE_FILEATTRIB_STRING("pyUnit", "#beta#gamma");
    WRITE_FILEATTRIB_STRING("pzUnit", "#beta#gamma");
    WRITE_FILEATTRIB_STRING("qUnit", "C");
    WRITE_FILEATTRIB_STRING("mUnit", "MeV");

    WRITE_FILEATTRIB_STRING("turnUnit", "1");
    WRITE_FILEATTRIB_STRING("bunchNumUnit", "1");

    WRITE_FILEATTRIB_STRING("timeUnit", "s");
    WRITE_FILEATTRIB_STRING("meanTimeUnit", "s");
    WRITE_FILEATTRIB_STRING("rmsTimeUnit", "s");

    if (Options::computePercentiles) {
        WRITE_FILEATTRIB_STRING("68-percentileUnit", "m");
        WRITE_FILEATTRIB_STRING("95-percentileUnit", "m");
        WRITE_FILEATTRIB_STRING("99-percentileUnit", "m");
        WRITE_FILEATTRIB_STRING("99_99-percentileUnit", "m");
        WRITE_FILEATTRIB_STRING("normalizedEps68PercentileUnit", "m rad");
        WRITE_FILEATTRIB_STRING("normalizedEps95PercentileUnit", "m rad");
        WRITE_FILEATTRIB_STRING("normalizedEps99PercentileUnit", "m rad");
        WRITE_FILEATTRIB_STRING("normalizedEps99_99PercentileUnit", "m rad");
    }

    if (collectionType_m == CollectionType::TEMPORAL) {
        WRITE_FILEATTRIB_STRING("type", "temporal");
    } else {
        WRITE_FILEATTRIB_STRING("type", "spatial");
    }
}


void LossDataSink::writeHeaderASCII() {
    bool hasTurn = hasTurnInformations();
    if (Ippl::myNode() == 0) {
        os_m << "# x (m),  y (m),  z (m),  px ( ),  py ( ),  pz ( ), id";
        if (hasTurn) {
            os_m << ",  turn ( ), bunchNumber ( )";
        }
        os_m << ", time (s)" << std::endl;
    }
}

void LossDataSink::addReferenceParticle(const Vector_t& x,
                                        const Vector_t& p,
                                        double time,
                                        double spos,
                                        long long globalTrackStep) {
    RefPartR_m.push_back(x);
    RefPartP_m.push_back(p);
    globalTrackStep_m.push_back((h5_int64_t)globalTrackStep);
    spos_m.push_back(spos);
    refTime_m.push_back(time);
}

void LossDataSink::addParticle(const OpalParticle& particle,
                               const boost::optional<std::pair<int, short>>& turnBunchNumPair) {
    if (turnBunchNumPair) {
        if (!particles_m.empty() && turnNumber_m.empty()) {
            throw GeneralClassicException("LossDataSink::addParticle",
                                          "Either no particle or all have turn number and bunch number");
        }
        turnNumber_m.push_back(turnBunchNumPair.get().first);
        bunchNumber_m.push_back(turnBunchNumPair.get().second);
    }

    particles_m.push_back(particle);
}

void LossDataSink::save(unsigned int numSets, OpalData::OpenMode openMode) {

    if (outputName_m.empty()) return;
    if (hasNoParticlesToDump()) return;

    if (openMode == OpalData::OpenMode::UNDEFINED) {
        openMode = OpalData::getInstance()->getOpenMode();
    }

    namespace fs = std::filesystem;
    if (h5hut_mode_m) {
        fileName_m = outputName_m + std::string(".h5");
        if (openMode == OpalData::OpenMode::WRITE || !fs::exists(fileName_m)) {
            openH5();
            writeHeaderH5();
        } else {
            openH5(H5_O_APPENDONLY);
            GET_NUM_STEPS ();
        }

        for (unsigned int i = 0; i < numSets; ++ i) {
            saveH5(i);
        }
        CLOSE_FILE ();
        H5file_m = 0;

    } else {
        fileName_m = outputName_m + std::string(".loss");
        if (openMode == OpalData::OpenMode::WRITE || !fs::exists(fileName_m)) {
            openASCII();
            writeHeaderASCII();
        } else {
            appendASCII();
        }
        saveASCII();
        closeASCII();
    }
    *gmsg << level2 << "Save '" << fileName_m << "'" << endl;

    Ippl::Comm->barrier();

    particles_m.clear();
    turnNumber_m.clear();
    bunchNumber_m.clear();
    spos_m.clear();
    refTime_m.clear();
    RefPartR_m.clear();
    RefPartP_m.clear();
    globalTrackStep_m.clear();
}

// Note: This was changed to calculate the global number of dumped particles
// because there are two cases to be considered:
// 1. ALL nodes have 0 lost particles -> nothing to be done.
// 2. Some nodes have 0 lost particles, some not -> H5 can handle that but all
// nodes HAVE to participate, otherwise H5 waits endlessly for a response from
// the nodes that didn't enter the saveH5 function. -DW
bool LossDataSink::hasNoParticlesToDump() const {

    size_t nLoc = particles_m.size();

    reduce(nLoc, nLoc, OpAddAssign());

    return nLoc == 0;
}

bool LossDataSink::hasTurnInformations() const {

    bool hasTurnInformation = !turnNumber_m.empty();

    allreduce(hasTurnInformation, 1, std::logical_or<bool>());

    return hasTurnInformation > 0;
}

void LossDataSink::saveH5(unsigned int setIdx) {
    size_t nLoc = particles_m.size();
    size_t startIdx = 0, endIdx = nLoc;
    if (setIdx + 1 < startSet_m.size()) {
        startIdx = startSet_m[setIdx];
        endIdx = startSet_m[setIdx + 1];
        nLoc = endIdx - startIdx;
    }

    std::unique_ptr<size_t[]>  locN(new size_t[Ippl::getNodes()]);
    std::unique_ptr<size_t[]> globN(new size_t[Ippl::getNodes()]);

    for (int i = 0; i < Ippl::getNodes(); i++) {
        globN[i] = locN[i] = 0;
    }

    locN[Ippl::myNode()] = nLoc;
    reduce(locN.get(), locN.get() + Ippl::getNodes(), globN.get(), OpAddAssign());

    DistributionMoments engine;
    engine.compute(particles_m.begin() + startIdx, particles_m.begin() + endIdx);

    /// Set current record/time step.
    SET_STEP ();
    SET_NUM_PARTICLES (nLoc);

    if (setIdx < spos_m.size()) {
        WRITE_STEPATTRIB_FLOAT64 ("SPOS", &(spos_m[setIdx]), 1);
        WRITE_STEPATTRIB_FLOAT64 ("TIME", &(refTime_m[setIdx]), 1);
        WRITE_STEPATTRIB_FLOAT64 ("RefPartR", (h5_float64_t *)&(RefPartR_m[setIdx]), 3);
        WRITE_STEPATTRIB_FLOAT64 ("RefPartP", (h5_float64_t *)&(RefPartP_m[setIdx]), 3);
        WRITE_STEPATTRIB_INT64("GlobalTrackStep", &(globalTrackStep_m[setIdx]), 1);
    }

    Vector_t tmpVector;
    double tmpDouble;
    WRITE_STEPATTRIB_FLOAT64("centroid", (tmpVector = engine.getMeanPosition(), &tmpVector[0]), 3);
    WRITE_STEPATTRIB_FLOAT64("RMSX", (tmpVector = engine.getStandardDeviationPosition(), &tmpVector[0]), 3);
    WRITE_STEPATTRIB_FLOAT64("MEANP", (tmpVector = engine.getMeanMomentum(), &tmpVector[0]), 3);
    WRITE_STEPATTRIB_FLOAT64("RMSP", (tmpVector = engine.getStandardDeviationMomentum(), &tmpVector[0]), 3);
    WRITE_STEPATTRIB_FLOAT64("#varepsilon", (tmpVector = engine.getNormalizedEmittance(), &tmpVector[0]), 3);
    WRITE_STEPATTRIB_FLOAT64("#varepsilon-geom", (tmpVector = engine.getGeometricEmittance(), &tmpVector[0]), 3);
    WRITE_STEPATTRIB_FLOAT64("ENERGY", (tmpDouble = engine.getMeanKineticEnergy(), &tmpDouble), 1);
    WRITE_STEPATTRIB_FLOAT64("dE", (tmpDouble = engine.getStdKineticEnergy(), &tmpDouble), 1);
    WRITE_STEPATTRIB_FLOAT64("TotalCharge", (tmpDouble = engine.getTotalCharge(), &tmpDouble), 1);
    WRITE_STEPATTRIB_FLOAT64("TotalMass", (tmpDouble = engine.getTotalMass(), &tmpDouble), 1);
    WRITE_STEPATTRIB_FLOAT64("meanTime", (tmpDouble = engine.getMeanTime(), &tmpDouble), 1);
    WRITE_STEPATTRIB_FLOAT64("rmsTime", (tmpDouble = engine.getStdTime(), &tmpDouble), 1);
    if (Options::computePercentiles) {
        WRITE_STEPATTRIB_FLOAT64("68-percentile", (tmpVector = engine.get68Percentile(), &tmpVector[0]), 3);
        WRITE_STEPATTRIB_FLOAT64("95-percentile", (tmpVector = engine.get95Percentile(), &tmpVector[0]), 3);
        WRITE_STEPATTRIB_FLOAT64("99-percentile", (tmpVector = engine.get99Percentile(), &tmpVector[0]), 3);
        WRITE_STEPATTRIB_FLOAT64("99_99-percentile", (tmpVector = engine.get99_99Percentile(), &tmpVector[0]), 3);

        WRITE_STEPATTRIB_FLOAT64("normalizedEps68Percentile",
                                 (tmpVector = engine.getNormalizedEmittance68Percentile(), &tmpVector[0]), 3);
        WRITE_STEPATTRIB_FLOAT64("normalizedEps95Percentile",
                                 (tmpVector = engine.getNormalizedEmittance95Percentile(), &tmpVector[0]), 3);
        WRITE_STEPATTRIB_FLOAT64("normalizedEps99Percentile",
                                 (tmpVector = engine.getNormalizedEmittance99Percentile(), &tmpVector[0]), 3);
        WRITE_STEPATTRIB_FLOAT64("normalizedEps99_99Percentile",
                                 (tmpVector = engine.getNormalizedEmittance99_99Percentile(), &tmpVector[0]), 3);
    }

    WRITE_STEPATTRIB_FLOAT64("maxR", (tmpVector = engine.getMaxR(), &tmpVector[0]), 3);

    // Write all data
    std::vector<char> buffer(nLoc * sizeof(h5_float64_t));
    h5_float64_t *f64buffer = nLoc > 0? reinterpret_cast<h5_float64_t*>(&buffer[0]) : nullptr;
    h5_int64_t *i64buffer = nLoc > 0? reinterpret_cast<h5_int64_t*>(&buffer[0]) : nullptr;

    ::i64transform(particles_m, startIdx, nLoc, i64buffer,
                   [](const OpalParticle &particle) { return particle.getId(); });
    WRITE_DATA_INT64 ("id", i64buffer);
    ::f64transform(particles_m, startIdx, nLoc, f64buffer,
                   [](const OpalParticle &particle) { return particle.getX(); });
    WRITE_DATA_FLOAT64 ("x", f64buffer);
    ::f64transform(particles_m, startIdx, nLoc, f64buffer,
                   [](const OpalParticle &particle) { return particle.getY(); });
    WRITE_DATA_FLOAT64 ("y", f64buffer);
    ::f64transform(particles_m, startIdx, nLoc, f64buffer,
                   [](const OpalParticle &particle) { return particle.getZ(); });
    WRITE_DATA_FLOAT64 ("z", f64buffer);
    ::f64transform(particles_m, startIdx, nLoc, f64buffer,
                   [](const OpalParticle &particle) { return particle.getPx(); });
    WRITE_DATA_FLOAT64 ("px", f64buffer);
    ::f64transform(particles_m, startIdx, nLoc, f64buffer,
                   [](const OpalParticle &particle) { return particle.getPy(); });
    WRITE_DATA_FLOAT64 ("py", f64buffer);
    ::f64transform(particles_m, startIdx, nLoc, f64buffer,
                   [](const OpalParticle &particle) { return particle.getPz(); });
    WRITE_DATA_FLOAT64 ("pz", f64buffer);
    ::f64transform(particles_m, startIdx, nLoc, f64buffer,
                   [](const OpalParticle &particle) { return particle.getCharge(); });
    WRITE_DATA_FLOAT64 ("q", f64buffer);
    ::f64transform(particles_m, startIdx, nLoc, f64buffer,
                   [](const OpalParticle &particle) { return particle.getMass(); });
    WRITE_DATA_FLOAT64 ("m", f64buffer);

    if (hasTurnInformations()) {
        std::copy(turnNumber_m.begin() + startIdx,
                  turnNumber_m.begin() + startIdx + nLoc,
                  i64buffer);
        WRITE_DATA_INT64 ("turn", i64buffer);

        std::copy(bunchNumber_m.begin() + startIdx,
                  bunchNumber_m.begin() + startIdx + nLoc,
                  i64buffer);
        WRITE_DATA_INT64 ("bunchNumber", i64buffer);
    }

    ::f64transform(particles_m, startIdx, nLoc, f64buffer,
                   [](const OpalParticle &particle) { return particle.getTime(); });
    WRITE_DATA_FLOAT64 ("time", f64buffer);

    ++ H5call_m;
}

void LossDataSink::saveASCII() {
    /*
      ASCII output
    */
    int tag = Ippl::Comm->next_tag(IPPL_APP_TAG3, IPPL_APP_CYCLE);
    bool hasTurn = hasTurnInformations();
    if (Ippl::Comm->myNode() == 0) {
        const unsigned partCount = particles_m.size();

        for (unsigned i = 0; i < partCount; i++) {
            const OpalParticle& particle = particles_m[i];
            os_m << particle.getX() << "   ";
            os_m << particle.getY() << "   ";
            os_m << particle.getZ() << "   ";
            os_m << particle.getPx() << "   ";
            os_m << particle.getPy() << "   ";
            os_m << particle.getPz() << "   ";
            os_m << particle.getId() << "   ";
            if (hasTurn) {
                os_m << turnNumber_m[i] << "   ";
                os_m << bunchNumber_m[i] << "   ";
            }
            os_m << particle.getTime()
                 << std::endl;
        }

        int notReceived =  Ippl::getNodes() - 1;
        while (notReceived > 0) {
            unsigned dataBlocks = 0;
            int node = COMM_ANY_NODE;
            Message *rmsg =  Ippl::Comm->receive_block(node, tag);
            if (rmsg == 0) {
                ERRORMSG("LossDataSink: Could not receive from client nodes output." << endl);
            }
            notReceived--;
            rmsg->get(&dataBlocks);
            rmsg->get(&hasTurn);
            for (unsigned i = 0; i < dataBlocks; i++) {
                long id;
                size_t bunchNum, turn;
                double rx, ry, rz, px, py, pz, time;

                os_m << (rmsg->get(&rx), rx) << "   ";
                os_m << (rmsg->get(&ry), ry) << "   ";
                os_m << (rmsg->get(&rz), rz) << "   ";
                os_m << (rmsg->get(&px), px) << "   ";
                os_m << (rmsg->get(&py), py) << "   ";
                os_m << (rmsg->get(&pz), pz) << "   ";
                os_m << (rmsg->get(&id), id) << "   ";
                if (hasTurn) {
                    os_m << (rmsg->get(&turn), turn) << "   ";
                    os_m << (rmsg->get(&bunchNum), bunchNum) << "   ";
                }
                os_m << (rmsg->get(&time), time)
                     << std::endl;
            }
            delete rmsg;
        }
    } else {
        Message *smsg = new Message();
        const unsigned msgsize = particles_m.size();
        smsg->put(msgsize);
        smsg->put(hasTurn);
        for (unsigned i = 0; i < msgsize; i++) {
            const OpalParticle& particle = particles_m[i];
            smsg->put(particle.getX());
            smsg->put(particle.getY());
            smsg->put(particle.getZ());
            smsg->put(particle.getPx());
            smsg->put(particle.getPy());
            smsg->put(particle.getPz());
            smsg->put(particle.getId());
            if (hasTurn) {
                smsg->put(turnNumber_m[i]);
                smsg->put(bunchNumber_m[i]);
            }
            smsg->put(particle.getTime());
        }
        bool res = Ippl::Comm->send(smsg, 0, tag);
        if (!res) {
            ERRORMSG("LossDataSink Ippl::Comm->send(smsg, 0, tag) failed " << endl);
        }
    }
}

/**
 * In Opal-T monitors can be traversed several times. We know how
* many times the bunch passes because we register the passage of
* the reference particle. This code tries to determine to which
* bunch (same bunch but different times) a particle belongs. For
* this we could use algorithms from data science such as k-means
* or dbscan. But they are an overkill for this application because
* the bunches are well separated.
*
* In a first step we a assign to each bunch the same number of
* particles and compute the mean time of passage and with it a time
* range. Of course this is only an approximation. So we reassign
* the particles to the bunches using the time ranges compute a
* better approximation. Two iterations should be sufficient for
* Opal-T where the temporal separation is large.
*
* @param numSets number of passes of the reference particle
*
*/
void LossDataSink::splitSets(unsigned int numSets) {
    if (numSets <= 1 ||
        particles_m.size() == 0) return;

    const size_t nLoc = particles_m.size();
    size_t avgNumPerSet = nLoc / numSets;
    std::vector<size_t> numPartsInSet(numSets, avgNumPerSet);
    for (unsigned int j = 0; j < (nLoc - numSets * avgNumPerSet); ++ j) {
        ++ numPartsInSet[j];
    }
    startSet_m.resize(numSets + 1, 0);

    std::vector<double> data(2 * numSets, 0.0);
    double* meanT = &data[0];
    double* numParticles = &data[numSets];
    std::vector<double> timeRange(numSets, 0.0);
    double maxT = particles_m[0].getTime();

    for (unsigned int iteration = 0; iteration < 2; ++ iteration) {
        size_t partIdx = 0;
        for (unsigned int j = 0; j < numSets; ++ j) {

            const size_t &numThisSet = numPartsInSet[j];
            for (size_t k = 0; k < numThisSet; ++ k, ++ partIdx) {
                meanT[j] += particles_m[partIdx].getTime();
                maxT = std::max(maxT, particles_m[partIdx].getTime());
            }
            numParticles[j] = numThisSet;
        }

        allreduce(&(data[0]), 2 * numSets, std::plus<double>());

        for (unsigned int j = 0; j < numSets; ++ j) {
            meanT[j] /= numParticles[j];
        }

        for (unsigned int j = 0; j < numSets - 1; ++ j) {
            timeRange[j] = 0.5 * (meanT[j] + meanT[j + 1]);
        }
        timeRange[numSets - 1] = maxT;

        std::fill(numPartsInSet.begin(),
                  numPartsInSet.end(),
                  0);

        size_t setNum = 0;
        size_t idxPrior = 0;
        for (size_t idx = 0; idx < nLoc; ++ idx) {
            if (particles_m[idx].getTime() > timeRange[setNum]) {
                numPartsInSet[setNum] = idx - idxPrior;
                idxPrior = idx;
                ++ setNum;
            }
        }
        numPartsInSet[numSets - 1] = nLoc - idxPrior;
    }

    for (unsigned int i = 0; i < numSets; ++ i) {
        startSet_m[i + 1] = startSet_m[i] + numPartsInSet[i];
    }
}

SetStatistics LossDataSink::computeSetStatistics(unsigned int setIdx) {
    SetStatistics stat;
    double part[6];

    const unsigned int totalSize = 45;
    double plainData[totalSize];
    double rminmax[6] = {};

    Util::KahanAccumulation data[totalSize];
    Util::KahanAccumulation *localCentroid = data + 1;
    Util::KahanAccumulation *localMoments = data + 7;
    Util::KahanAccumulation *localOthers = data + 43;

    size_t startIdx = 0;
    size_t nLoc = particles_m.size();
    if (setIdx + 1 < startSet_m.size()) {
        startIdx = startSet_m[setIdx];
        nLoc = startSet_m[setIdx + 1] - startSet_m[setIdx];
    }

    data[0].sum = nLoc;

    unsigned int idx = startIdx;
    for (unsigned long k = 0; k < nLoc; ++ k, ++ idx) {
        const OpalParticle& particle = particles_m[idx];

        part[0] = particle.getX();
        part[1] = particle.getPx();
        part[2] = particle.getY();
        part[3] = particle.getPy();
        part[4] = particle.getZ();
        part[5] = particle.getPz();

        for (int i = 0; i < 6; i++) {
            localCentroid[i] += part[i];
            for (int j = 0; j <= i; j++) {
                localMoments[i * 6 + j] += part[i] * part[j];
            }
        }
        localOthers[0] += particle.getTime();
        localOthers[1] += std::pow(particle.getTime(), 2);

        ::cminmax(rminmax[0], rminmax[1], particle.getX());
        ::cminmax(rminmax[2], rminmax[3], particle.getY());
        ::cminmax(rminmax[4], rminmax[5], particle.getZ());
    }

    for (int i = 0; i < 6; i++) {
        for (int j = 0; j < i; j++) {
            localMoments[j * 6 + i] = localMoments[i * 6 + j];
        }
    }

    for (unsigned int i = 0; i < totalSize; ++ i) {
        plainData[i] = data[i].sum;
    }

    allreduce(plainData, totalSize, std::plus<double>());
    allreduce(rminmax, 6, std::greater<double>());

    if (plainData[0] == 0.0) return stat;

    double *centroid = plainData + 1;
    double *moments = plainData + 7;
    double *others = plainData + 43;

    stat.outputName_m = outputName_m;
    stat.spos_m = spos_m[setIdx];
    stat.refTime_m = refTime_m[setIdx];
    stat.RefPartR_m = RefPartR_m[setIdx];
    stat.RefPartP_m = RefPartP_m[setIdx];
    stat.nTotal_m = (unsigned long)std::round(plainData[0]);

    for (unsigned int i = 0 ; i < 3u; i++) {
        stat.rmean_m(i) = centroid[2 * i] / stat.nTotal_m;
        stat.pmean_m(i) = centroid[(2 * i) + 1] / stat.nTotal_m;
        stat.rsqsum_m(i) = (moments[2 * i * 6 + 2 * i] -
                            stat.nTotal_m * std::pow(stat.rmean_m(i), 2));
        stat.psqsum_m(i) = std::max(0.0,
                                    moments[(2 * i + 1) * 6 + (2 * i) + 1] -
                                    stat.nTotal_m * std::pow(stat.pmean_m(i), 2));
        stat.rpsum_m(i) = (moments[(2 * i) * 6 + (2 * i) + 1] -
                           stat.nTotal_m * stat.rmean_m(i) * stat.pmean_m(i));
    }
    stat.tmean_m = others[0] / stat.nTotal_m;
    stat.trms_m = std::sqrt(std::max(0.0, (others[1] / stat.nTotal_m - std::pow(stat.tmean_m, 2))));

    stat.eps2_m = ((stat.rsqsum_m * stat.psqsum_m - stat.rpsum_m * stat.rpsum_m) /
                   (1.0 * stat.nTotal_m * stat.nTotal_m));

    stat.rpsum_m /= stat.nTotal_m;

    for (unsigned int i = 0 ; i < 3u; i++) {
        stat.rrms_m(i) = std::sqrt(std::max(0.0, stat.rsqsum_m(i)) / stat.nTotal_m);
        stat.prms_m(i) = std::sqrt(std::max(0.0, stat.psqsum_m(i)) / stat.nTotal_m);
        stat.eps_norm_m(i)  =  std::sqrt(std::max(0.0, stat.eps2_m(i)));
        double tmp = stat.rrms_m(i) * stat.prms_m(i);
        stat.fac_m(i) = (tmp == 0) ? 0.0 : 1.0 / tmp;
        stat.rmin_m(i) = -rminmax[2 * i];
        stat.rmax_m(i) = rminmax[2 * i + 1];
    }

    stat.rprms_m = stat.rpsum_m * stat.fac_m;

    return stat;
}
