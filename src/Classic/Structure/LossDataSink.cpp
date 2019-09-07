#include "Structure/LossDataSink.h"

#include "Ippl.h"
#include "Utilities/Options.h"
#include "AbstractObjects/OpalData.h"
#include "Utilities/GeneralClassicException.h"
#include "Utilities/Util.h"

#include <boost/filesystem.hpp>
#include "OPALconfig.h"

#include <cassert>

#define ADD_ATTACHMENT( fname ) {             \
    h5_int64_t h5err = H5AddAttachment (H5file_m, fname); \
    if (h5err <= H5_ERR) {                                              \
        std::stringstream ss;                                               \
        ss << "failed to add attachment " << fname << " to file " << fn_m; \
        throw GeneralClassicException(std::string(__func__), ss.str()); \
    }\
}
#define WRITE_FILEATTRIB_STRING( attribute, value ) {             \
    h5_int64_t h5err = H5WriteFileAttribString (H5file_m, attribute, value); \
    if (h5err <= H5_ERR) {                                              \
        std::stringstream ss;                                               \
        ss << "failed to write string attribute " << attribute << " to file " << fn_m; \
        throw GeneralClassicException(std::string(__func__), ss.str()); \
    }\
}
#define WRITE_STEPATTRIB_FLOAT64( attribute, value, size ) { \
        h5_int64_t h5err = H5WriteStepAttribFloat64 (H5file_m, attribute, value, size); \
    if (h5err <= H5_ERR) { \
        std::stringstream ss; \
        ss << "failed to write float64 attribute " << attribute << " to file " << fn_m; \
        throw GeneralClassicException(std::string(__func__), ss.str()); \
    }\
}
#define WRITE_STEPATTRIB_INT64( attribute, value, size ) { \
        h5_int64_t h5err = H5WriteStepAttribInt64 (H5file_m, attribute, value, size); \
    if (h5err <= H5_ERR) { \
        std::stringstream ss; \
        ss << "failed to write int64 attribute " << attribute << " to file " << fn_m; \
        throw GeneralClassicException(std::string(__func__), ss.str()); \
    }\
}
#define WRITE_DATA_FLOAT64( name, value ) { \
        h5_int64_t h5err = H5PartWriteDataFloat64 (H5file_m, name, value); \
    if (h5err <= H5_ERR) { \
        std::stringstream ss; \
        ss << "failed to write float64 data " << name << " to file " << fn_m; \
        throw GeneralClassicException(std::string(__func__), ss.str()); \
    }\
}
#define WRITE_DATA_INT64( name, value ) { \
        h5_int64_t h5err = H5PartWriteDataInt64 (H5file_m, name, value); \
    if (h5err <= H5_ERR) { \
        std::stringstream ss; \
        ss << "failed to write int64 data " << name << " to file " << fn_m; \
        throw GeneralClassicException(std::string(__func__), ss.str()); \
    }\
}
#define SET_STEP() {               \
    h5_int64_t h5err = H5SetStep (H5file_m, H5call_m); \
    if (h5err <= H5_ERR) {                                              \
        std::stringstream ss;                                               \
        ss << "failed to set step " << H5call_m << " in file " << fn_m; \
        throw GeneralClassicException(std::string(__func__), ss.str()); \
    }\
}
#define GET_NUM_STEPS() {                             \
    H5call_m = H5GetNumSteps( H5file_m );                        \
    if (H5call_m <= H5_ERR) {                                             \
        std::stringstream ss;                                               \
        ss << "failed to get number of steps of file " << fn_m; \
        throw GeneralClassicException(std::string(__func__), ss.str()); \
    }\
}
#define SET_NUM_PARTICLES( num ) {             \
    h5_int64_t h5err = H5PartSetNumParticles (H5file_m, num); \
    if (h5err <= H5_ERR) {                                              \
        std::stringstream ss;                                               \
        ss << "failed to set number of particles to " << num << " in step " << H5call_m << " in file " << fn_m; \
        throw GeneralClassicException(std::string(__func__), ss.str()); \
    }\
}

#define OPEN_FILE( fname, mode, props ) {        \
    H5file_m = H5OpenFile (fname, mode, props); \
    if(H5file_m == (h5_file_t)H5_ERR) { \
        std::stringstream ss;                                               \
        ss << "failed to open file " << fn_m; \
        throw GeneralClassicException(std::string(__func__), ss.str()); \
    }\
}
#define CLOSE_FILE() {               \
    h5_int64_t h5err = H5CloseFile (H5file_m); \
    if (h5err <= H5_ERR) {                                              \
        std::stringstream ss;                                               \
        ss << "failed to close file " << fn_m; \
        throw GeneralClassicException(std::string(__func__), ss.str()); \
    }\
}


SetStatistics::SetStatistics():
    element_m(""),
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

LossDataSink::LossDataSink(std::string elem, bool hdf5Save, ElementBase::ElementType type):
    h5hut_mode_m(hdf5Save),
    H5file_m(0),
    element_m(elem),
    H5call_m(0),
    type_m(type)
{
    x_m.clear();
    y_m.clear();
    z_m.clear();
    px_m.clear();
    py_m.clear();
    pz_m.clear();
    id_m.clear();
    turn_m.clear();
    bunchNum_m.clear();
    time_m.clear();
}

LossDataSink::LossDataSink(const LossDataSink &rsh):
    h5hut_mode_m(rsh.h5hut_mode_m),
    H5file_m(rsh.H5file_m),
    element_m(rsh.element_m),
    H5call_m(rsh.H5call_m),
    RefPartR_m(rsh.RefPartR_m),
    RefPartP_m(rsh.RefPartP_m),
    globalTrackStep_m(rsh.globalTrackStep_m),
    refTime_m(rsh.refTime_m),
    spos_m(rsh.spos_m),
    type_m(rsh.type_m)
{
    x_m.clear();
    y_m.clear();
    z_m.clear();
    px_m.clear();
    py_m.clear();
    pz_m.clear();
    id_m.clear();
    turn_m.clear();
    bunchNum_m.clear();
    time_m.clear();
}

LossDataSink::LossDataSink() {
    LossDataSink(std::string("NULL"), false);
}

LossDataSink::~LossDataSink() noexcept(false) {
    if (H5file_m) {
        CLOSE_FILE ();
        H5file_m = 0;
    }
    Ippl::Comm->barrier();
}

void LossDataSink::openH5(h5_int32_t mode) {
    h5_prop_t props = H5CreateFileProp ();
    MPI_Comm comm = Ippl::getComm();
    H5SetPropFileMPIOCollective (props, &comm);
    OPEN_FILE (fn_m.c_str(), mode, props);
    H5CloseProp (props);
}


void LossDataSink::writeHeaderH5() {

    // Write file attributes to describe phase space to H5 file.
    std::stringstream OPAL_version;
    OPAL_version << OPAL_PROJECT_NAME << " " << OPAL_PROJECT_VERSION << " # git rev. " << Util::getGitRevision();
    WRITE_FILEATTRIB_STRING ("OPAL_version", OPAL_version.str().c_str());
    WRITE_FILEATTRIB_STRING ("tUnit", "s");
    WRITE_FILEATTRIB_STRING ("xUnit", "m");
    WRITE_FILEATTRIB_STRING ("yUnit", "m");
    WRITE_FILEATTRIB_STRING ("zUnit", "m");
    WRITE_FILEATTRIB_STRING ("pxUnit", "#beta#gamma");
    WRITE_FILEATTRIB_STRING ("pyUnit", "#beta#gamma");
    WRITE_FILEATTRIB_STRING ("pzUnit", "#beta#gamma");
    WRITE_FILEATTRIB_STRING ("idUnit", "1");

    /// in case of circular machines
    if (hasTimeAttribute()) {
        WRITE_FILEATTRIB_STRING ("turnUnit", "1");
        WRITE_FILEATTRIB_STRING ("timeUnit", "s");
    }
    WRITE_FILEATTRIB_STRING ("SPOSUnit", "mm");
    WRITE_FILEATTRIB_STRING ("TIMEUnit", "s");

    WRITE_FILEATTRIB_STRING ("mpart", "GeV");
    WRITE_FILEATTRIB_STRING ("qi", "C");

    ADD_ATTACHMENT (OpalData::getInstance()->getInputFn().c_str());
}

void LossDataSink::addReferenceParticle(const Vector_t &x,
                                        const  Vector_t &p,
                                        double time,
                                        double spos,
                                        long long globalTrackStep
                                        ) {
    RefPartR_m.push_back(x);
    RefPartP_m.push_back(p);
    globalTrackStep_m.push_back((h5_int64_t)globalTrackStep);
    spos_m.push_back(spos);
    refTime_m.push_back(time);
}

void LossDataSink::addParticle(const Vector_t &x,const  Vector_t &p, const size_t  id) {
    x_m.push_back(x(0));
    y_m.push_back(x(1));
    z_m.push_back(x(2));
    px_m.push_back(p(0));
    py_m.push_back(p(1));
    pz_m.push_back(p(2));
    id_m.push_back(id);
}

// For ring type simulation, dump the time and turn number
void LossDataSink::addParticle(const Vector_t &x, const Vector_t &p, const size_t id,
                               const double time, const size_t turn,
                               const size_t& bunchNum) {
    addParticle(x, p, id);
    turn_m.push_back(turn);
    time_m.push_back(time);
    bunchNum_m.push_back(bunchNum);
}

void LossDataSink::save(unsigned int numSets, OpalData::OPENMODE openMode) {

    if (element_m == std::string("")) return;
    if (hasNoParticlesToDump()) return;

    if (openMode == OpalData::OPENMODE::UNDEFINED) {
        openMode = OpalData::getInstance()->getOpenMode();
    }

    namespace fs = boost::filesystem;
    if (h5hut_mode_m) {
        if (!Options::enableHDF5) return;

        fn_m = element_m + std::string(".h5");
        INFOMSG(level2 << "Save " << fn_m << endl);
        if (openMode == OpalData::OPENMODE::WRITE || !fs::exists(fn_m)) {
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
    }
    else {
        fn_m = element_m + std::string(".loss");
        INFOMSG(level2 << "Save " << fn_m << endl);
        if (openMode == OpalData::OPENMODE::WRITE || !fs::exists(fn_m)) {
            openASCII();
            writeHeaderASCII();
        } else {
            appendASCII();
        }
        saveASCII();
        closeASCII();
    }
    Ippl::Comm->barrier();

    x_m.clear();
    y_m.clear();
    z_m.clear();
    px_m.clear();
    py_m.clear();
    pz_m.clear();
    id_m.clear();
    turn_m.clear();
    bunchNum_m.clear();
    time_m.clear();
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
bool LossDataSink::hasNoParticlesToDump() {

    size_t nLoc = x_m.size();

    reduce(nLoc, nLoc, OpAddAssign());

    return nLoc == 0;
}

bool LossDataSink::hasTimeAttribute() {

    size_t tLoc = time_m.size();

    reduce(tLoc, tLoc, OpAddAssign());

    return tLoc > 0;
}

void LossDataSink::saveH5(unsigned int setIdx) {
    size_t startIdx = 0;
    size_t nLoc = x_m.size();
    if (setIdx + 1 < startSet_m.size()) {
        startIdx = startSet_m[setIdx];
        nLoc = startSet_m[setIdx + 1] - startSet_m[setIdx];
    }

    std::unique_ptr<size_t[]>  locN(new size_t[Ippl::getNodes()]);
    std::unique_ptr<size_t[]> globN(new size_t[Ippl::getNodes()]);

    for(int i = 0; i < Ippl::getNodes(); i++) {
        globN[i] = locN[i] = 0;
    }

    locN[Ippl::myNode()] = nLoc;
    reduce(locN.get(), locN.get() + Ippl::getNodes(), globN.get(), OpAddAssign());

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
    // Write all data
    WRITE_DATA_FLOAT64 ("x", &x_m[startIdx]);
    WRITE_DATA_FLOAT64 ("y", &y_m[startIdx]);
    WRITE_DATA_FLOAT64 ("z", &z_m[startIdx]);
    WRITE_DATA_FLOAT64 ("px", &px_m[startIdx]);
    WRITE_DATA_FLOAT64 ("py", &py_m[startIdx]);
    WRITE_DATA_FLOAT64 ("pz", &pz_m[startIdx]);

    /// Write particle id numbers.
    std::vector<h5_int64_t> larray(id_m.begin() + startIdx, id_m.end() );
    WRITE_DATA_INT64 ("id", &larray[0]);
    if (hasTimeAttribute()) {
        WRITE_DATA_FLOAT64 ("time", &time_m[startIdx]);
        larray.assign (turn_m.begin() + startIdx, turn_m.end() );
        WRITE_DATA_INT64 ("turn", &larray[0]);
        larray.assign (bunchNum_m.begin() + startIdx, bunchNum_m.end() );
        WRITE_DATA_INT64 ("bunchNumber", &larray[0]);
    }

    ++ H5call_m;
}

void LossDataSink::saveASCII() {

    /*
      ASCII output
    */
    int tag = Ippl::Comm->next_tag(IPPL_APP_TAG3, IPPL_APP_CYCLE);
    bool hasTime = hasTimeAttribute(); // reduce needed for the case when node 0 has no particles
    if(Ippl::Comm->myNode() == 0) {
        const unsigned partCount = x_m.size();

        for(unsigned i = 0; i < partCount; i++) {
            os_m << element_m   << "   ";
            os_m << x_m[i] << "   ";
            os_m << y_m[i] << "   ";
            os_m << z_m[i] << "   ";
            os_m << px_m[i] << "   ";
            os_m << py_m[i] << "   ";
            os_m << pz_m[i] << "   ";
            os_m << id_m[i]   << "   ";
            if (hasTime) {
                os_m << turn_m[i] << "   ";
                os_m << bunchNum_m[i] << "   ";
                os_m << time_m[i];
            }
            os_m << std::endl;
        }

        int notReceived =  Ippl::getNodes() - 1;
        while(notReceived > 0) {
            unsigned dataBlocks = 0;
            int node = COMM_ANY_NODE;
            Message *rmsg =  Ippl::Comm->receive_block(node, tag);
            if(rmsg == 0) {
                ERRORMSG("LossDataSink: Could not receive from client nodes output." << endl);
            }
            notReceived--;
            rmsg->get(&dataBlocks);
            for(unsigned i = 0; i < dataBlocks; i++) {
                long id;
                size_t bunchNum, turn;
                double rx, ry, rz, px, py, pz, time;
                rmsg->get(&id);
                rmsg->get(&rx);
                rmsg->get(&ry);
                rmsg->get(&rz);
                rmsg->get(&px);
                rmsg->get(&py);
                rmsg->get(&pz);
                if (hasTime) {
                    rmsg->get(&turn);
                    rmsg->get(&bunchNum);
                    rmsg->get(&time);
                }
                os_m << element_m << "   ";
                os_m << rx << "   ";
                os_m << ry << "   ";
                os_m << rz << "   ";
                os_m << px << "   ";
                os_m << py << "   ";
                os_m << pz << "   ";
                os_m << id << "   ";
                if (hasTime) {
                    os_m << turn << "   ";
                    os_m << bunchNum << "   ";
                    os_m << time;
                }
                os_m << std::endl;
            }
            delete rmsg;
        }
    } else {
        Message *smsg = new Message();
        const unsigned msgsize = x_m.size();
        smsg->put(msgsize);
        for(unsigned i = 0; i < msgsize; i++) {
            smsg->put(id_m[i]);
            smsg->put(x_m[i]);
            smsg->put(y_m[i]);
            smsg->put(z_m[i]);
            smsg->put(px_m[i]);
            smsg->put(py_m[i]);
            smsg->put(pz_m[i]);
            if (hasTime) {
                smsg->put(turn_m[i]);
                smsg->put(bunchNum_m[i]);
                smsg->put(time_m[i]);
            }
        }
        bool res = Ippl::Comm->send(smsg, 0, tag);
        if(! res)
            ERRORMSG("LossDataSink Ippl::Comm->send(smsg, 0, tag) failed " << endl;);
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
        x_m.size() == 0 ||
        time_m.size() != x_m.size()) return;

    const size_t nLoc = x_m.size();
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
    double maxT = time_m[0];

    for (unsigned int iteration = 0; iteration < 2; ++ iteration) {
        size_t partIdx = 0;
        for (unsigned int j = 0; j < numSets; ++ j) {

            const size_t &numThisSet = numPartsInSet[j];
            for (size_t k = 0; k < numThisSet; ++ k, ++ partIdx) {
                meanT[j] += time_m[partIdx];
                maxT = std::max(maxT, time_m[partIdx]);
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
            if (time_m[idx] > timeRange[setNum]) {
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

namespace {
    void cminmax(double &min, double &max, double val) {
        if (-val > min) {
            min = -val;
        } else if (val > max) {
            max = val;
        }
    }
}

SetStatistics LossDataSink::computeSetStatistics(unsigned int setIdx) {
    SetStatistics stat;
    double part[6];

    const unsigned int totalSize = 45;
    double plainData[totalSize];
    double rminmax[6];

    Util::KahanAccumulation data[totalSize];
    Util::KahanAccumulation *localCentroid = data + 1;
    Util::KahanAccumulation *localMoments = data + 7;
    Util::KahanAccumulation *localOthers = data + 43;

    size_t startIdx = 0;
    size_t nLoc = x_m.size();
    if (setIdx + 1 < startSet_m.size()) {
        startIdx = startSet_m[setIdx];
        nLoc = startSet_m[setIdx + 1] - startSet_m[setIdx];
    }

    data[0].sum = nLoc;

    unsigned int idx = startIdx;
    for(unsigned long k = 0; k < nLoc; ++ k, ++ idx) {
        part[1] = px_m[idx];
        part[3] = py_m[idx];
        part[5] = pz_m[idx];
        part[0] = x_m[idx];
        part[2] = y_m[idx];
        part[4] = z_m[idx];

        for(int i = 0; i < 6; i++) {
            localCentroid[i] += part[i];
            for(int j = 0; j <= i; j++) {
                localMoments[i * 6 + j] += part[i] * part[j];
            }
        }
        localOthers[0] += time_m[idx];
        localOthers[1] += std::pow(time_m[idx], 2);

        ::cminmax(rminmax[0], rminmax[1], x_m[idx]);
        ::cminmax(rminmax[2], rminmax[3], y_m[idx]);
        ::cminmax(rminmax[4], rminmax[5], z_m[idx]);
    }

    for(int i = 0; i < 6; i++) {
        for(int j = 0; j < i; j++) {
            localMoments[j * 6 + i] = localMoments[i * 6 + j];
        }
    }

    for (unsigned int i = 0; i < totalSize; ++ i) {
        plainData[i] = data[i].sum;
    }

    new_reduce(plainData, totalSize, std::plus<double>());
    new_reduce(rminmax, 6, std::greater<double>());

    if (Ippl::myNode() != 0) return stat;
    if (plainData[0] == 0.0) return stat;

    double *centroid = plainData + 1;
    double *moments = plainData + 7;
    double *others = plainData + 43;

    stat.element_m = element_m;
    stat.spos_m = spos_m[setIdx];
    stat.refTime_m = refTime_m[setIdx];
    stat.RefPartR_m = RefPartR_m[setIdx];
    stat.RefPartP_m = RefPartP_m[setIdx];
    stat.nTotal_m = (unsigned long)floor(plainData[0] + 0.5);

    for(unsigned int i = 0 ; i < 3u; i++) {
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
    stat.trms_m = sqrt(std::max(0.0, (others[1] / stat.nTotal_m - std::pow(stat.tmean_m, 2))));

    stat.eps2_m = ((stat.rsqsum_m * stat.psqsum_m - stat.rpsum_m * stat.rpsum_m) /
                   (1.0 * stat.nTotal_m * stat.nTotal_m));

    stat.rpsum_m /= stat.nTotal_m;

    for(unsigned int i = 0 ; i < 3u; i++) {
        stat.rrms_m(i) = sqrt(std::max(0.0, stat.rsqsum_m(i)) / stat.nTotal_m);
        stat.prms_m(i) = sqrt(std::max(0.0, stat.psqsum_m(i)) / stat.nTotal_m);
        stat.eps_norm_m(i)  =  std::sqrt(std::max(0.0, stat.eps2_m(i)));
        double tmp = stat.rrms_m(i) * stat.prms_m(i);
        stat.fac_m(i) = (tmp == 0) ? 0.0 : 1.0 / tmp;
        stat.rmin_m(i) = -rminmax[2 * i];
        stat.rmax_m(i) = rminmax[2 * i + 1];
    }

    stat.rprms_m = stat.rpsum_m * stat.fac_m;

    return stat;
}


// vi: set et ts=4 sw=4 sts=4:
// Local Variables:
// mode:c++
// c-basic-offset: 4
// indent-tabs-mode:nil
// End: