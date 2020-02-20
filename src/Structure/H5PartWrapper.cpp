//
//  Copyright & License: See Copyright.readme in src directory
//

#include "Structure/H5PartWrapper.h"

#include "OPALconfig.h"
#include "AbstractObjects/OpalData.h"
#include "Utilities/Options.h"
#include "Physics/Physics.h"

#include <boost/filesystem.hpp>

#include <fstream>

namespace {
    const h5_int64_t H5TypesCHAR = H5_STRING_T;
    const h5_int64_t H5TypesFLOAT = H5_FLOAT32_T;
    const h5_int64_t H5TypesDOUBLE = H5_FLOAT64_T;
    const h5_int64_t H5TypesINT32 = H5_INT32_T;
    const h5_int64_t H5TypesINT64 = H5_INT64_T;
}

std::string H5PartWrapper::copyFilePrefix_m = ".copy";

H5PartWrapper::H5PartWrapper(const std::string &fileName, h5_int32_t flags):
    file_m(0),
    fileName_m(fileName),
    predecessorOPALFlavour_m("NOT SET"),
    numSteps_m(0),
    startedFromExistingFile_m(false)
{
    open(flags);
}

H5PartWrapper::H5PartWrapper(const std::string &fileName, int restartStep, std::string sourceFile, h5_int32_t flags):
    file_m(0),
    fileName_m(fileName),
    predecessorOPALFlavour_m("NOT SET"),
    numSteps_m(0),
    startedFromExistingFile_m(true)
{
    namespace fs = boost::filesystem;

    if (sourceFile == "") sourceFile = fileName_m;

    copyFile(sourceFile, restartStep, flags);

    open(H5_O_RDWR);
}

H5PartWrapper::~H5PartWrapper() {
    close();
}

void H5PartWrapper::close() {
    if (file_m) {
        Ippl::Comm->barrier();

        REPORTONERROR(H5CloseFile(file_m));

        file_m = 0;
    }
}

void H5PartWrapper::open(h5_int32_t flags) {
    close();

    h5_prop_t props = H5CreateFileProp ();
    MPI_Comm comm = Ippl::getComm();
    h5_err_t h5err = H5SetPropFileMPIOCollective (props, &comm);
#if defined (NDEBUG)
    (void)h5err;
#endif
    assert (h5err != H5_ERR);
    file_m = H5OpenFile (fileName_m.c_str(), flags, props);
    assert (file_m != (h5_file_t)H5_ERR);
    H5CloseProp (props);
}

void H5PartWrapper::storeCavityInformation() {
    /// Write number of Cavities with autophase information
    h5_int64_t nAutoPhaseCavities = OpalData::getInstance()->getNumberOfMaxPhases();
    h5_int64_t nFormerlySavedAutoPhaseCavities = 0;
    bool fileWasClosed = (file_m == 0);

    if (nAutoPhaseCavities == 0) return;
    if (fileWasClosed) open(H5_O_APPENDONLY);
    if (!H5HasFileAttrib(file_m, "nAutoPhaseCavities") ||
        H5ReadFileAttribInt64(file_m, "nAutoPhaseCavities", &nFormerlySavedAutoPhaseCavities) != H5_SUCCESS) {
        nFormerlySavedAutoPhaseCavities = 0;
    }
    if (nFormerlySavedAutoPhaseCavities == nAutoPhaseCavities) {
        if (fileWasClosed) close();
        return;
    }

    WRITEFILEATTRIB(Int64, file_m, "nAutoPhaseCavities", &nAutoPhaseCavities, 1);

    unsigned int elementNumber = 1;
    std::vector<MaxPhasesT>::iterator it = OpalData::getInstance()->getFirstMaxPhases();
    std::vector<MaxPhasesT>::iterator end = OpalData::getInstance()->getLastMaxPhases();
    for(; it < end; ++ it, ++ elementNumber) {
        if (elementNumber <= nFormerlySavedAutoPhaseCavities) continue;

        std::string nameAttributeName = "Cav-" + std::to_string(elementNumber) + "-name";
        std::string valueAttributeName  = "Cav-" + std::to_string(elementNumber) + "-value";

        std::string elementName   = (*it).first;
        h5_float64_t elementPhase = (*it).second;

        WRITESTRINGFILEATTRIB(file_m, nameAttributeName.c_str(), elementName.c_str());
        WRITEFILEATTRIB(Float64, file_m, valueAttributeName.c_str(), &elementPhase, 1);

        INFOMSG("Saved phases in the h5 file: "
                << nameAttributeName << " -> " << elementName << " --- "
                << valueAttributeName << " -> " << elementPhase << endl);
    }

    if (fileWasClosed) close();
}

void H5PartWrapper::copyFile(const std::string &sourceFile, int lastStep, h5_int32_t flags) {

    namespace fs = boost::filesystem;
    if (!fs::exists(sourceFile)) {
        throw OpalException("H5PartWrapper::copyFile",
                            "source file '" + sourceFile + "' does not exist");
    }

    if (sourceFile == fileName_m) {
        h5_prop_t props = H5CreateFileProp ();
        MPI_Comm comm = Ippl::getComm();
        h5_err_t h5err = H5SetPropFileMPIOCollective (props, &comm);
#if defined (NDEBUG)
        (void)h5err;
#endif
        assert (h5err != H5_ERR);
        h5_file_t source = H5OpenFile (sourceFile.c_str(), H5_O_RDONLY, props);
        assert (source != (h5_file_t)H5_ERR);
        H5CloseProp (props);
        h5_ssize_t numStepsInSource = H5GetNumSteps(source);

        if (lastStep == -1 || lastStep >= numStepsInSource) {
            REPORTONERROR(H5SetStep(source, numStepsInSource - 1));

            char opalFlavour[128];
            READSTEPATTRIB(String, source, "OPAL_flavour", opalFlavour);
            predecessorOPALFlavour_m = std::string(opalFlavour);

            REPORTONERROR(H5CloseFile(source));

            numSteps_m = numStepsInSource;
            return;
        }

        REPORTONERROR(H5CloseFile(source));

        Ippl::Comm->barrier();

        std::string sourceFileName = copyFilePrefix_m + fileName_m;
        if (Ippl::myNode() == 0) {
            fs::rename(fileName_m, sourceFileName);
        }

        Ippl::Comm->barrier();

        open(flags);
        props = H5CreateFileProp ();
        comm = Ippl::getComm();
        h5err = H5SetPropFileMPIOCollective (props, &comm);
        assert (h5err != H5_ERR);
        source = H5OpenFile (sourceFileName.c_str(), H5_O_RDONLY, props);
        assert (source != (h5_file_t)H5_ERR);
        H5CloseProp (props);
        copyHeader(source);

        if (lastStep < 0) {
            if (-lastStep > numStepsInSource) {
                lastStep = 0;
            } else {
                lastStep = numStepsInSource + lastStep;
            }
        }

        // don't copy the whole file, it takes very long
        copyStep(source, lastStep);
        ++ numSteps_m;

        REPORTONERROR(H5CloseFile(source));

        if (Ippl::myNode() == 0) {
            fs::remove(sourceFileName);
        }

        close();
    } else {

        open(flags);

        h5_prop_t props = H5CreateFileProp ();
        MPI_Comm comm = Ippl::getComm();
        h5_err_t h5err = H5SetPropFileMPIOCollective (props, &comm);
#if defined (NDEBUG)
        (void)h5err;
#endif
        assert (h5err != H5_ERR);
        h5_file_t source = H5OpenFile (sourceFile.c_str(), H5_O_RDONLY, props);
        assert (source != (h5_file_t)H5_ERR);
        H5CloseProp (props);
        h5_ssize_t numStepsInSource = H5GetNumSteps(source);

        if (lastStep == -1 || lastStep >= numStepsInSource) {
            REPORTONERROR(H5SetStep(source, numStepsInSource - 1));

            char opalFlavour[128];
            READSTEPATTRIB(String, source, "OPAL_flavour", opalFlavour);
            predecessorOPALFlavour_m = std::string(opalFlavour);

            REPORTONERROR(H5CloseFile(source));
            close();

            copyFileSystem(sourceFile);

            numSteps_m = numStepsInSource;

        } else {
            copyHeader(source);

            if (lastStep < 0) {
                if (-lastStep > numStepsInSource) {
                    lastStep = 0;
                } else {
                    lastStep = numStepsInSource + lastStep;
                }
            }
            // don't copy the whole file, it takes very long
            copyStep(source, lastStep);
            ++ numSteps_m;

            REPORTONERROR(H5CloseFile(source));
        }

        close();
    }
}

void H5PartWrapper::copyFileSystem(const std::string &sourceFile) {
    // namespace fs = boost::filesystem;

    if (sourceFile == fileName_m) return;

    int sourceNode = 0;
    if (Ippl::myNode() == sourceNode) {

        // copy_file not working due to bug in boost, see
        // https://svn.boost.org/trac/boost/ticket/10038
        // try {
            // fs::copy_file(sourceFile, fileName_m, fs::copy_option::none);
        // } catch (fs::filesystem_error &m) {

            // ERRORMSG(m.what() << endl);

        std::ifstream source(sourceFile, std::ios::binary);
        std::ofstream dest(fileName_m, std::ios::binary);

        std::istreambuf_iterator<char> begin_source(source);
        std::istreambuf_iterator<char> end_source;
        std::ostreambuf_iterator<char> begin_dest(dest);
        std::copy(begin_source, end_source, begin_dest);

        source.close();

        sendFailureMessage(dest.bad(),
                           "H5PartWrapper::copyFile",
                           "could not copy file " + sourceFile);
        dest.close();
    } else {
        receiveFailureMessage(sourceNode,
                              "H5PartWrapper::copyFile",
                              "received message to throw exception from node 0");
    }
}

void H5PartWrapper::copyHeader(
    h5_file_t source
    ) {
    h5_int64_t numFileAttributes = H5GetNumFileAttribs(source);

    const h5_size_t lengthAttributeName = 256;
    char attributeName[lengthAttributeName];
    h5_int64_t attributeType;
    h5_size_t numAttributeElements;
    std::vector<char> buffer(256);
    h5_float32_t *f32buffer = reinterpret_cast<h5_float32_t*>(&buffer[0]);
    h5_float64_t *f64buffer = reinterpret_cast<h5_float64_t*>(&buffer[0]);
    h5_int32_t *i32buffer = reinterpret_cast<h5_int32_t*>(&buffer[0]);
    h5_int64_t *i64buffer = reinterpret_cast<h5_int64_t*>(&buffer[0]);

    for (h5_int64_t i = 0; i < numFileAttributes; ++ i) {
        REPORTONERROR(H5GetFileAttribInfo(source,
                                          i,
                                          attributeName,
                                          lengthAttributeName,
                                          &attributeType,
                                          &numAttributeElements));

        if (attributeType == H5_STRING_T) {
            if (buffer.size() < numAttributeElements) {
                buffer.resize(numAttributeElements);
            }

            READFILEATTRIB(String, source, attributeName, &buffer[0]);
            WRITESTRINGFILEATTRIB(file_m, attributeName, &buffer[0]);

        } else if (attributeType == H5_FLOAT32_T) {
            if (buffer.size() < numAttributeElements * sizeof(h5_float32_t)) {
                buffer.resize(numAttributeElements * sizeof(h5_float32_t));
            }

            READFILEATTRIB(Float32, source, attributeName, f32buffer);
            WRITEFILEATTRIB(Float32, file_m, attributeName, f32buffer, numAttributeElements);

        } else if (attributeType == H5_FLOAT64_T) {
            if (buffer.size() < numAttributeElements * sizeof(h5_float64_t)) {
                buffer.resize(numAttributeElements * sizeof(h5_float64_t));
            }

            READFILEATTRIB(Float64, source, attributeName, f64buffer);
            WRITEFILEATTRIB(Float64, file_m, attributeName, f64buffer, numAttributeElements);

        } else if (attributeType == H5_INT32_T) {
            if (buffer.size() < numAttributeElements * sizeof(h5_int32_t)) {
                buffer.resize(numAttributeElements * sizeof(h5_int32_t));
            }

            READFILEATTRIB(Int32, source, attributeName, i32buffer);
            WRITEFILEATTRIB(Int32, file_m, attributeName, i32buffer, numAttributeElements);

        } else if (attributeType == H5_INT64_T) {
            if (buffer.size() < numAttributeElements * sizeof(h5_int64_t)) {
                buffer.resize(numAttributeElements * sizeof(h5_int64_t));
            }

            READFILEATTRIB(Int64, source, attributeName, i64buffer);
            WRITEFILEATTRIB(Int64, file_m, attributeName, i64buffer, numAttributeElements);

        } else {
            throw OpalException("H5PartWrapper::copyHeader",
                                "unknown data type: " + std::to_string(attributeType));
        }
    }
}

void H5PartWrapper::copyStep(
    h5_file_t source,
    int step
    ) {
    REPORTONERROR(H5SetStep(file_m, numSteps_m));
    REPORTONERROR(H5SetStep(source, step));

    copyStepHeader(source);
    copyStepData(source);
}

void H5PartWrapper::copyStepHeader(
    h5_file_t source
    ) {
    h5_int64_t numStepAttributes = H5GetNumStepAttribs(source);

    h5_size_t lengthAttributeName = 256;
    char attributeName[lengthAttributeName];
    h5_int64_t attributeType;
    h5_size_t numAttributeElements;

    std::vector<char> buffer(256);
    h5_float32_t *f32buffer = reinterpret_cast<h5_float32_t*>(&buffer[0]);
    h5_float64_t *f64buffer = reinterpret_cast<h5_float64_t*>(&buffer[0]);
    h5_int32_t *i32buffer = reinterpret_cast<h5_int32_t*>(&buffer[0]);
    h5_int64_t *i64buffer = reinterpret_cast<h5_int64_t*>(&buffer[0]);

    READSTEPATTRIB(String, source, "OPAL_flavour", &buffer[0]);
    predecessorOPALFlavour_m = std::string(&buffer[0]);

    for (h5_int64_t i = 0; i < numStepAttributes; ++ i) {
        REPORTONERROR(H5GetStepAttribInfo(source,
                                          i,
                                          attributeName,
                                          lengthAttributeName,
                                          &attributeType,
                                          &numAttributeElements));

        if (attributeType == H5TypesCHAR) {
            if (buffer.size() < numAttributeElements) {
                buffer.resize(numAttributeElements);
            }

            READSTEPATTRIB(String, source, attributeName, &buffer[0]);
            WRITESTRINGSTEPATTRIB(file_m, attributeName, &buffer[0]);

        } else if (attributeType == H5TypesFLOAT) {
            if (buffer.size() < numAttributeElements * sizeof(h5_float32_t)) {
                buffer.resize(numAttributeElements * sizeof(h5_float32_t));
            }

            READSTEPATTRIB(Float32, source, attributeName, f32buffer);
            WRITESTEPATTRIB(Float32, file_m, attributeName, f32buffer, numAttributeElements);

        } else if (attributeType == H5TypesDOUBLE) {
            if (buffer.size() < numAttributeElements * sizeof(h5_float64_t)) {
                buffer.resize(numAttributeElements * sizeof(h5_float64_t));
            }

            READSTEPATTRIB(Float64, source, attributeName, f64buffer);
            WRITESTEPATTRIB(Float64, file_m, attributeName, f64buffer, numAttributeElements);

        } else if (attributeType == H5TypesINT32) {
            if (buffer.size() < numAttributeElements * sizeof(h5_int32_t)) {
                buffer.resize(numAttributeElements * sizeof(h5_int32_t));
            }

            READSTEPATTRIB(Int32, source, attributeName, i32buffer);
            WRITESTEPATTRIB(Int32, file_m, attributeName, i32buffer, numAttributeElements);

        } else if (attributeType == H5TypesINT64) {
            if (buffer.size() < numAttributeElements * sizeof(h5_int64_t)) {
                buffer.resize(numAttributeElements * sizeof(h5_int64_t));
            }

            READSTEPATTRIB(Int64, source, attributeName, i64buffer);
            WRITESTEPATTRIB(Int64, file_m, attributeName, i64buffer, numAttributeElements);

        } else {
            throw OpalException("H5PartWrapper::copyStepHeader",
                                "unknown data type: " + std::to_string(attributeType));
        }
    }
}

void H5PartWrapper::copyStepData(
    h5_file_t source
    ) {
    h5_size_t lengthSetName = 256;
    char setName[lengthSetName];
    h5_int64_t setType;
    h5_size_t numSetElements;

    h5_ssize_t numParticles = H5PartGetNumParticles(source);
    h5_ssize_t numParticlesPerNode = numParticles / Ippl::getNodes();

    h5_ssize_t firstParticle = numParticlesPerNode * Ippl::myNode();
    h5_ssize_t lastParticle = firstParticle + numParticlesPerNode - 1;
    if (Ippl::myNode() == Ippl::getNodes() - 1)
        lastParticle = numParticles - 1;

    REPORTONERROR(H5PartSetView(source, firstParticle, lastParticle));

    numParticles = lastParticle - firstParticle + 1;
    REPORTONERROR(H5PartSetNumParticles(file_m, numParticles));

    std::vector<char> buffer(numParticles * sizeof(h5_float64_t));
    h5_float32_t *f32buffer = reinterpret_cast<h5_float32_t*>(&buffer[0]);
    h5_float64_t *f64buffer = reinterpret_cast<h5_float64_t*>(&buffer[0]);
    h5_int32_t *i32buffer = reinterpret_cast<h5_int32_t*>(&buffer[0]);
    h5_int64_t *i64buffer = reinterpret_cast<h5_int64_t*>(&buffer[0]);

   h5_ssize_t numDatasets = H5PartGetNumDatasets(source);

    for (h5_ssize_t i = 0; i < numDatasets; ++ i) {
        REPORTONERROR(H5PartGetDatasetInfo(source, i, setName, lengthSetName, &setType, &numSetElements));

        if (setType == H5TypesFLOAT) {
            READDATA(Float32, source, setName, f32buffer);
            WRITEDATA(Float32, file_m, setName, f32buffer);
        } else if (setType == H5TypesDOUBLE) {
            READDATA(Float64, source, setName, f64buffer);
            WRITEDATA(Float64, file_m, setName, f64buffer);
        } else if (setType == H5TypesINT32) {
            READDATA(Int32, source, setName, i32buffer);
            WRITEDATA(Int32, file_m, setName, i32buffer);
        } else if (setType == H5TypesINT64) {
            READDATA(Int64, source, setName, i64buffer);
            WRITEDATA(Int64, file_m, setName, i64buffer);
        } else {
            throw OpalException("H5PartWrapper::copyStepData",
                                "unknown data type: " + std::to_string(setType));
        }
    }

    numParticles = H5PartGetNumParticles(file_m);
    REPORTONERROR(H5PartSetView(source, -1, -1));
}

void H5PartWrapper::sendFailureMessage(bool failed,
                                       const std::string &where,
                                       const std::string &what) {
    int tag = 101;
    Message *mess = new Message();
    putMessage(*mess, failed);
    Ippl::Comm->broadcast_all(mess, tag);
    delete mess;

    if (failed) throw OpalException(where, what);
}

void H5PartWrapper::receiveFailureMessage(int sourceNode,
                                          const std::string &where,
                                          const std::string &what) {
    int tag = 101;
    bool failed;
    Message *mess = Ippl::Comm->receive_block(sourceNode, tag);
    getMessage(*mess, failed);
    delete mess;

    if (failed) throw OpalException(where, what);
}

size_t H5PartWrapper::getNumParticles() const {
    if (!file_m) {
        throw OpalException("H5PartWrapper::getNumParticles",
                            "no file opened");
    }

    REPORTONERROR(H5SetStep(file_m, numSteps_m - 1));
    h5_ssize_t numParticles = H5PartGetNumParticles(file_m);

    return numParticles;
}

// vi: set et ts=4 sw=4 sts=4:
// Local Variables:
// mode:c++
// c-basic-offset: 4
// indent-tabs-mode: nil
// require-final-newline: nil
// End:
