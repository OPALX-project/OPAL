#ifndef OPAL_H5PARTWRAPPER_H
#define OPAL_H5PARTWRAPPER_H

//
//  Copyright & License: See Copyright.readme in src directory
//

/*!
  H5PartWrapper: a class that manages calls to H5Part
*/

template <class T, unsigned Dim>
class PartBunchBase;

#include "Utility/IpplInfo.h"

#include "H5hut.h"

#include <map>
#include <string>

#define REPORTONERROR(rc) H5PartWrapper::reportOnError(rc, __FILE__, __LINE__)
#define READFILEATTRIB(type, file, name, value) REPORTONERROR(H5ReadFileAttrib##type(file, name, value));
#define WRITESTRINGFILEATTRIB(file, name, value) REPORTONERROR(H5WriteFileAttribString(file, name, value));
#define WRITEFILEATTRIB(type, file, name, value, length) REPORTONERROR(H5WriteFileAttrib##type(file, name, value, length));

#define READSTEPATTRIB(type, file, name, value) REPORTONERROR(H5ReadStepAttrib##type(file, name, value));
#define WRITESTRINGSTEPATTRIB(file, name, value) REPORTONERROR(H5WriteStepAttribString(file, name, value));
#define WRITESTEPATTRIB(type, file, name, value, length) REPORTONERROR(H5WriteStepAttrib##type(file, name, value, length));

#define READDATA(type, file, name, value) REPORTONERROR(H5PartReadData##type(file, name, value));
#define WRITEDATA(type, file, name, value) REPORTONERROR(H5PartWriteData##type(file, name, value));

class H5PartWrapper {
public:
    virtual ~H5PartWrapper();

    void close();

    double getLastPosition();
    virtual void readHeader() = 0;
    virtual void readStep(PartBunchBase<double, 3>*, h5_ssize_t firstParticle, h5_ssize_t lastParticle) = 0;

    virtual void writeHeader() = 0;
    virtual void writeStep(PartBunchBase<double, 3>*, const std::map<std::string, double> &) = 0;

    virtual bool predecessorIsSameFlavour() const = 0;

    void storeCavityInformation();

    size_t getNumParticles() const;
protected:
    H5PartWrapper(const std::string &fileName, h5_int32_t flags = H5_O_WRONLY);
    H5PartWrapper(const std::string &fileName, int restartStep, std::string sourceFile, h5_int32_t flags = H5_O_RDWR);

    void open(h5_int32_t flags);

    void copyFile(const std::string &sourceFile, int lastStep = -1, h5_int32_t flags = H5_O_WRONLY);
    void copyFileSystem(const std::string &sourceFile);
    void copyHeader(h5_file_t source);
    void copyStep(h5_file_t source, int step);
    void copyStepHeader(h5_file_t source);
    void copyStepData(h5_file_t source);
    void sendFailureMessage(bool failed,
                            const std::string &where,
                            const std::string &what);
    void receiveFailureMessage(int sourceNode,
                               const std::string &where,
                               const std::string &what);

    static void reportOnError(h5_int64_t rc, const char* file, int line);

    h5_file_t file_m;
    std::string fileName_m;
    std::string predecessorOPALFlavour_m;
    h5_int64_t numSteps_m;
    bool startedFromExistingFile_m;

    static std::string copyFilePrefix_m;
};

inline
void H5PartWrapper::reportOnError(h5_int64_t rc, const char* file, int line) {
    if (rc != H5_SUCCESS)
        ERRORMSG("H5 rc= " << rc << " in " << file << " @ line " << line << endl);
}

inline
double H5PartWrapper::getLastPosition() {
    if (!file_m) open(H5_O_RDONLY);

    h5_ssize_t numStepsInSource = H5GetNumSteps(file_m);
    h5_ssize_t readStep = numStepsInSource - 1;
    REPORTONERROR(H5SetStep(file_m, readStep));

    h5_float64_t pathLength;
    READSTEPATTRIB(Float64, file_m, "SPOS", &pathLength);

    return pathLength;
}

#endif //OPAL_H5PARTWRAPPER_H
