//
//  Copyright & License: See Copyright.readme in src directory
//

#include "Structure/IpplInfoWrapper.h"

#include <cstring>

IpplInfoWrapper::IpplInfoWrapper(const std::string &inputFileName, int infoLevel, int warnLevel, MPI_Comm comm) {
    std::string infoLevelStr = std::to_string(infoLevel);
    std::string warnLevelStr = std::to_string(warnLevel);

    exeName_m = 0;
    inputFileName_m = exeName_m + 5;
    noComm_m = inputFileName_m + inputFileName.size() + 1;
    info_m = noComm_m + 13;
    infoLevel_m = info_m + 7;
    warn_m = infoLevel_m + infoLevelStr.size() + 1;
    warnLevel_m = warn_m + 7;

    unsigned int totalSize = warnLevel_m + warnLevelStr.size() + 1;
    buffer_m = new char[totalSize];

    strcpy(buffer_m + exeName_m, "opal");
    strcpy(buffer_m + inputFileName_m, inputFileName.c_str());
    strcpy(buffer_m + noComm_m, "--nocomminit");
    strcpy(buffer_m + info_m, "--info");
    strcpy(buffer_m + infoLevel_m, infoLevelStr.c_str());
    strcpy(buffer_m + warn_m, "--warn");
    strcpy(buffer_m + warnLevel_m, warnLevelStr.c_str());

    arg_m = new char*[7];
    arg_m[0] = buffer_m + exeName_m;
    arg_m[1] = buffer_m + inputFileName_m;
    arg_m[2] = buffer_m + noComm_m;
    arg_m[3] = buffer_m + info_m;
    arg_m[4] = buffer_m + infoLevel_m;
    arg_m[5] = buffer_m + warn_m;
    arg_m[6] = buffer_m + warnLevel_m;

    int narg = 5;
    instance_m = new Ippl(narg, arg_m, Ippl::KEEP, comm);
}

IpplInfoWrapper::~IpplInfoWrapper() {
    delete instance_m;
    delete[] buffer_m;
    delete[] arg_m;

    Ippl::deleteGlobals();
}