#ifndef IPPLINFOWRAPPER_H
#define IPPLINFOWRAPPER_H

//
//  Copyright & License: See Copyright.readme in src directory
//

/*!
  Class documentation
*/

#include "Utility/IpplInfo.h"

class IpplInfoWrapper {
public:
    IpplInfoWrapper(const std::string &inputFileName, int infoLevel, int warnLevel, MPI_Comm comm);
    ~IpplInfoWrapper();

private:
    unsigned int exeName_m;
    unsigned int inputFileName_m;
    unsigned int noComm_m;
    unsigned int info_m;
    unsigned int infoLevel_m;
    unsigned int warn_m;
    unsigned int warnLevel_m;

    char *buffer_m;
    char **arg_m;

    Ippl *instance_m;
};

#endif