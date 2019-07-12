#ifndef OPAL_MULTI_BUNCH_DUMP_H
#define OPAL_MULTI_BUNCH_DUMP_H

#include "SDDSWriter.h"

class MultiBunchDump : public SDDSWriter {

public:
    struct beaminfo_t {
        double time;
        double azimuth;
        long unsigned int nParticles;
        double ekin;
        double dEkin;
        double rrms[3];
        double prms[3];
        double emit[3];
        double mean[3];
        double halo[3];
    };

public:
    MultiBunchDump(const std::string& fname, bool restart,
                   const short& bunch);

    void fillHeader();

    using SDDSWriter::write;
    void write(PartBunchBase<double, 3>* beam, const double& azimuth);

    bool calcBunchBeamParameters(PartBunchBase<double, 3>* beam);

private:
    bool isFirst_m;

    // the bunch number
    short bunch_m;

    beaminfo_t binfo_m;
};

#endif
