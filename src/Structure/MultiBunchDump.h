#ifndef OPAL_MULTI_BUNCH_DUMP_H
#define OPAL_MULTI_BUNCH_DUMP_H

#include <string>
#include <fstream>

class MultiBunchDump {

public:
    struct beaminfo_t {
        double time;
        double azimuth;
        unsigned int nParticles;
        double ekin;
        double dEkin;
        double rrms[3];
        double prms[3];
        double emit[3];
        double mean[3];
        double halo[3];
    };

public:
    MultiBunchDump();

    void writeHeader(const std::string& fname) const;

    void writeData(const beaminfo_t& binfo, short bunch);

private:
    void open_m(std::ofstream& out, const std::string& fname) const;

    void close_m(std::ofstream& out) const;

    // in non-restart mode we delete all *.smb files initially
    void remove_m() const;

private:
    std::string fbase_m;
    std::string fext_m;
};

#endif
