#ifndef OPAL_MULTI_BUNCH_HANDLER_H
#define OPAL_MULTI_BUNCH_HANDLER_H

#include "Algorithms/PartBunchBase.h"

#include <vector>

/* Helper class that stores bunch injection
 * information like azimuth, radius etc. of first
 * bunch in multi-bunch mode of ParallelCyclotronTracker.
 */
class MultiBunchHandler {
public:
    struct injection_t {
        injection_t()
            : time(0.0)
            , pathlength(0.0)
            , azimuth(0.0)
            , radius(0.0)
        { };

        double time;            // ns
        double pathlength;      // m
        double azimuth;         // deg
        double radius;          // mm
    };

    struct beaminfo_t {
        beaminfo_t(const injection_t& injection = injection_t())
            : time(injection.time)
            , azimuth(injection.azimuth)
            , radius(injection.radius)
            , prevAzimuth(-1.0)
            , pathlength(injection.pathlength)
            , nParticles(0)
            , ekin(0.0)
            , dEkin(0.0)
            , rrms{0.0}
            , prms{0.0}
            , emit{0.0}
            , mean{0.0}
            , correlation{0.0}
            , halo{0.0}
        { };

        double time;
        double azimuth;
        double radius;
        double prevAzimuth;
        double pathlength;
        long unsigned int nParticles;
        double ekin;
        double dEkin;
        double rrms[3];
        double prms[3];
        double emit[3];
        double mean[3];
        double correlation[3];
        double halo[3];
    };

    // multi-bunch modes
    enum class MB_MODE {
        FORCE  = 0,
        AUTO   = 1
    };

    // multi-bunch binning type
    enum class MB_BINNING {
        GAMMA = 0,
        BUNCH = 1
    };

    /*
     * @param numBunch > 1 --> multi bunch mode
     * @param eta binning value
     * @param para only for MB_MODE::AUTO
     * @param mode of multi-bunch
     * @param binning type of particle binning
     */
    MultiBunchHandler(PartBunchBase<double, 3> *beam,
                      const int& numBunch,
                      const double& eta,
                      const double& para,
                      const std::string& mode,
                      const std::string& binning);

    void saveBunch(PartBunchBase<double, 3> *beam);

    bool readBunch(PartBunchBase<double, 3> *beam,
                   const PartData& ref);

    /* Returns:
     *     0 - if nothing happened
     *     1 - if bunch got saved
     *     2 - if bunch got injected
     */
    short injectBunch(PartBunchBase<double, 3> *beam,
                      const PartData& ref,
                      bool& flagTransition);

    void updateParticleBins(PartBunchBase<double, 3> *beam);
    
    /// set the working sub-mode for multi-bunch mode: "FORCE" or "AUTO"
    void setMode(const std::string& mbmode);

    // set binning type
    void setBinning(std::string binning);

    void setRadiusTurns(const double& radius);

    void setNumBunch(short n);

    short getNumBunch() const;

    bool isForceMode() const;

    bool calcBunchBeamParameters(PartBunchBase<double, 3>* beam, short bunchNr);

    beaminfo_t& getBunchInfo(short bunchNr);

    const beaminfo_t& getBunchInfo(short bunchNr) const;

    injection_t& getInjectionValues();

    void updateTime(const double& dt);

    void updatePathLength(const std::vector<double>& lpaths);

private:
    // store the data of the beam which are required for injecting a
    // new bunch for multibunch filename
    std::string onebunch_m;

    /// The number of bunches specified in TURNS of RUN commond
    short numBunch_m;

    // parameter for reset bin in multi-bunch run
    double eta_m;

    // 0 for single bunch (default),
    // 1 for FORCE,
    // 2 for AUTO
    MB_MODE mode_m; //multiBunchMode_m;

    // 0 for GAMMA (default),
    // 1 for BUNCH
    MB_BINNING binning_m; //binningType_m;

    // control parameter for AUTO multi-bunch mode
    double coeffDBunches_m;

    // used for automatic injection in multi-bunch mode
    double radiusLastTurn_m;
    double radiusThisTurn_m;

    // record how many bunches has already been injected.
    short bunchCount_m;

    // each list entry belongs to a bunch
    std::vector<beaminfo_t> binfo_m;

    // global attributes of injection
    injection_t injection_m;
};


inline
void MultiBunchHandler::setNumBunch(short n) {
    bunchCount_m = n;
}


inline
short MultiBunchHandler::getNumBunch() const {
    return bunchCount_m;
}


inline
bool MultiBunchHandler::isForceMode() const {
    return (mode_m == MB_MODE::FORCE);
}


inline
MultiBunchHandler::beaminfo_t& MultiBunchHandler::getBunchInfo(short bunchNr) {
    PAssert_GE(bunchNr, 0);
    PAssert_LT(bunchNr, (short)binfo_m.size());
    return binfo_m[bunchNr];
}


inline
const MultiBunchHandler::beaminfo_t& MultiBunchHandler::getBunchInfo(short bunchNr) const {
    PAssert_GE(bunchNr, 0);
    PAssert_LT(bunchNr, (short)binfo_m.size());
    return binfo_m[bunchNr];
}


inline
MultiBunchHandler::injection_t& MultiBunchHandler::getInjectionValues() {
    return injection_m;
}

#endif
