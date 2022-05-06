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
#ifndef LOSSDATASINK_H_
#define LOSSDATASINK_H_

//////////////////////////////////////////////////////////////
#include "Algorithms/Vektor.h"
#include "Algorithms/OpalParticle.h"
#include "AbsBeamline/ElementBase.h"
#include "AbstractObjects/OpalData.h"

#include <boost/optional.hpp>

#include <string>
#include <fstream>
#include <vector>
#include <functional>
#include <set>

#include "H5hut.h"

struct SetStatistics {
    SetStatistics();

    std::string outputName_m;
    double spos_m;
    double refTime_m; // ns
    double tmean_m; // ns
    double trms_m; // ns
    unsigned long nTotal_m;
    Vector_t RefPartR_m;
    Vector_t RefPartP_m;
    Vector_t rmin_m;
    Vector_t rmax_m;
    Vector_t rmean_m;
    Vector_t pmean_m;
    Vector_t rrms_m;
    Vector_t prms_m;
    Vector_t rprms_m;
    Vector_t normEmit_m;
    Vector_t rsqsum_m;
    Vector_t psqsum_m;
    Vector_t rpsum_m;
    Vector_t eps2_m;
    Vector_t eps_norm_m;
    Vector_t fac_m;
};

namespace std {
    template<>
    struct less<SetStatistics> {
        bool operator() (const SetStatistics& x, const SetStatistics& y) const {
            return x.spos_m < y.spos_m;
        }
    };
}
enum class CollectionType: unsigned short {
                                           SPATIAL = 0,
                                           TEMPORAL
};

/*
  - In the destructor we do ALL the file handling
  - h5hut_mode_m defines h5hut or ASCII
 */
class LossDataSink {
 public:

    LossDataSink() = default;

    LossDataSink(std::string outfn, bool hdf5Save, CollectionType = CollectionType::TEMPORAL);

    LossDataSink(const LossDataSink &rsh);
    ~LossDataSink() noexcept(false);

    bool inH5Mode() { return h5hut_mode_m;}

    void save(unsigned int numSets = 1, OpalData::OpenMode openMode = OpalData::OpenMode::UNDEFINED);

    void addReferenceParticle(const Vector_t &x,
                              const Vector_t &p,
                              double time,
                              double spos,
                              long long globalTrackStep);

    void addParticle(const OpalParticle &, const boost::optional<std::pair<int, short int>> &turnBunchNumPair = boost::none);

    size_t size() const;

    std::set<SetStatistics> computeStatistics(unsigned int numSets);

private:
    void openASCII() {
        if(Ippl::myNode() == 0) {
            os_m.open(fileName_m.c_str(), std::ios::out);
        }
    }
    void openH5(h5_int32_t mode = H5_O_WRONLY);

    void appendASCII() {
        if (Ippl::myNode() == 0) {
            os_m.open(fileName_m.c_str(), std::ios::app);
        }
    }

    void writeHeaderASCII();
    void writeHeaderH5();

    void saveASCII();
    void saveH5(unsigned int setIdx);

    void closeASCII() {
        if (Ippl::myNode() == 0) {
            os_m.close();
        }
    }

    bool hasNoParticlesToDump() const;
    bool hasTurnInformations() const;

    void reportOnError(h5_int64_t rc, const char* file, int line);

    void splitSets(unsigned int numSets);
    SetStatistics computeSetStatistics(unsigned int setIdx);

    // filename without extension
    std::string fileName_m;

    // write either in ASCII or H5hut format
    bool h5hut_mode_m;

    // used to write out data in ASCII mode
    std::ofstream os_m;

    /// used to write out data in H5hut mode
    h5_file_t H5file_m;

    std::string outputName_m;

    /// Current record, or time step, of H5 file.
    h5_int64_t H5call_m;

    std::vector<OpalParticle> particles_m;
    std::vector<size_t> bunchNumber_m;
    std::vector<size_t> turnNumber_m;

    std::vector<Vector_t> RefPartR_m;
    std::vector<Vector_t> RefPartP_m;
    std::vector<h5_int64_t> globalTrackStep_m;
    std::vector<double> refTime_m;
    std::vector<double> spos_m;

    std::vector<unsigned long> startSet_m;

    CollectionType collectionType_m;
};

inline
size_t LossDataSink::size() const {
    return particles_m.size();
}

inline
std::set<SetStatistics> LossDataSink::computeStatistics(unsigned int numStatistics) {
    std::set<SetStatistics> stats;

    splitSets(numStatistics);

    for (unsigned int i = 0; i < numStatistics; ++ i) {
        auto setStats = computeSetStatistics(i);
        if (setStats.nTotal_m > 0) {
            stats.insert(setStats);
        }
    }

    return stats;
}

#endif
