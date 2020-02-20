//
//  Copyright & License: See Copyright.readme in src directory
//

#ifndef OPAL_LOSSOUTPUT_H_
#define OPAL_LOSSOUTPUT_H_

//////////////////////////////////////////////////////////////
#include "Utility/IpplInfo.h"
#include "Algorithms/Vektor.h"
#include "AbsBeamline/ElementBase.h"
#include "AbstractObjects/OpalData.h"

#include <string>
#include <fstream>
#include <vector>
#include <functional>
#include <set>

#include <hdf5.h>
#include "H5hut.h"

struct SetStatistics {
    SetStatistics();

    std::string element_m;
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
    struct less<SetStatistics> : binary_function<SetStatistics, SetStatistics, bool> {
        bool operator() (const SetStatistics& x, const SetStatistics& y) const {
            return x.spos_m < y.spos_m;
        }
    };
}
/*
  - In the destructor we do ALL the file handling
  - h5hut_mode_m defines h5hut or ASCII
 */
class LossDataSink {
 public:
    LossDataSink();

    LossDataSink(std::string elem, bool hdf5Save, ElementBase::ElementType type = ElementBase::ANY);

    LossDataSink(const LossDataSink &rsh);
    ~LossDataSink() noexcept(false);

    bool inH5Mode() { return h5hut_mode_m;}

    void save(unsigned int numSets = 1, OpalData::OPENMODE openMode = OpalData::OPENMODE::UNDEFINED);

    void addReferenceParticle(const Vector_t &x,
                              const Vector_t &p,
                              double time,
                              double spos,
                              long long globalTrackStep);

    void addParticle(const Vector_t &x, const Vector_t &p, const size_t id);

    void addParticle(const Vector_t &x, const Vector_t &p, const size_t  id,
                     const double time, const size_t turn, const size_t& bunchNum = 0);

    size_t size() const;

    std::set<SetStatistics> computeStatistics(unsigned int numSets);

private:
    void openASCII() {
        if(Ippl::myNode() == 0) {
            os_m.open(fn_m.c_str(), std::ios::out);
        }
    }
    void openH5(h5_int32_t mode = H5_O_WRONLY);

    void appendASCII() {
        if(Ippl::myNode() == 0) {
            os_m.open(fn_m.c_str(), std::ios::app);
        }
    }

    void writeHeaderASCII() {
        if(Ippl::myNode() == 0) {
            //FIXME Issue #45 (Cyclotron units)
            os_m << "# Element " << element_m << " x (m),  y (m),  z (m),  px ( ),  py ( ),  pz ( ), id";
            if (time_m.size() != 0) {
                os_m << ",  turn, bunchNumber, time (ns) ";
            }
            os_m << std::endl;
        }
    }
    void writeHeaderH5();

    void saveASCII();
    void saveH5(unsigned int setIdx);

    void closeASCII() {
        if(Ippl::myNode() == 0)
            os_m.close();
    }

    bool hasNoParticlesToDump();

    bool hasTimeAttribute();

    void reportOnError(h5_int64_t rc, const char* file, int line);

    void splitSets(unsigned int numSets);
    SetStatistics computeSetStatistics(unsigned int setIdx);

    // filename without extension
    std::string fn_m;

    // write either in ASCII or H5hut format
    bool h5hut_mode_m;

    // used to write out data in ASCII mode
    std::ofstream os_m;

    /// used to write out data in H5hut mode
    h5_file_t H5file_m;

    std::string element_m;

    /// Current record, or time step, of H5 file.
    h5_int64_t H5call_m;

    std::vector<long>   id_m;
    std::vector<double>  x_m;
    std::vector<double>  y_m;
    std::vector<double>  z_m;
    std::vector<double> px_m;
    std::vector<double> py_m;
    std::vector<double> pz_m;
    std::vector<size_t> bunchNum_m;

    std::vector<size_t> turn_m;
    std::vector<double> time_m;

    std::vector<Vector_t> RefPartR_m;
    std::vector<Vector_t> RefPartP_m;
    std::vector<h5_int64_t> globalTrackStep_m;
    std::vector<double> refTime_m;
    std::vector<double> spos_m;

    std::vector<unsigned long> startSet_m;

    ElementBase::ElementType type_m;
};

inline
size_t LossDataSink::size() const {
    return x_m.size();
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

// vi: set et ts=4 sw=4 sts=4:
// Local Variables:
// mode:c++
// c-basic-offset: 4
// indent-tabs-mode: nil
// require-final-newline: nil
// End:
