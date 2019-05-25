#ifndef OPAL_ORBITTHREADER_H
#define OPAL_ORBITTHREADER_H

#include "Algorithms/IndexMap.h"
#include "Algorithms/Vektor.h"
#include "Elements/OpalBeamline.h"
#include "Steppers/BorisPusher.h"

#include <string>
#include <map>

class OrbitThreader
{
public:

    OrbitThreader(const PartData &ref,
                  const Vector_t &r,
                  const Vector_t &p,
                  double s,
                  double maxDiffZBunch,
                  double t,
                  double dT,
                  size_t maxIntegSteps,
                  double zstop,
                  OpalBeamline &bl);

    void execute();

    IndexMap::value_t query(IndexMap::key_t::first_type step,
                            IndexMap::key_t::second_type length);

    std::pair<double, double> getRange(const IndexMap::value_t::value_type &element,
                                       double position) const;
    IndexMap::value_t getTouchingElements(const std::pair<double, double> &range);

private:
    Vector_t r_m;
    Vector_t p_m;
    double pathLength_m;
    double dZ_m;
    double time_m;
    double dt_m;
    const size_t maxIntegSteps_m;
    const double zstop_m;

    OpalBeamline &itsOpalBeamline_m;
    IndexMap imap_m;

    unsigned int errorFlag_m;

    BorisPusher integrator_m;
    const PartData &reference_m;

    std::ofstream logger_m;
    size_t loggingFrequency_m;

    struct elementPosition {
        double startField_m;
        double endField_m;
        double elementEdge_m;
    };

    struct elementPositionComp {
        bool operator() (const elementPosition &a, const elementPosition &b) const {
            return a.elementEdge_m < b.elementEdge_m;
        }
    };

    std::multimap<std::string, elementPosition> elementRegistry_m;

    void trackBack(double maxDrift);
    void integrate(const IndexMap::value_t &activeSet, size_t maxSteps, double maxDrift = 10.0);
    bool containsCavity(const IndexMap::value_t &activeSet);
    void autophaseCavities(const IndexMap::value_t &activeSet, const std::set<std::string> &visitedElements);
    double getMaxDesignEnergy(const IndexMap::value_t &elementSet) const;

    void registerElement(const IndexMap::value_t &elementSet, double, const Vector_t &r, const Vector_t &p);
    void processElementRegister();
    void setDesignEnergy(FieldList &allElements, const std::set<std::string> &visitedElements);
    double computeMaximalImplicitDrift();
};

inline
IndexMap::value_t OrbitThreader::query(IndexMap::key_t::first_type pathLength,
                                       IndexMap::key_t::second_type length) {
    return imap_m.query(pathLength, length);
}

inline
std::pair<double, double> OrbitThreader::getRange(const IndexMap::value_t::value_type &element,
                                                  double position) const {
    return imap_m.getRange(element, position);
}

inline
IndexMap::value_t OrbitThreader::getTouchingElements(const std::pair<double, double> &range) {
    return imap_m.getTouchingElements(range);
}

#endif