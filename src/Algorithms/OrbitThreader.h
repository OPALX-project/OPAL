//
// Class OrbitThreader
//
// This class determines the design path by tracking the reference particle through
// the 3D lattice.
//
// Copyright (c) 2016,       Christof Metzger-Kraus, Helmholtz-Zentrum Berlin, Germany
//               2017 - 2020 Christof Metzger-Kraus
//
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
#ifndef OPAL_ORBITTHREADER_H
#define OPAL_ORBITTHREADER_H

#include "Algorithms/IndexMap.h"
#include "Algorithms/StepSizeConfig.h"
#include "Algorithms/Vektor.h"
#include "Elements/OpalBeamline.h"
#include "Steppers/BorisPusher.h"
#include "Structure/BoundingBox.h"
#include "Structure/ValueRange.h"

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
                  StepSizeConfig stepSizes,
                  OpalBeamline &bl);

    void execute();

    IndexMap::value_t query(IndexMap::key_t::first_type step,
                            IndexMap::key_t::second_type length);

    IndexMap::key_t getRange(const IndexMap::value_t::value_type &element,
                             double position) const;
    IndexMap::value_t getTouchingElements(const IndexMap::key_t &range) const;

    BoundingBox getBoundingBox() const;

private:
    /// position of reference particle in lab coordinates
    Vector_t r_m;
    /// momentum of reference particle
    Vector_t p_m;
    /// position of reference particle in path length
    double pathLength_m;
    /// distance to track back before tracking forward
    /// (length of bunch but not beyond cathode)
    double distTrackBack_m;
    /// the simulated time
    double time_m;
    /// the time step
    double dt_m;
    ValueRange<long> stepRange_m;
    long currentStep_m{0};

    /// final position in path length
    StepSizeConfig stepSizes_m;
    const double zstop_m;
    ValueRange<double> pathLengthRange_m;

    OpalBeamline &itsOpalBeamline_m;
    IndexMap imap_m;

    unsigned int errorFlag_m;

    BorisPusher integrator_m;
    const PartData &reference_m;

    std::ofstream logger_m;
    size_t loggingFrequency_m;

    BoundingBox globalBoundingBox_m;

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

    void trackBack();
    void integrate(const IndexMap::value_t &activeSet, double maxDrift = 10.0);
    bool containsCavity(const IndexMap::value_t &activeSet);
    void autophaseCavities(const IndexMap::value_t &activeSet, const std::set<std::string> &visitedElements);
    double getMaxDesignEnergy(const IndexMap::value_t &elementSet) const;

    void registerElement(const IndexMap::value_t &elementSet, double, const Vector_t &r, const Vector_t &p);
    void processElementRegister();
    void setDesignEnergy(FieldList &allElements, const std::set<std::string> &visitedElements);
    void computeBoundingBox();
    void updateBoundingBoxWithCurrentPosition();
    double computeDriftLengthToBoundingBox(const std::set<std::shared_ptr<Component>> & elements,
                                           const Vector_t & position,
                                           const Vector_t & direction) const;

    void checkElementLengths(const std::set<std::shared_ptr<Component>>& elements);
};

inline
IndexMap::value_t OrbitThreader::query(IndexMap::key_t::first_type pathLength,
                                       IndexMap::key_t::second_type length) {
    return imap_m.query(pathLength, length);
}

inline
IndexMap::key_t OrbitThreader::getRange(const IndexMap::value_t::value_type &element,
                                        double position) const {
    return imap_m.getRange(element, position);
}

inline
IndexMap::value_t OrbitThreader::getTouchingElements(const IndexMap::key_t &range) const {
    return imap_m.getTouchingElements(range);
}

inline
BoundingBox OrbitThreader::getBoundingBox() const {
    return globalBoundingBox_m;
}
#endif