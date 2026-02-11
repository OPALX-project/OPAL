//
// Class TrackRun
//   The RUN command.
//
// Copyright (c) 200x - 2023, Paul Scherrer Institut, Villigen PSI, Switzerland
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
#ifndef OPAL_TrackRun_HH
#define OPAL_TrackRun_HH

#include "AbstractObjects/Action.h"

#include <array>
#include <memory>
#include <string>
#include <string_view>
#include <vector>

class Beam;
class OpalData;
class DataSink;
class Distribution;
class FieldSolver;
class H5PartWrapper;
class Inform;
class ParallelTTracker;
class Tracker;

class TrackRun: public Action {

public:
    /// Exemplar constructor.
    TrackRun();

    virtual ~TrackRun();

    /// Make clone.
    virtual TrackRun* clone(const std::string& name);

    /// Execute the command.
    virtual void execute();

    using Action::print;
    Inform& print(Inform& os) const;

    static std::shared_ptr<Tracker> getTracker();

private:
    enum class RunMethod: unsigned short {
        NONE,
        PARALLELT,
        CYCLOTRONT,
        THICK
    };

    static inline constexpr std::array<std::pair<RunMethod, std::string_view>, 5> runMethodMap {{
        {TrackRun::RunMethod::PARALLELT,  "PARALLEL-T"},
        {TrackRun::RunMethod::PARALLELT,  "OPAL-T"},
        {TrackRun::RunMethod::CYCLOTRONT, "CYCLOTRON-T"},
        {TrackRun::RunMethod::CYCLOTRONT, "OPAL-CYCL"},
        {TrackRun::RunMethod::THICK,      "THICK"}
    }};

    // Not implemented.
    TrackRun(const TrackRun&);
    void operator=(const TrackRun&);

    // Clone constructor.
    TrackRun(const std::string& name, TrackRun* parent);

    void setRunMethod();
    std::string getRunMethodName() const;

    void setupTTracker();
    void setupCyclotronTracker();
    void setupThickTracker();
    void setupFieldsolver();

    void initPhaseSpaceSink();

    void initDataSink(const int& numBunch = 1);

    void setBoundaryGeometry();

    double setDistributionParallelT(Beam* beam);

    /*  itsTracker_m is a static object; this enables access to the last executed
     *  tracker object without excessive gymnastics, e.g. for access to the 
     *  field maps in PyField
     */
    static std::shared_ptr<Tracker> itsTracker_m;

    Distribution* dist_m;
    std::vector<Distribution*> distrs_m;

    FieldSolver* fieldSolver_m;

    DataSink* dataSink_m;

    H5PartWrapper* phaseSpaceSink_m;

    OpalData* opalData_m;

    bool isFollowupTrack_m;

    static const std::string defaultDistribution;

    RunMethod method_m;

    // macromass and charge for simulation particles
    double macromass_m;
    double macrocharge_m;
};

inline Inform& operator<<(Inform& os, const TrackRun& b) {
    return b.print(os);
}

#endif // OPAL_TrackRun_HH
