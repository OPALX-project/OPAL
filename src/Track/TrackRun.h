//
// Class TrackRun
//   The RUN command.
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
#ifndef OPAL_TrackRun_HH
#define OPAL_TrackRun_HH

#include "AbstractObjects/Action.h"

#include <string>
#include <vector>

class Beam;
class OpalData;
class DataSink;
class Distribution;
class Tracker;
class ParallelTTracker;
class FieldSolver;
class H5PartWrapper;

class TrackRun: public Action {

public:

    /// Exemplar constructor.
    TrackRun();

    virtual ~TrackRun();

    /// Make clone.
    virtual TrackRun *clone(const std::string &name);

    /// Execute the command.
    virtual void execute();

private:

    // Not implemented.
    TrackRun(const TrackRun &);
    void operator=(const TrackRun &);

    // Clone constructor.
    TrackRun(const std::string &name, TrackRun *parent);

    void setupTTracker();
    void setupCyclotronTracker();
    void setupThickTracker();
    void setupFieldsolver();

    void initDataSink(const int& numBunch = 1);

    double setDistributionParallelT(Beam *beam);

    // Pointer to tracking algorithm.
    Tracker *itsTracker;

    Distribution *dist;

    std::vector<Distribution *> distrs_m;

    FieldSolver  *fs;

    DataSink *ds;

    H5PartWrapper *phaseSpaceSink_m;

    OpalData *opal;

    static const std::string defaultDistribution;
};

#endif // OPAL_TrackRun_HH