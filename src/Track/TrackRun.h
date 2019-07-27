#ifndef OPAL_TrackRun_HH
#define OPAL_TrackRun_HH

// ------------------------------------------------------------------------
// $RCSfile: TrackRun.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1.4.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: TrackRun
//
// ------------------------------------------------------------------------
//
// $Date: 2004/11/12 20:10:12 $
// $Author: adelmann $
//
// ------------------------------------------------------------------------

#include "AbstractObjects/Action.h"

#include <string>

class Beam;
class OpalData;
class DataSink;
class Distribution;
class Tracker;
class ParallelTTracker;
class FieldSolver;
class H5PartWrapper;

// Class TrackRun
// ------------------------------------------------------------------------
/// The RUN command.

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

    void setupSliceTracker();
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