#ifndef OPAL_TrackCmd_HH
#define OPAL_TrackCmd_HH

// ------------------------------------------------------------------------
// $RCSfile: TrackCmd.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: TrackCmd
//
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:33:46 $
// $Author: Andreas Adelmann $
//
// ------------------------------------------------------------------------

#include "AbstractObjects/Action.h"


// Class TrackCmd
// ------------------------------------------------------------------------
/// The TRACK command.

class TrackCmd: public Action {

public:

    /// Exemplar constructor.
    TrackCmd();

    virtual ~TrackCmd();

    /// Make clone.
    virtual TrackCmd *clone(const std::string &name);

    /// Execute the command.
    virtual void execute();

    /// Return the timestep in seconds
    std::vector<double> getDT() const;

    double getDTSCINIT() const;

    double getDTAU() const;

    /// Return the elapsed time (sec) of the bunch
    double getT0() const;

    /// Return the maximum timsteps we integrate the system
    std::vector<unsigned long long> getMAXSTEPS() const;

    /// Return the timsteps per revolution period. ONLY available for OPAL-cycl.
    /// In OPAL-cycl, timestep is calculated by STEPSPERTURN, rather than given in TRACK command.
    int getSTEPSPERTURN() const;

    /// location at which the simulation starts
    double getZSTART() const;

    /// location at which the simulation stops
    std::vector<double> getZSTOP() const;

    /// return the name of time integrator
    int getTIMEINTEGRATOR() const;

private:

    // Not implemented.
    TrackCmd(const TrackCmd &);
    void operator=(const TrackCmd &);

    // Clone constructor.
    TrackCmd(const std::string &name, TrackCmd *parent);
};

#endif // OPAL_TrackCmd_HH