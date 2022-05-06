//
// Class TrackCmd
//   The class for the OPAL TRACK command.
//
// Copyright (c) 200x - 2022, Paul Scherrer Institut, Villigen PSI, Switzerland
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
#ifndef OPAL_TrackCmd_HH
#define OPAL_TrackCmd_HH

#include "AbstractObjects/Action.h"
#include "Steppers/Steppers.h"

#include <map>


class TrackCmd: public Action {

public:
    TrackCmd();

    virtual ~TrackCmd();

    virtual TrackCmd* clone(const std::string& name);

    virtual void execute();

    /// Return the timestep in seconds
    std::vector<double> getDT() const;

    double getDTSCINIT() const;

    double getDTAU() const;

    /// Return the elapsed time (sec) of the bunch
    double getT0() const;

    /// Return the maximum timsteps we integrate the system
    std::vector<unsigned long long> getMaxSteps() const;

    /// Return the timsteps per revolution period. ONLY available for OPAL-cycl.
    /// In OPAL-cycl, timestep is calculated by STEPSPERTURN, rather than given in TRACK command.
    int getStepsPerTurn() const;

    /// location at which the simulation starts
    double getZStart() const;

    /// location at which the simulation stops
    std::vector<double> getZStop() const;

    /// return the name of time integrator
    Steppers::TimeIntegrator getTimeIntegrator();

private:
    // Not implemented.
    TrackCmd(const TrackCmd&);
    void operator=(const TrackCmd&);

    // Clone constructor.
    TrackCmd(const std::string& name, TrackCmd* parent);

    static const std::map<std::string, Steppers::TimeIntegrator> stringTimeIntegrator_s;
};

#endif // OPAL_TrackCmd_HH
