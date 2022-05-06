//
// Class StepSizeConfig
//
// This class stores tuples of time step sizes, path length range limits and limit of number of step sizes.
//
// Copyright (c) 2019 - 2021, Christof Metzger-Kraus
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
#ifndef STEPSIZECONFIG_H
#define STEPSIZECONFIG_H

#include "Structure/ValueRange.h"
#include "Utility/Inform.h"

#include <list>
#include <tuple>

class StepSizeConfig {
public:
    StepSizeConfig();

    StepSizeConfig(const StepSizeConfig &right);

    void operator=(const StepSizeConfig &) = delete;

    void push_back(double dt,
                   double zstop,
                   unsigned long numSteps);

    void sortAscendingZStop();

    void resetIterator();

    bool reachedStart() const;

    bool reachedEnd() const;

    void clear();

    void reverseDirection();

    StepSizeConfig& advanceToPos(double spos);

    StepSizeConfig& operator++();

    StepSizeConfig& operator--();

    void shiftZStopRight(double front);
    void shiftZStopLeft(double back);

    double getdT() const;

    double getZStop() const;

    unsigned long getNumSteps() const;

    unsigned long long getMaxSteps() const;

    unsigned long long getNumStepsFinestResolution() const;

    double getMinTimeStep() const;

    double getFinalZStop() const;

    void print(Inform &out) const;

    ValueRange<double> getPathLengthRange() const;

private:
    typedef std::tuple<double, double, unsigned long> entry_t;
    typedef std::list<entry_t> container_t;

    container_t configurations_m;
    container_t::iterator it_m;
};

inline
StepSizeConfig::StepSizeConfig():
    configurations_m(),
    it_m(configurations_m.begin())
{ }

inline
StepSizeConfig::StepSizeConfig(const StepSizeConfig &right):
    configurations_m(right.configurations_m),
    it_m(configurations_m.begin())
{ }

inline
void StepSizeConfig::push_back(double dt,
                               double zstop,
                               unsigned long numSteps) {
    configurations_m.push_back(std::make_tuple(dt, zstop, numSteps));
}

inline
void StepSizeConfig::resetIterator() {
    it_m = configurations_m.begin();
}

inline
bool StepSizeConfig::reachedStart() const {
    return (it_m == configurations_m.begin());
}

inline
bool StepSizeConfig::reachedEnd() const {
    return (it_m == configurations_m.end());
}

inline
void StepSizeConfig::clear() {
    configurations_m.clear();
    it_m = configurations_m.begin();
}

#endif