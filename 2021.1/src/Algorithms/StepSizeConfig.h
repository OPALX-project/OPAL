#ifndef STEPSIZECONFIG_H
#define STEPSIZECONFIG_H

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