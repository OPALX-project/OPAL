#include "Algorithms/StepSizeConfig.h"
#include "Utilities/OpalException.h"

#include <algorithm>
#include <numeric>
#include <iterator>
#include <cmath>

void StepSizeConfig::sortAscendingZStop() {
    configurations_m.sort([] (const entry_t &a,
                              const entry_t &b) -> bool
                          {
                              return std::get<1>(a) < std::get<1>(b);
                          });
}

void StepSizeConfig::reverseDirection() {
    unsigned int posIterator = std::distance(it_m, configurations_m.end()) - 1;
    configurations_m.reverse();
    it_m = configurations_m.begin();
    std::advance(it_m, posIterator);
}

StepSizeConfig& StepSizeConfig::advanceToPos(double spos) {
    while (getZStop() < spos && std::next(it_m) != configurations_m.end()) {
        ++ it_m;
    }

    return *this;
}

StepSizeConfig& StepSizeConfig::operator++() {
    if (reachedEnd()) {
        throw OpalException("StepSizeConfig::operator++",
                            "iterator is at end of list of configurations");
    }

    ++ it_m;

    return *this;
}

StepSizeConfig& StepSizeConfig::operator--() {
    if (reachedStart()) {
        throw OpalException("StepSizeConfig::operator--",
                            "iterator is at begin of list of configurations");
    }

    -- it_m;

    return *this;
}

void StepSizeConfig::shiftZStopRight(double front) {
    auto it = configurations_m.begin();
    while (std::get<1>(*it) < front &&
           std::next(it) != configurations_m.end()) {
        ++ it;
    }

    double zstop = std::get<1>(*it);
    if (zstop < front) return;

    std::get<1>(*it) = front;
    for (++ it ;
         it != configurations_m.end();
         ++ it) {
        std::swap(zstop, std::get<1>(*it));
    }
}

void StepSizeConfig::shiftZStopLeft(double back) {
    auto it = configurations_m.rbegin();
    while (std::get<1>(*it) > back &&
           std::next(it) != configurations_m.rend()) {
        ++ it;
    }

    double zstop = std::get<1>(*it);
    if (zstop > back) return;

    std::get<1>(*it) = back;
    for (++ it ;
         it != configurations_m.rend();
         ++ it) {
        std::swap(zstop, std::get<1>(*it));
    }
}

double StepSizeConfig::getdT() const {
    if (reachedEnd()) {
        throw OpalException("StepSizeConfig::getdT",
                            "iterator is at end of list of configurations");
    }

    return std::get<0>(*it_m);
}

double StepSizeConfig::getZStop() const {
    if (reachedEnd()) {
        throw OpalException("StepSizeConfig::getZStop",
                            "iterator is at end of list of configurations");
    }

    return std::get<1>(*it_m);
}

unsigned long StepSizeConfig::getNumSteps() const {
    if (reachedEnd()) {
        throw OpalException("StepSizeConfig::getNumSteps",
                            "iterator is at end of list of configurations");
    }

    return std::get<2>(*it_m);
}

unsigned long long StepSizeConfig::getMaxSteps() const {
    unsigned long long maxSteps = 0;
    for (const auto config: configurations_m) {
        maxSteps += std::get<2>(config);
    }

    return maxSteps;
}

unsigned long long StepSizeConfig::getNumStepsFinestResolution() const {
    double minTimeStep = std::get<0>(configurations_m.front());
    unsigned long long totalNumSteps = 0;

    for (const auto config: configurations_m) {
        const double &dt = std::get<0>(config);
        const unsigned long &numSteps = std::get<2>(config);

        if (minTimeStep > dt) {
            totalNumSteps = std::ceil(totalNumSteps * minTimeStep / dt);
            minTimeStep = dt;
        }

        totalNumSteps += std::ceil(numSteps * dt / minTimeStep);
    }

    return totalNumSteps;
}

double StepSizeConfig::getMinTimeStep() const {
    double minTimeStep = std::get<0>(configurations_m.front());
    for (const auto config: configurations_m) {
        if (minTimeStep > std::get<0>(config)) {
            minTimeStep = std::get<0>(config);
        }
    }

    return minTimeStep;
}

double StepSizeConfig::getFinalZStop() const {
    return std::get<1>(configurations_m.back());
}

void StepSizeConfig::print(Inform &out) const {
    out << std::scientific << "   "
        << std::setw(20) << "dt [ns] "
        << std::setw(20) << "zStop [m] "
        << std::setw(20) << "num Steps [1]"
        << endl;

    for (auto it = configurations_m.begin();
         it != configurations_m.end();
         ++ it) {
        if (it_m == it) {
            out << "-> ";
        } else {
            out << "   ";
        }

        out << std::setw(20) << std::get<0>(*it)
            << std::setw(20) << std::get<1>(*it)
            << std::setw(20) << std::get<2>(*it)
            << endl;
    }
}