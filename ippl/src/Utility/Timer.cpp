//
// Class Timer
//   This class is used in IpplTimings.
//
// Copyright (c) 2019, Matthias Frey, Paul Scherrer Institut, Villigen PSI, Switzerland
// All rights reserved
//
// Implemented as part of the PhD thesis
// "Precise Simulations of Multibunches in High Intensity Cyclotrons"
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
#include "Timer.h"

Timer::Timer() {
    this->clear();
}

void Timer::clear() {
    wall_m = 0.0;
    cpu_m = 0.0;
}

void Timer::start() {
    wall_start_m = std::chrono::steady_clock::now();
    cpu_start_m = std::clock();
}

void Timer::stop() {
    const auto wall_end = std::chrono::steady_clock::now();
    const auto cpu_end = std::clock();

    wall_m += std::chrono::duration<double>(
        wall_end - wall_start_m).count();

    cpu_m += static_cast<double>(cpu_end - cpu_start_m) /
             CLOCKS_PER_SEC;
}

double Timer::clock_time() {
    return wall_m;
}

double Timer::cpu_time() {
    return cpu_m;
}