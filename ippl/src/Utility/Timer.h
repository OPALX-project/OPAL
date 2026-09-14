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
#ifndef TIMER_H
#define TIMER_H

#include <chrono>
#include <ctime>

class Timer {
public:
    Timer();

    void clear();               // Set all accumulated times to 0
    void start();               // Start timer
    void stop();                // Stop timer

    double clock_time();        // Report clock time accumulated in seconds
    double cpu_time();          // Report total cpu_time which is just user_time + system_time

private:
    std::chrono::steady_clock::time_point wall_start_m;
    std::clock_t cpu_start_m;

    double wall_m;
    double cpu_m;
};

#endif