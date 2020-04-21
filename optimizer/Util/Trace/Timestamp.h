//
// Class Timestamp
//
// Copyright (c) 2010 - 2013, Yves Ineichen, ETH Zürich
// All rights reserved
//
// Implemented as part of the PhD thesis
// "Toward massively parallel multi-objective optimization with application to
// particle accelerators" (https://doi.org/10.3929/ethz-a-009792359)
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
#ifndef __TIMESTAMP_H__
#define __TIMESTAMP_H__

#include <sstream>
#include <boost/chrono.hpp>
#include <ctime>

#include "Util/Trace/TraceComponent.h"

#include "mpi.h"

class Timestamp : public TraceComponent {

public:


    Timestamp()
        : TraceComponent("Timestamp")
    {}

    void execute(std::ostringstream &dump) {

        boost::chrono::time_point<boost::chrono::system_clock> now;
        now = boost::chrono::system_clock::now();
        std::time_t now_time = boost::chrono::system_clock::to_time_t(now);

        std::ostringstream timestamp;
        timestamp << std::ctime(&now_time);

        prepend(dump, timestamp);
    }

    virtual ~Timestamp()
      {}

};

#endif