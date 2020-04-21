//
// Class Optimizer
//   An abstract class defining the interface for all optimizer
//   components.
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
#ifndef __OPTIMIZER_H__
#define __OPTIMIZER_H__

#include <vector>

#include "Util/Types.h"
#include "Pilot/Poller.h"

class Optimizer : protected Poller {

public:

    Optimizer(MPI_Comm comm) : Poller(comm) {}
    virtual ~Optimizer() {}

    /// type of bounds for design variables
    typedef std::vector< std::pair<double, double> > bounds_t;

    /// entry point for optimizer
    virtual void initialize() = 0;

protected:

    // propagate poller hooks
    virtual void setupPoll() = 0;
    virtual void prePoll() = 0;
    virtual void postPoll() = 0;
    virtual void onStop() = 0;
    virtual bool onMessage(MPI_Status status, size_t length) = 0;

private:

};

#endif
