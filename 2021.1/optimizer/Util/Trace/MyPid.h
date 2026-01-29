//
// Class MyPid
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
#include "Util/Trace/TraceComponent.h"

#include "mpi.h"

class MyPid : public TraceComponent {

public:


    MyPid(std::string name, MPI_Comm comm)
        : TraceComponent(name)
    {
        mypid_ = 0;
        MPI_Comm_rank(comm, &mypid_);
    }

    void execute(std::ostringstream &dump) {
        dump << mypid_;
    }

private:

    int mypid_;

};
