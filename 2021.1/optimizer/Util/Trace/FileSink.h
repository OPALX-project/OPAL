//
// Class FileSink
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
#ifndef __TRACE_FILE_SINK_H__
#define __TRACE_FILE_SINK_H__

#include <sstream>
#include <iostream>
#include <fstream>

#include "Util/Trace/TraceComponent.h"

class FileSink : public TraceComponent {

public:


    FileSink(std::string filename)
        : TraceComponent("FileSink")
        , filename_(filename)
    {}

    virtual ~FileSink()
    {}

    void execute(std::ostringstream &dump) {
        std::ofstream file;
        file.open(filename_.c_str(), std::ios::app);
        file << dump.str() << std::flush;
        file.close();
    }

private:

    std::string filename_;

};

#endif