//
// Class CoutSink
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
#ifndef __TRACE_COUT_SINK_H__
#define __TRACE_COUT_SINK_H__

#include <iostream>
#include <string>

#include "Util/Trace/TraceComponent.h"

class CoutSink : public TraceComponent {

public:


    CoutSink(std::string prefix = "")
        : TraceComponent("CoutSink")
        , prefix_(prefix) {

        clear_color_ = "\e[0m";
    }

    ~CoutSink()
    {}


    void setColor(std::string color)      { color_ = color; }
    void setClearColor(std::string color) { clear_color_ = color; }


    void execute(std::ostringstream &dump) {
        std::cout << color_ << prefix_
                  << dump.str()
                  << clear_color_ << std::flush;
    }

private:

    std::string prefix_;
    std::string color_;
    std::string clear_color_;

};

#endif
