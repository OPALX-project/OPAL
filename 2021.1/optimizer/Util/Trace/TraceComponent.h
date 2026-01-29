//
// Class TraceComponent
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
#ifndef __TRACE_COMPONENT_H__
#define __TRACE_COMPONENT_H__

#include <string>
#include <sstream>

class TraceComponent {

public:

    TraceComponent(std::string name) : name_(name)
    {}

    ~TraceComponent()
    {}

    virtual void execute(std::ostringstream &dump) = 0;

    void prepend(std::ostringstream &dump, std::ostringstream &prepender) {

        prepender << dump.str();
        dump.str("");
        dump.clear();
        dump << prepender.str();
    }

    void prepend(std::ostringstream &dump, std::string prepender) {

        std::ostringstream tmp;
        tmp << prepender << dump.str();
        dump.str("");
        dump.clear();
        dump << tmp.str();
    }

private:

    std::string name_;

};

#endif
