//
// Class Trace
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
#ifndef __TRACE_H__
#define __TRACE_H__

#include <string>
#include <sstream>
#include <vector>
#include <map>

#include "boost/smart_ptr.hpp"

#include "Util/Trace/TraceComponent.h"

class Trace {

public:

    Trace(std::string name)
        : name_(name)
    {}

    ~Trace()
    {}

    void registerComponent(std::string name,
            boost::shared_ptr<TraceComponent> component) {
        nameToIdx_.insert(
            std::pair<std::string, size_t>(name, pipeline_.size()));
        pipeline_.push_back(component);
    }

    void unregisterComponent(std::string /*name*/) {
        //TODO: set null @ idx
    }

    void log(std::ostringstream &dump) {
        for(boost::shared_ptr<TraceComponent> component : pipeline_) {
            component->execute(dump);
        }
    }

private:

    std::string name_;

    std::vector< boost::shared_ptr<TraceComponent> > pipeline_;
    std::map< std::string, size_t > nameToIdx_;

};

#endif
