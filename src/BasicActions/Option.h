//
// Class Option
//   The OPTION command.
//   The user interface allowing setting of OPAL options.
//   The actual option flags are contained in namespace Options.
//
// Copyright (c) 200x - 2020, Paul Scherrer Institut, Villigen PSI, Switzerland
// All rights reserved
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
#ifndef OPAL_Option_HH
#define OPAL_Option_HH

#include "AbstractObjects/Action.h"
#include "Utilities/Options.h"

#include <boost/bimap.hpp>

#include <string>

class Option: public Action {

public:

    Option();
    virtual ~Option();

    /// Make clone.
    virtual Option* clone(const std::string& name);

    /// Execute the command.
    virtual void execute();

private:

    void handlePsDumpFrame(const std::string& dumpFrame);

    static std::string getDumpFrameString(const DumpFrame& df);

    using Object::update;
    void update(const std::vector<Attribute>&);

    // Not implemented.
    Option(const Option&);
    void operator=(const Option&);

    // Clone constructor.
    Option(const std::string& name, Option* parent);

    static const boost::bimap<DumpFrame, std::string> bmDumpFrameString_s;
};

#endif // OPAL_Option_HH
