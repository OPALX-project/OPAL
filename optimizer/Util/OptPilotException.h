//
// Class OptPilotException
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
#ifndef __OPTPILOTEXCEPTION_H__
#define __OPTPILOTEXCEPTION_H__

#include <string>

class OptPilotException {

public:

    OptPilotException(const std::string &meth, const std::string &descr) {
        descr_ = descr;
        meth_ = meth;
    }

    virtual const char* what() const throw() {
        return descr_.c_str();
    }

    virtual const char* where() const throw() {
        return meth_.c_str();
    }

private:

    std::string descr_;
    std::string meth_;

};

#endif