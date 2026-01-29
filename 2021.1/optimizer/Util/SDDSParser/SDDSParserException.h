//
// Class SDDSParserException
//
// Copyright (c) 2018, Christof Metzger-Kraus, Open Sourcerer
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
#ifndef __SDDSPARSEREXCEPTION_H__
#define __SDDSPARSEREXCEPTION_H__

#include <string>

class SDDSParserException {

public:

    SDDSParserException(const std::string &meth, const std::string &descr) {
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