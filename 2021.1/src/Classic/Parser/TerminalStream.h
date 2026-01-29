//
// Class TerminalStream
//   A stream of input tokens.
//   The source of tokens is the terminal.
//
// Copyright (c) 2012-2019, Paul Scherrer Institut, Villigen PSI, Switzerland
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
#ifndef CLASSIC_TerminalStream_HH
#define CLASSIC_TerminalStream_HH

#include "Parser/AbsFileStream.h"

class TerminalStream: public AbsFileStream {

public:

    /// Constructor.
    TerminalStream(const char program[]);

    virtual ~TerminalStream();

    /// Read next input line.
    virtual bool fillLine();

private:

    // Not implemented.
    TerminalStream(const TerminalStream &);
    void operator=(const TerminalStream &);
};

#endif // CLASSIC_TerminalStream_HH
