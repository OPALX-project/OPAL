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
#include "Parser/TerminalStream.h"

#include <iomanip>
#include <iostream>

TerminalStream::TerminalStream(
    const char[]
    ): AbsFileStream("standard input") {
}


TerminalStream::~TerminalStream() {
}


bool TerminalStream::fillLine() {
    if(std::cin.eof()) {
        // We must test for end of file, even on terminal, as it may be rerouted.
        return false;
    } else {
        std::cerr << "==>";
        std::getline(std::cin, line, '\n');
        line += "\n";
        curr_line++;
        std::cerr.width(5);
        std::cerr << curr_line << " " << line;
        curr_char = 0;
        return true;
    }
}