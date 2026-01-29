//
// Class Format
// Format class to allow serializing message objects into plain buffers
// to send directly with mpi calls or similar means
//
// Copyright (c) 2008 - 2020, Paul Scherrer Institut, Villigen PSI, Switzerland
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
#ifndef FORMATTER_H
#define FORMATTER_H

#include <algorithm>
#include <cstring>
#include <vector>
#include "Message/Message.h"

class Format {
public:
    Format(Message*);
    unsigned int getItemCount() {
        return items;
    }
    unsigned int getSize() {
        return size;
    }
    unsigned int getFormatSize() {
        return 2 * items * sizeof(int);
    }
    unsigned int getItemElems(int i) {
        return format_array[2 * i + 0];
    }
    unsigned int getItemBytes(int i) {
        return format_array[2 * i + 1];
    }

    void print();

private:
    unsigned int items, size;
    std::vector<unsigned int> format_array;
};

#endif