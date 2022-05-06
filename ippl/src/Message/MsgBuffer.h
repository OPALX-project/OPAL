//
// Class MsgBuffer
// MsgBuffer class to allow serializing message objects into plain buffers
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
#ifndef MSGBUFFER_H
#define MSGBUFFER_H

#include <algorithm>
#include <cstring>
#include <vector>
#include "Message/Message.h"

class Format;

class MsgBuffer {
public:
    // creates buffer with space to hold count messages of format f
    MsgBuffer(Format* f, int count, int offset = 0);
    MsgBuffer(Format* f, char* d, int size);

    bool add(Message*);
    Message* get();

    template <class T>
    void get(T& v) {
        std::memcpy(&v, data.data() + readpos, sizeof(T));
        readpos += sizeof(T);
    }

    template <class T>
    void put(T& v) {
        std::memcpy(data.data() + writepos, &v, sizeof(T));
        writepos += sizeof(T);
    }

    int getSize() {
        return writepos;
    }

    void* getBuffer() {
        char* data_ptr = data.empty() ? static_cast<char*>(0) : &(data[0]);
        return data_ptr;
    }

    Format* getFormat() {
        return format;
    }

    ~MsgBuffer();

private:
    Format* format;
    unsigned int datasize, writepos, readpos;
    std::vector<char> data;
};

#endif