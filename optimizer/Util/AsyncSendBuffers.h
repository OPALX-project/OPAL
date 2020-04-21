//
// Class AsyncSendBuffer and AsyncSendBuffers
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
#include <vector>
#include <algorithm>
#include <string.h>
#include "mpi.h"
#include "boost/smart_ptr.hpp"
#include "boost/bind.hpp"

class AsyncSendBuffer {

public:

    AsyncSendBuffer(std::ostringstream &os) {
        this->size_req   = new MPI_Request();
        this->buffer_req = new MPI_Request();
        this->buf_size   = os.str().length();
        buffer           = new char[buf_size];
        memcpy(buffer, os.str().c_str(), buf_size);
    }

    ~AsyncSendBuffer() {
        delete   size_req;
        delete   buffer_req;
        delete[] buffer;
    }

    bool hasCompleted() {
        int bufferflag = 0;
        MPI_Test(this->buffer_req, &bufferflag, MPI_STATUS_IGNORE);
        if(bufferflag) {
            int sizeflag = 0;
            MPI_Test(this->buffer_req, &sizeflag, MPI_STATUS_IGNORE);
            if(sizeflag) {
                return true;
            }
        }
        return false;
    }

    void send(int recv_rank, int size_tag, int data_tag, MPI_Comm comm) {
        MPI_Isend(&buf_size, 1, MPI_LONG, recv_rank, size_tag, comm, size_req);
        MPI_Isend(buffer, buf_size, MPI_CHAR, recv_rank, data_tag, comm, buffer_req);
    }


private:

    // can't use smart pointers because MPI will hold last valid reference to
    // pointer
    MPI_Request *size_req;
    MPI_Request *buffer_req;
    char        *buffer;

    size_t buf_size;

};


class AsyncSendBuffers {

public:
    AsyncSendBuffers() {}

    void insert(boost::shared_ptr<AsyncSendBuffer> buf) {
        collection_.push_back(buf);
    }

    void cleanup() {
        collection_.erase(std::remove_if(
            collection_.begin(), collection_.end(), boost::bind(&AsyncSendBuffer::hasCompleted, _1)),
            collection_.end());
    }

    size_t size() {
        return collection_.size();
    }

private:
    std::vector< boost::shared_ptr<AsyncSendBuffer> > collection_;
};

