//
// Global functions MPIHelper
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
#include <iomanip>
#include <limits>
#include <string>

#include "Util/MPIHelper.h"

namespace {

/// length-prefixed so keys may contain arbitrary characters, including spaces
void writeString(std::ostream& os, const std::string& s) {
    os << s.size() << ' ' << s;
}

void readString(std::istream& is, std::string& s) {
    std::size_t n = 0;
    is >> n;
    is.get(); // consume the single separator space
    s.resize(n);
    is.read(&s[0], n);
}

} // namespace

void serialize(Param_t params, std::ostringstream& os) {

    os << std::setprecision(std::numeric_limits<double>::max_digits10);
    os << params.size();
    for (const auto& kv : params) {
        os << ' ';
        writeString(os, kv.first);
        os << ' ' << kv.second;
    }
}

void serialize(reqVarContainer_t reqvars, std::ostringstream& os) {

    os << reqvars.size();
    for (const auto& kv : reqvars) {
        os << ' ';
        writeString(os, kv.first);
        os << ' ';
        kv.second.writeState(os);
    }
}

void deserialize(const char* buffer, std::size_t buf_size, Param_t& params) {

    params.clear();
    std::istringstream is(std::string(buffer, buf_size));

    std::size_t n = 0;
    is >> n;
    for (std::size_t i = 0; i < n; ++i) {
        std::string key;
        readString(is, key);
        double value = 0.0;
        is >> value;
        params[key] = value;
    }
}

void deserialize(const char* buffer, std::size_t buf_size, reqVarContainer_t& reqvars) {

    reqvars.clear();
    std::istringstream is(std::string(buffer, buf_size));

    std::size_t n = 0;
    is >> n;
    for (std::size_t i = 0; i < n; ++i) {
        std::string key;
        readString(is, key);
        reqVarInfo_t info;
        info.readState(is);
        reqvars[key] = info;
    }
}

void MPI_Bcast_params(Param_t& params, std::size_t root, MPI_Comm comm) {

    int pid = 0;
    std::size_t my_pid = 0;
    MPI_Comm_rank(comm, &pid);
    my_pid = static_cast<std::size_t>(pid);

    std::size_t buf_size = 0;
    std::ostringstream os;

    if (my_pid == root) {
        serialize(params, os);
        buf_size = os.str().length();
    }

    MPI_Bcast(&buf_size, 1, MPI_UNSIGNED_LONG, root, comm);

    char *buffer = new char[buf_size];
    if (my_pid == root) memcpy(buffer, os.str().c_str(), buf_size);

    MPI_Bcast(buffer, buf_size, MPI_CHAR, root, comm);
    if (my_pid != root) deserialize(buffer, buf_size, params);

    delete[] buffer;
}

void MPI_Send_params(Param_t params, std::size_t pid, MPI_Comm comm) {

    std::ostringstream os;
    serialize(params, os);
    std::size_t buf_size = os.str().length();

    MPI_Send(&buf_size, 1, MPI_UNSIGNED_LONG, pid,
             MPI_EXCHANGE_SERIALIZED_DATA_TAG, comm);

    char *buffer = new char[buf_size];
    memcpy(buffer, os.str().c_str(), buf_size);

    MPI_Send(buffer, buf_size, MPI_CHAR, pid,
             MPI_EXCHANGE_SERIALIZED_DATA_TAG, comm);

    delete[] buffer;
}

std::pair<std::size_t*, char*> MPI_ISend_params(Param_t params, std::size_t pid,
                                           MPI_Comm comm, MPI_Request *req) {

    std::ostringstream os;
    serialize(params, os);
    std::size_t* buf_size = new std::size_t();
    *buf_size = os.str().length();

    MPI_Isend(buf_size, 1, MPI_UNSIGNED_LONG, pid,
              MPI_EXCHANGE_SERIALIZED_DATA_TAG, comm, req);

    char *buffer = new char[*buf_size];
    memcpy(buffer, os.str().c_str(), *buf_size);

    MPI_Isend(buffer, *buf_size, MPI_CHAR, pid,
              MPI_EXCHANGE_SERIALIZED_DATA_TAG, comm, req);

    std::pair<std::size_t*, char*> p(buf_size, buffer);

    return p;
}

void MPI_Recv_params(Param_t& params, std::size_t pid, MPI_Comm comm) {

    MPI_Status status;
    std::size_t buf_size = 0;
    MPI_Recv(&buf_size, 1, MPI_UNSIGNED_LONG, pid,
             MPI_EXCHANGE_SERIALIZED_DATA_TAG, comm, &status);

    char *buffer = new char[buf_size]();

    MPI_Recv(buffer, buf_size, MPI_CHAR, pid,
             MPI_EXCHANGE_SERIALIZED_DATA_TAG, comm, &status);

    deserialize(buffer, buf_size, params);

    delete[] buffer;
}

void MPI_Send_reqvars(reqVarContainer_t reqvars, std::size_t pid, MPI_Comm comm) {

    std::ostringstream os;
    serialize(reqvars, os);
    std::size_t buf_size = os.str().length();

    MPI_Send(&buf_size, 1, MPI_UNSIGNED_LONG, pid,
             MPI_EXCHANGE_SERIALIZED_DATA_TAG, comm);

    char *buffer = new char[buf_size];
    memcpy(buffer, os.str().c_str(), buf_size);

    MPI_Send(buffer, buf_size, MPI_CHAR, pid,
             MPI_EXCHANGE_SERIALIZED_DATA_TAG, comm);

    delete[] buffer;
}

void MPI_Recv_reqvars(reqVarContainer_t& reqvars, std::size_t pid, MPI_Comm comm) {

    MPI_Status status;
    std::size_t buf_size = 0;
    MPI_Recv(&buf_size, 1, MPI_UNSIGNED_LONG, pid,
             MPI_EXCHANGE_SERIALIZED_DATA_TAG, comm, &status);

    char *buffer = new char[buf_size]();

    MPI_Recv(buffer, buf_size, MPI_CHAR, pid,
             MPI_EXCHANGE_SERIALIZED_DATA_TAG, comm, &status);

    deserialize(buffer, buf_size, reqvars);

    delete[] buffer;
}