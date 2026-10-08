//
// Class FieldArray
//   See FieldArray.h.
//
// Copyright (c) 2026, Paul Scherrer Institut, Villigen PSI, Switzerland
// All rights reserved.
//
// This file is part of OPAL.
//
// OPAL is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// You should have received a copy of the GNU General Public License
// along with OPAL.  If not, see <https://www.gnu.org/licenses/>.
//

#include "Fields/FieldArray.h"

#include <hdf5.h>

#include <fcntl.h>
#include <sys/mman.h>
#include <unistd.h>

#include <algorithm>
#include <cstdint>

namespace {
    // The access pattern during tracking is an 8-corner gather per particle, i.e.
    // effectively random over a multi-GB region, so tell the kernel not to bother with
    // readahead.
    constexpr int trackingAdvice = MADV_RANDOM;
}

void FieldArray::allocate(std::size_t n) {
    reset();
    owned_m.assign(n, 0.0);
    data_m = owned_m.data();
    size_m = n;
}

void FieldArray::reset() {
    if (mapping_m != nullptr) {
        munmap(mapping_m, mapBytes_m);
        mapping_m  = nullptr;
        mapBytes_m = 0;
    }
    std::vector<double>().swap(owned_m);
    data_m = nullptr;
    size_m = 0;
}

bool FieldArray::tryMap(const std::string& filename,
                        const std::string& dataset,
                        std::size_t expectedSamples,
                        std::string& why) {
    // Everything here is best effort: any failure means "read it the old way", never
    // an exception. HDF5's own error stack is silenced so a missing dataset does not
    // print a wall of text on a path we expect to fail sometimes.
    H5E_auto2_t oldFunc = nullptr;
    void* oldData = nullptr;
    H5Eget_auto2(H5E_DEFAULT, &oldFunc, &oldData);
    H5Eset_auto2(H5E_DEFAULT, nullptr, nullptr);
    struct ErrRestore {
        H5E_auto2_t f; void* d;
        ~ErrRestore() { H5Eset_auto2(H5E_DEFAULT, f, d); }
    } errRestore{oldFunc, oldData};

    // Serial (POSIX) VFD on purpose: we only want metadata, every rank asks
    // independently, and this must not be collective over the tracking communicator.
    hid_t file = H5Fopen(filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
    if (file < 0) { why = "cannot open serially"; return false; }

    hid_t dset = H5Dopen2(file, dataset.c_str(), H5P_DEFAULT);
    if (dset < 0) { H5Fclose(file); why = "no dataset " + dataset; return false; }

    bool ok = false;
    haddr_t offset = HADDR_UNDEF;
    do {
        hid_t dcpl = H5Dget_create_plist(dset);
        if (dcpl < 0) { why = "no create plist"; break; }
        const H5D_layout_t layout = H5Pget_layout(dcpl);
        const int nfilters = H5Pget_nfilters(dcpl);
        H5Pclose(dcpl);
        if (layout != H5D_CONTIGUOUS) { why = "not contiguous"; break; }
        if (nfilters != 0)            { why = "filtered"; break; }

        // The in-file representation must be bit-identical to a double in memory,
        // otherwise H5Dread's type conversion is doing real work and we cannot skip it.
        hid_t dtype = H5Dget_type(dset);
        const htri_t sameType = H5Tequal(dtype, H5T_IEEE_F64LE);
        H5Tclose(dtype);
        if (sameType <= 0) { why = "not little-endian float64"; break; }

        if (H5Dget_storage_size(dset) != expectedSamples * sizeof(double)) {
            why = "storage size mismatch"; break;
        }

        offset = H5Dget_offset(dset);
        if (offset == HADDR_UNDEF) { why = "no file offset"; break; }
        // mmap can start at any page, but the samples themselves must land on a
        // natural double boundary. HDF5 does not promise that.
        if (offset % sizeof(double) != 0) { why = "offset not 8-byte aligned"; break; }
        ok = true;
    } while (false);

    H5Dclose(dset);
    H5Fclose(file);
    if (!ok) return false;

    const int fd = open(filename.c_str(), O_RDONLY);
    if (fd < 0) { why = "cannot open for mapping"; return false; }

    const std::size_t pageSize = static_cast<std::size_t>(sysconf(_SC_PAGESIZE));
    const off_t base  = static_cast<off_t>((offset / pageSize) * pageSize);
    const std::size_t delta = static_cast<std::size_t>(offset - static_cast<haddr_t>(base));
    const std::size_t bytes = expectedSamples * sizeof(double) + delta;

    // MAP_SHARED so every process on the node faults the same page-cache pages; MAP_PRIVATE
    // would be copy-on-write and defeat the entire purpose. PROT_READ keeps it read-only,
    // so a stray write is a fault rather than silent corruption of the file.
    void* p = mmap(nullptr, bytes, PROT_READ, MAP_SHARED, fd, base);
    close(fd);   // the mapping keeps its own reference to the file
    if (p == MAP_FAILED) { why = "mmap failed"; return false; }

    madvise(p, bytes, trackingAdvice);

    reset();
    mapping_m  = p;
    mapBytes_m = bytes;
    data_m     = reinterpret_cast<const double*>(static_cast<const char*>(p) + delta);
    size_m     = expectedSamples;
    return true;
}

bool FieldArray::anyNonZero(std::size_t begin, std::size_t end) const {
    end = std::min(end, size_m);
    if (begin >= end) {
        return false;
    }

    void* adviseBase = nullptr;
    std::size_t adviseBytes = 0;
    if (mapping_m != nullptr) {
        // madvise wants a page-aligned start; rounding down stays inside the mapping,
        // whose base is page aligned.
        const std::uintptr_t pageMask = static_cast<std::uintptr_t>(sysconf(_SC_PAGESIZE)) - 1;
        const std::uintptr_t lo = reinterpret_cast<std::uintptr_t>(data_m + begin) & ~pageMask;
        const std::uintptr_t hi = reinterpret_cast<std::uintptr_t>(data_m + end);
        adviseBase  = reinterpret_cast<void*>(lo);
        adviseBytes = hi - lo;
        madvise(adviseBase, adviseBytes, MADV_SEQUENTIAL);
    }

    // A block at a time, so the inner loop has no early exit and vectorises.
    // x != 0.0 is also true for NaN, which must not be mistaken for an empty field.
    constexpr std::size_t block = 4096;
    bool found = false;
    for (std::size_t i = begin; i < end && !found; i += block) {
        const std::size_t stop = std::min(i + block, end);
        for (std::size_t j = i; j < stop; ++j) {
            found |= (data_m[j] != 0.0);
        }
    }

    if (adviseBase != nullptr) {
        madvise(adviseBase, adviseBytes, trackingAdvice);
    }
    return found;
}
