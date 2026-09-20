//
// Class FieldArray
//   Storage for one component of a 3D field map.
//
//   Holds the samples either in its own heap buffer (the historical behaviour) or as a
//   read-only view of the HDF5 file mapped into memory. The mapped form exists because
//   every MPI rank otherwise allocates a private copy of every map: two maps of ~5 GB
//   each cost ~11.5 GB per rank, which caps a job by memory long before it runs out of
//   cores. Mapping the file instead lets the page cache hold one physical copy per node,
//   shared by every rank on it, with no MPI involvement.
//
//   The mapped form is only possible when the dataset is stored contiguously, unfiltered,
//   as little-endian float64, at a byte offset the host can align - i.e. when the bytes on
//   disk are already exactly what we would have read into memory. tryMap() checks all of
//   that and returns false if anything does not hold, in which case the caller reads the
//   data the way it always did.
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
#ifndef CLASSIC_FIELDARRAY_H
#define CLASSIC_FIELDARRAY_H

#include <cstddef>
#include <string>
#include <vector>

class FieldArray {
public:
    FieldArray() = default;
    ~FieldArray() { reset(); }

    FieldArray(const FieldArray&) = delete;
    FieldArray& operator=(const FieldArray&) = delete;

    /// Number of samples. Zero for an array that has not been filled.
    std::size_t size() const { return size_m; }
    bool empty() const { return size_m == 0; }

    /// Read access. This is the hot path: one load through a pointer the caller
    /// already has, exactly as indexing a vector was.
    double operator[](std::size_t i) const { return data_m[i]; }

    const double* data() const { return data_m; }

    /// True when the samples are a view of a mapped file rather than our own buffer.
    /// Callers that modify the samples in place must check this.
    bool isMapped() const { return mapping_m != nullptr; }

    /// Allocate a writable buffer of n samples, zero filled. Releases any previous
    /// contents, including a mapping.
    void allocate(std::size_t n);

    /// Writable pointer, for handing to a reader. Null when the array is mapped;
    /// a mapped array is read-only by construction.
    double* mutableData() { return owned_m.empty() ? nullptr : owned_m.data(); }

    /// Try to back this array by a read-only mapping of one dataset of an HDF5 file.
    /// Returns false, leaving the array untouched, if the file cannot be opened, the
    /// dataset is absent, or its layout is not one we can map directly. `why` receives
    /// a short reason on failure, for the log.
    bool tryMap(const std::string& filename,
                const std::string& dataset,
                std::size_t expectedSamples,
                std::string& why);

    /// Release the buffer or the mapping.
    void reset();

private:
    const double* data_m    = nullptr;  ///< what operator[] reads; into owned_m or the mapping
    std::size_t   size_m    = 0;        ///< samples, not bytes
    std::vector<double> owned_m;        ///< used when not mapped
    void*         mapping_m = nullptr;  ///< page-aligned mmap base, null when not mapped
    std::size_t   mapBytes_m = 0;       ///< length passed to mmap, for munmap
};

#endif
