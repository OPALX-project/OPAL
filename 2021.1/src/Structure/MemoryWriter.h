//
// Class MemoryWriter
//   This class writes a SDDS file with virtual memory usage information.
//
// Copyright (c) 2019, Matthias Frey, Paul Scherrer Institut, Villigen PSI, Switzerland
// All rights reserved
//
// Implemented as part of the PhD thesis
// "Precise Simulations of Multibunches in High Intensity Cyclotrons"
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
#ifndef OPAL_MEMORY_WRITER_H
#define OPAL_MEMORY_WRITER_H

#include "SDDSWriter.h"

class MemoryWriter : public SDDSWriter {

public:
    MemoryWriter(const std::string& fname, bool restart);

    void write(const PartBunchBase<double, 3> *beam) override;

private:
    void fillHeader();
};

#endif