//
// Class GridLBalWriter
//   This class writes a SDDS file with AMR grid load balancing information.
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
#ifndef OPAL_GRID_LBAL_WRITER_H
#define OPAL_GRID_LBAL_WRITER_H

#include "SDDSWriter.h"

class GridLBalWriter : public SDDSWriter {

public:
    GridLBalWriter(const std::string& fname, bool restart);

    void write(PartBunchBase<double, 3> *beam);

private:
    void fillHeader(const PartBunchBase<double, 3> *beam);
};

#endif