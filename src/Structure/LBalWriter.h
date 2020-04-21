//
// Class LBalWriter
//   This class writes a SDDS file with MPI load balancing information.
//
// Copyright (c) 2019, Matthias Frey, Paul Scherrer Institut, Villigen PSI, Switzerland
//                     Christof Metzger-Kraus, Open Sourcerer
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
#ifndef OPAL_LBAL_WRITER_H
#define OPAL_LBAL_WRITER_H

#include "SDDSWriter.h"

class LBalWriter : public SDDSWriter {

public:
    LBalWriter(const std::string& fname, bool restart);

    void write(PartBunchBase<double, 3> *beam) override;

private:
    void fillHeader(PartBunchBase<double, 3> *beam);
};


#endif