//
// Class MultiBunchDump
//   This class writes a SDDS file of single bunch statistics in multibunch simulations.
//
// Copyright (c) 2018 - 2019, Matthias Frey, Paul Scherrer Institut, Villigen PSI, Switzerland
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
#ifndef OPAL_MULTI_BUNCH_DUMP_H
#define OPAL_MULTI_BUNCH_DUMP_H

#include "StatBaseWriter.h"

#include "Algorithms/MultiBunchHandler.h"


class MultiBunchDump : public StatBaseWriter {

public:
    typedef MultiBunchHandler::beaminfo_t beaminfo_t;

    MultiBunchDump(const std::string& fname, bool restart);

    void fillHeader();

    using SDDSWriter::write;
    void write(PartBunchBase<double, 3>* beam, const beaminfo_t& binfo);
};

#endif
