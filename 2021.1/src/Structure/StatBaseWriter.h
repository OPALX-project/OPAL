//
// Class StatBaseWriter
//   This common base class for the StatWriter and MultiBunchDump class.
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
#ifndef OPAL_STAT_BASE_WRITER_H
#define OPAL_STAT_BASE_WRITER_H

#include "SDDSWriter.h"
#include "Utilities/Util.h"

class StatBaseWriter : public SDDSWriter {

public:
    StatBaseWriter(const std::string& fname, bool restart);

    /** \brief
     *  delete the last 'numberOfLines' lines of the statistics file
     */
    unsigned int rewindToSpos(double maxSpos);
};


inline
unsigned int StatBaseWriter::rewindToSpos(double maxSPos) {
    if (Ippl::myNode() == 0) {
        return Util::rewindLinesSDDS(this->fname_m, maxSPos);
    }
    return 0;
}

#endif