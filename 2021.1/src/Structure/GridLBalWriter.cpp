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
#include "GridLBalWriter.h"

#include "Utilities/Timer.h"

#include <sstream>

#include "Algorithms/AmrPartBunch.h"

GridLBalWriter::GridLBalWriter(const std::string& fname, bool restart)
    : SDDSWriter(fname, restart)
{ }


void GridLBalWriter::fillHeader(const PartBunchBase<double, 3> *beam) {

    if (this->hasColumns()) {
        return;
    }

    columns_m.addColumn("t", "double", "ns", "Time");

    const AmrPartBunch* amrbeam = dynamic_cast<const AmrPartBunch*>(beam);

    if ( !amrbeam )
        throw OpalException("GridLBalWriter::fillHeader()",
                            "Can not write grid load balancing for non-AMR runs.");
    int nLevel = (amrbeam->getAmrObject())->maxLevel() + 1;

    for (int lev = 0; lev < nLevel; ++lev) {
        std::stringstream tmp1;
        tmp1 << "\"level-" << lev << "\"";

        std::stringstream tmp2;
        tmp2 << "Number of boxes at level " << lev;

        columns_m.addColumn(tmp1.str(), "long", "1", tmp2.str());
    }

    for (int p = 0; p < Ippl::getNodes(); ++p) {
        std::stringstream tmp1;
        tmp1 << "\"processor-" << p << "\"";

        std::stringstream tmp2;
        tmp2 << "Number of grid points per processor " << p;

        columns_m.addColumn(tmp1.str(), "long", "1", tmp2.str());
    }

    if ( mode_m == std::ios::app )
        return;

    OPALTimer::Timer simtimer;

    std::string dateStr(simtimer.date());
    std::string timeStr(simtimer.time());

    std::stringstream ss;
    ss << "Grid load balancing statistics '"
       << OpalData::getInstance()->getInputFn() << "' "
       << dateStr << "" << timeStr;

    this->addDescription(ss.str(), "grid lbal parameters");

    this->addDefaultParameters();

    this->addInfo("ascii", 1);
}


void GridLBalWriter::write(PartBunchBase<double, 3> *beam) {
    AmrPartBunch* amrbeam = dynamic_cast<AmrPartBunch*>(beam);

    if ( !amrbeam )
        throw OpalException("GridLBalWriter::write()",
                            "Can not write grid load balancing for non-AMR runs.");

    std::map<int, long> gridPtsPerCore;
    std::vector<int> gridsPerLevel;

    amrbeam->getAmrObject()->getGridStatistics(gridPtsPerCore, gridsPerLevel);

    if ( Ippl::myNode() != 0 )
        return;

    this->fillHeader(beam);

    this->open();

    this->writeHeader();

    columns_m.addColumnValue("t", beam->getT() * 1e9); // 1

    int nLevel = (amrbeam->getAmrObject())->maxLevel() + 1;

    for (int lev = 0; lev < nLevel; ++lev) {
        std::stringstream ss;
        ss << "\"level-" << lev << "\"";
        // we need to cast due to boost::variant
        columns_m.addColumnValue(ss.str(), toString(gridsPerLevel[lev]));
    }

    int nProcs = Ippl::getNodes();
    for (int p = 0; p < nProcs; ++p) {
        std::stringstream ss;
        ss << "\"processor-" << p << "\"";
        // we need to cast due to boost::variant
        columns_m.addColumnValue(ss.str(), toString(gridPtsPerCore[p]));
    }

    this->writeRow();

    this->close();
}