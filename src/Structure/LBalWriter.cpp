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
#include "LBalWriter.h"

#include "OPALconfig.h"
#include "AbstractObjects/OpalData.h"
#include "Utilities/Util.h"
#include "Utilities/Timer.h"
#include "Algorithms/PartBunchBase.h"

#ifdef ENABLE_AMR
#include "Algorithms/AmrPartBunch.h"
#endif

LBalWriter::LBalWriter(const std::string& fname, bool restart)
    : SDDSWriter(fname, restart)
{ }


void LBalWriter::fillHeader(PartBunchBase<double, 3> */*beam*/) {

    if (this->hasColumns()) {
        return;
    }

    columns_m.addColumn("t", "double", "ns", "Time");

    for (int p = 0; p < Ippl::getNodes(); ++p) {
        std::stringstream tmp1;
        tmp1 << "\"processor-" << p << "\"";

        std::stringstream tmp2;
        tmp2 << "Number of particles of processor " << p;

        columns_m.addColumn(tmp1.str(), "long", "1", tmp2.str());
    }

#ifdef ENABLE_AMR
    if ( AmrPartBunch* amrbeam = dynamic_cast<AmrPartBunch*>(beam) ) {

        int nLevel = (amrbeam->getAmrObject())->maxLevel() + 1;

        for (int lev = 0; lev < nLevel; ++lev) {
            std::stringstream tmp1;
            tmp1 << "\"level-" << lev << "\"";

            std::stringstream tmp2;
            tmp2 << "Number of particles at level " << lev;
            columns_m.addColumn(tmp1.str(), "long", "1", tmp2.str());
        }
    }
#endif

    if ( mode_m == std::ios::app )
        return;

    OPALTimer::Timer simtimer;

    std::string dateStr(simtimer.date());
    std::string timeStr(simtimer.time());

    std::stringstream ss;
    ss << "Processor statistics '"
       << OpalData::getInstance()->getInputFn() << "' "
       << dateStr << "" << timeStr;

    this->addDescription(ss.str(), "lbal parameters");

    this->addDefaultParameters();


    this->addInfo("ascii", 1);
}


void LBalWriter::write(PartBunchBase<double, 3> *beam) {

    beam->gatherLoadBalanceStatistics();

#ifdef ENABLE_AMR
    if ( AmrPartBunch* amrbeam = dynamic_cast<AmrPartBunch*>(beam) ) {
        amrbeam->gatherLevelStatistics();
    }
#endif

    if ( Ippl::myNode() != 0 )
        return;


    this->fillHeader(beam);

    this->open();

    this->writeHeader();

    columns_m.addColumnValue("t", beam->getT() * 1e9); // 1

    size_t nProcs = Ippl::getNodes();
    for (size_t p = 0; p < nProcs; ++ p) {
        std::stringstream ss;
        ss << "\"processor-" << p << "\"";
        columns_m.addColumnValue(ss.str(), beam->getLoadBalance(p));
    }

#ifdef ENABLE_AMR
    if ( AmrPartBunch* amrbeam = dynamic_cast<AmrPartBunch*>(beam) ) {
        int nLevel = (amrbeam->getAmrObject())->maxLevel() + 1;
        for (int lev = 0; lev < nLevel; ++lev) {
            std::stringstream ss;
            ss << "\"level-" << lev << "\"";
            columns_m.addColumnValue(ss.str(), amrbeam->getLevelStatistics(lev));
        }
    }
#endif

    this->writeRow();

    this->close();
}
