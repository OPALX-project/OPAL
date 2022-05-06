//
// Class MonitorStatisticsWriter
//   This class writes statistics of monitor element.
//
// Copyright (c) 2019, Christof Metzger-Kraus, Open Sourcerer
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
#ifndef OPAL_MONITOR_STATISTICS_WRITER_H
#define OPAL_MONITOR_STATISTICS_WRITER_H

#include "Structure/SDDSWriter.h"

struct SetStatistics;

class MonitorStatisticsWriter : public SDDSWriter {

public:
    MonitorStatisticsWriter(const std::string& fname, bool restart);

    void addRow(const SetStatistics& set);

private:
    void fillHeader();
};

#endif
