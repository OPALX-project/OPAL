#ifndef OPAL_MONITOR_STATISTICS_WRITER_H
#define OPAL_MONITOR_STATISTICS_WRITER_H

#include "SDDSWriter.h"

struct SetStatistics;

class MonitorStatisticsWriter : public SDDSWriter {

public:
    MonitorStatisticsWriter(const std::string& fname, bool restart);

    void addRow(const SetStatistics &set);

private:
    void fillHeader();
};


#endif