#include "MemoryWriter.h"

#include "AbstractObjects/OpalData.h"
#include "Algorithms/PartBunchBase.h"
#include "Utilities/Timer.h"
#include "Ippl.h"

MemoryWriter::MemoryWriter(const std::string& fname, bool restart)
    : SDDSWriter(fname, restart)
{ }


void MemoryWriter::fillHeader() {

    if (this->hasColumns()) {
        return;
    }

    columns_m.addColumn("t", "double", "ns", "Time");

    columns_m.addColumn("s", "double", "m", "Path length");

    IpplMemoryUsage::IpplMemory_p memory = IpplMemoryUsage::getInstance();
    columns_m.addColumn("memory", "double", memory->getUnit(), "Total Memory");

    for (int p = 0; p < Ippl::getNodes(); ++p) {
        std::stringstream tmp1;
        tmp1 << "\"processor-" << p << "\"";

        std::stringstream tmp2;
        tmp2 << "Memory per processor " << p;
        columns_m.addColumn(tmp1.str(), "double", memory->getUnit(), tmp2.str());
    }

    if ( mode_m == std::ios::app )
        return;

    OPALTimer::Timer simtimer;

    std::string dateStr(simtimer.date());
    std::string timeStr(simtimer.time());


    std::stringstream ss;

    ss << "Memory statistics '"
       << OpalData::getInstance()->getInputFn() << "' "
       << dateStr << "" << timeStr;

    this->addDescription(ss.str(), "memory parameters");

    this->addDefaultParameters();


    this->addInfo("ascii", 1);
}


void MemoryWriter::write(PartBunchBase<double, 3> *beam)
{
    IpplMemoryUsage::IpplMemory_p memory = IpplMemoryUsage::getInstance();
    memory->sample();

    if (Ippl::myNode() != 0) {
        return;
    }

    double  pathLength = 0.0;
    if (OpalData::getInstance()->isInOPALCyclMode())
        pathLength = beam->getLPath();
    else
        pathLength = beam->get_sPos();

    fillHeader();

    this->open();

    this->writeHeader();

    columns_m.addColumnValue("t", beam->getT() * 1e9);             // 1
    columns_m.addColumnValue("s", pathLength);                     // 2

    int nProcs = Ippl::getNodes();
    double total = 0.0;
    for (int p = 0; p < nProcs; ++p) {
        total += memory->getMemoryUsage(p);
    }

    columns_m.addColumnValue("memory", total);

    for (int p = 0; p < nProcs; p++) {
        std::stringstream ss;
        ss << "\"processor-" << p << "\"";
        columns_m.addColumnValue(ss.str(),  memory->getMemoryUsage(p));
    }

    this->writeRow();

    this->close();
}
