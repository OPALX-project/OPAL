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
