#ifndef OPAL_MEMORY_WRITER_H
#define OPAL_MEMORY_WRITER_H

#include "SDDSWriter.h"

class MemoryWriter : public SDDSWriter {

public:
    MemoryWriter(const std::string& fname, bool restart);

    void write(PartBunchBase<double, 3> *beam) override;

private:
    void fillHeader();
};

#endif