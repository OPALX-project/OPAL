#ifndef OPAL_GRID_LBAL_WRITER_H
#define OPAL_GRID_LBAL_WRITER_H

#include "SDDSWriter.h"

class GridLBalWriter : public SDDSWriter {

public:
    GridLBalWriter(const std::string& fname, bool restart);

    void write(PartBunchBase<double, 3> *beam) override;

private:
    void fillHeader(PartBunchBase<double, 3> *beam);
};

#endif