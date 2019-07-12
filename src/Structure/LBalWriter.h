#ifndef OPAL_LBAL_WRITER_H
#define OPAL_LBAL_WRITER_H

#include "SDDSWriter.h"

class LBalWriter : public SDDSWriter {

public:
    LBalWriter(const std::string& fname, bool restart);

    void write(PartBunchBase<double, 3> *beam) override;

private:
    void fillHeader(PartBunchBase<double, 3> *beam);
};


#endif