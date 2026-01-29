#ifndef __SDDSREADER_H__
#define __SDDSREADER_H__

#include "SDDSParser.h"

class SDDSReader: public SDDS::SDDSParser
{
 public:
    SDDSReader(const std::string &fname):
        SDDSParser(fname)
    { }

    inline void parseFile()
    {
        run();
    }
};

#endif