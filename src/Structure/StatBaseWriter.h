#ifndef OPAL_STAT_BASE_WRITER_H
#define OPAL_STAT_BASE_WRITER_H

#include "SDDSWriter.h"
#include "Utilities/Util.h"

class StatBaseWriter : public SDDSWriter {

public:
    StatBaseWriter(const std::string& fname, bool restart);

    /** \brief
     *  delete the last 'numberOfLines' lines of the statistics file
     */
    unsigned int rewindToSpos(double maxSpos);
};


inline
unsigned int StatBaseWriter::rewindToSpos(double maxSPos) {
    if (Ippl::myNode() == 0) {
        return Util::rewindLinesSDDS(this->fname_m, maxSPos);
    }
    return 0;
}

#endif