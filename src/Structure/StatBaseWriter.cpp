#include "StatBaseWriter.h"

StatBaseWriter::StatBaseWriter(const std::string& fname, bool restart)
    : SDDSWriter(fname, restart)
{ }
