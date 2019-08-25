#include "Structure/SDDSColumnSet.h"

void SDDSColumnSet::addColumn(const std::string& name,
                              const std::string& type,
                              const std::string& unit,
                              const std::string& desc,
                              std::ios_base::fmtflags flags,
                              unsigned short prec) {

    if (name2idx_m.find(name) != name2idx_m.end()) {
        throw OpalException("SDDSColumnSet::addColumn",
                            "column name '" + name + "' already exists");
    }

    name2idx_m.insert(std::make_pair(name, columns_m.size()));
    columns_m.emplace_back(SDDSColumn(name, type, unit, desc, flags, prec));
}


void SDDSColumnSet::writeHeader(std::ostream& os,
                                const std::string& indent) const {
    for (unsigned int i = 0; i < columns_m.size(); ++ i) {
        auto & col = columns_m[i];
        col.writeHeader(os, i + 1, indent);
    }
}


void SDDSColumnSet::writeRow(std::ostream& os) const {
    for (auto & col: columns_m) {
        os << col;
    }
    os << std::endl;
}