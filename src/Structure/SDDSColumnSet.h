#ifndef SDDSWRITERCOLUMNSET_H
#define SDDSWRITERCOLUMNSET_H

#include "Structure/SDDSColumn.h"
#include "Utilities/OpalException.h"

#include <string>
#include <ostream>
#include <vector>
#include <map>

class SDDSColumnSet {
public:
    SDDSColumnSet();

    void addColumn(const std::string& name,
                   const std::string& type,
                   const std::string& unit,
                   const std::string& desc);

    template<typename T>
    void addColumnValue(const std::string& name,
                        const T& val);

    void writeHeader(std::ostream& os,
                     const std::string& indent) const;

    void writeRow(std::ostream& os) const;

private:
    std::vector<SDDSColumn> columns_m;
    std::map<std::string, unsigned int> name2idx_m;
};


inline
SDDSColumnSet::SDDSColumnSet()
{ }


template<typename T>
void SDDSColumnSet::addColumnValue(const std::string& name,
                                   const T& val) {

    auto it = name2idx_m.find(name);
    if (it == name2idx_m.end()) {
        throw OpalException("SDDSColumnSet::addColumnValue",
                            "column name '" + name + "' doesn't exists");
    }

    auto & col = columns_m[it->second];
    col.addValue(val);
}


#endif