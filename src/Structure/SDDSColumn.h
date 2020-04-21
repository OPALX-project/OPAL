//
// Class SDDSColumn
//   This class writes column entries of SDDS files.
//
// Copyright (c) 2019, Christof Metzger-Kraus, Open Sourcerer
// All rights reserved
//
// This file is part of OPAL.
//
// OPAL is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// You should have received a copy of the GNU General Public License
// along with OPAL. If not, see <https://www.gnu.org/licenses/>.
//
#ifndef SDDSWRITERCOLUMN_H
#define SDDSWRITERCOLUMN_H

#include <boost/variant.hpp>

#include <ostream>
#include <tuple>
#include <string>

class SDDSColumn {
public:
    SDDSColumn(const std::string& name,
               const std::string& type,
               const std::string& unit,
               const std::string& desc,
               std::ios_base::fmtflags flags,
               unsigned short precision);

    template<typename T>
    void addValue(const T& val);

    void writeHeader(std::ostream& os,
                     unsigned int colNr,
                     const std::string& indent) const;

protected:

    void writeValue(std::ostream& os) const;

private:
    friend
    std::ostream& operator<<(std::ostream& os,
                             const SDDSColumn& col);

    typedef std::tuple<std::string,
                       std::string,
                       std::string> desc_t;

    typedef boost::variant<float,
                           double,
                           long unsigned int,
                           char,
                           std::string> variant_t;
    std::string name_m;
    desc_t description_m;
    variant_t value_m;

    std::ios_base::fmtflags writeFlags_m;
    unsigned short writePrecision_m;

    mutable bool set_m;
};


template<typename T>
void SDDSColumn::addValue(const T& val) {
    value_m = val;
    set_m = true;
}


std::ostream& operator<<(std::ostream& os,
                             const SDDSColumn& col);

#endif