#include "Structure/SDDSColumn.h"
#include "Utilities/OpalException.h"

#include <iomanip>

SDDSColumn::SDDSColumn(const std::string& name,
                       const std::string& type,
                       const std::string& unit,
                       const std::string& desc):
    name_m(name),
    description_m(std::make_tuple(type, unit, desc)),
    set_m(false)
{ }


void SDDSColumn::writeHeader(std::ostream& os,
                             unsigned int colNr,
                             const std::string& indent) const {
    os << "&column\n"
       << indent << "name=" << name_m << ",\n"
       << indent << "type=" << std::get<0>(description_m) << ",\n";

    if (std::get<1>(description_m) != "")
        os << indent << "units=" << std::get<1>(description_m) << ",\n";

    os << indent << "description=\"" << colNr << " " << std::get<2>(description_m) << "\"\n"
       << "&end\n";
}


void SDDSColumn::writeValue(std::ostream& os) const {
    if (!set_m) {
        throw OpalException("SDDSColumn::writeValue",
                            "value for column '" + name_m + "' isn't set");
    }

    os << value_m << std::setw(10) << "\t";
    set_m = false;
}

std::ostream& operator<<(std::ostream& os,
                         const SDDSColumn& col) {
    col.writeValue(os);

    return os;
}