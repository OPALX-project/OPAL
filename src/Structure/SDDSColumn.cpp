#include "Structure/SDDSColumn.h"
#include "Utilities/OpalException.h"

#include <iomanip>
#include <list>

SDDSColumn::SDDSColumn(const std::string& name,
                       const std::string& type,
                       const std::string& unit,
                       const std::string& desc,
                       std::ios_base::fmtflags flags,
                       unsigned short prec):
    name_m(name),
    description_m(std::make_tuple(type, unit, desc)),
    writeFlags_m(flags),
    writePrecision_m(prec),
    set_m(false)
{
    std::list<std::ios_base::fmtflags> numericalBase({std::ios_base::dec,
                                                      std::ios_base::hex,
                                                      std::ios_base::oct});
    std::list<std::ios_base::fmtflags> floatFormat({std::ios_base::fixed,
                                                    std::ios_base::scientific});
    std::list<std::ios_base::fmtflags> adjustmentFlags({std::ios_base::internal,
                                                        std::ios_base::left,
                                                        std::ios_base::right});

    // This code ensures that for each group of flags only one flag is given
    for (std::ios_base::fmtflags flag: numericalBase) {
        if (writeFlags_m & flag) {
            writeFlags_m = (flag | (writeFlags_m & ~std::ios_base::basefield));
            break;
        }
    }
    for (std::ios_base::fmtflags flag: floatFormat) {
        if (writeFlags_m & flag) {
            writeFlags_m = (flag | (writeFlags_m & ~std::ios_base::floatfield));
            break;
        }
    }
    for (std::ios_base::fmtflags flag: adjustmentFlags) {
        if (writeFlags_m & flag) {
            writeFlags_m = (flag | (writeFlags_m & ~std::ios_base::adjustfield));
            break;
        }
    }

}


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

    os.flags(writeFlags_m);
    os.precision(writePrecision_m);
    os << value_m << std::setw(10) << "\t";
    set_m = false;
}

std::ostream& operator<<(std::ostream& os,
                         const SDDSColumn& col) {
    col.writeValue(os);

    return os;
}