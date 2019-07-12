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
               const std::string& desc);

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