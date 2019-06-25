#ifndef OPAL_SDDS_WRITER_H
#define OPAL_SDDS_WRITER_H

#include <fstream>
#include <string>
#include <queue>
#include <map>
#include <vector>
#include <tuple>
#include <utility>
#include <ostream>
#include <iomanip>

#include <boost/filesystem.hpp>
#include <boost/variant.hpp>

#include "Algorithms/PartBunchBase.h"

class SDDSDataRow {
public:
    SDDSDataRow(const std::vector<std::string>& cN, bool header):
        columnNames_m(cN),
        header_m(header),
        indent_m("        ")
    { }

    template<typename T>
    void addColumn(const std::string& name,
                   const std::string& type,
                   const std::string& unit,
                   const std::string& desc,
                   const T& value);

    std::ostream& write(std::ostream &os) const;

private:
    typedef boost::variant<float,
                           double,
                           long unsigned int,
                           char,
                           std::string> variant_t;
    
    typedef std::tuple<std::string,
                       std::string,
                       std::string> info_t;

    std::vector<std::string> columnNames_m;
    std::map<std::string, variant_t> columnValues_m;
    std::map<std::string, info_t>    columnInfo_m;
    bool header_m;
    std::string indent_m;
};

class SDDSWriter {

public:
    // order: text, content
    typedef std::pair<std::string, std::string> desc_t;

    // order: name, type, description
    typedef std::tuple<std::string,
                       std::string,
                       std::string> param_t;

    // order: mode, no row counts
    typedef std::pair<std::string, size_t> data_t;

    // order: name, type, unit, description
    typedef std::tuple<std::string,
                       std::string,
                       std::string,
                       std::string> cols_t;


    SDDSWriter(const std::string& fname, bool restart);

    virtual void write(PartBunchBase<double, 3>* beam) { };

    /** \brief
     *  delete the last 'numberOfLines' lines of the file 'fileName'
     */
    void rewindLines(size_t numberOfLines);

    void replaceVersionString();


    bool exists();
    
    bool isModeOut();

protected:

    void addDescription(const std::string& text,
                        const std::string& content);

    template<typename T>
    void addParameter(const std::string& name,
                      const std::string& type,
                      const std::string& desc,
                      const T& value);

    void addDefaultParameters();

    void addColumn(const std::string& name,
                   const std::string& type,
                   const std::string& unit,
                   const std::string& desc);

    const std::vector<std::string>& getColumnNames() const;

    void addInfo(const std::string& mode,
                 const size_t& no_row_counts);

    void writeRow(const SDDSDataRow& row);

    template<typename T>
    void writeValue(const T& value);

    void newline();

    void open();

    void close();

    /** \brief Write SDDS header.
     *
     * Writes the appropriate SDDS format header information, The SDDS tools can be used
     * for plotting data.
     */
    void writeHeader();

    std::string fname_m;

    /** \brief First write to the statistics output file.
     *
     * Initially set to std::ios::out so that SDDS format header information is written to file
     * during the first write call to the statistics output file. Variable is then
     * reset to std::ios::app so that header information is only written once.
     */
    std::ios_base::openmode mode_m;

    void writeDescription_m();

    void writeParameters_m();

private:


    void writeColumns_m();

    void writeInfo_m();

    std::ofstream os_m;

    std::string indent_m;

    desc_t desc_m;
    std::queue<param_t> params_m;
    std::queue<std::string> paramValues_m;
    std::queue<cols_t>  cols_m;
    std::vector<std::string> columnNames_m;
    data_t info_m;
};




inline
bool SDDSWriter::exists() {
    return boost::filesystem::exists(fname_m);
}


inline
bool SDDSWriter::isModeOut() {
    return (mode_m == std::ios::out);
}

template<typename T>
void SDDSWriter::writeValue(const T& value) {
    os_m << value << std::setw(10) << "\t";
}


inline
void SDDSWriter::newline() {
    os_m << std::endl;
}


inline
void SDDSWriter::addDescription(const std::string& text,
                                const std::string& content)
{
    desc_m = std::make_pair(text, content);
}


template<typename T>
void SDDSWriter::addParameter(const std::string& name,
                              const std::string& type,
                              const std::string& desc,
                              const T& value)
{
    params_m.push(std::make_tuple(name, type, desc));
    std::stringstream ss;
    ss << value;
    paramValues_m.push(ss.str());
}


inline
const std::vector<std::string>& SDDSWriter::getColumnNames() const {
    return columnNames_m;
}


inline
void SDDSWriter::addInfo(const std::string& mode,
                         const size_t& no_row_counts) {
    info_m = std::make_pair(mode, no_row_counts);
}


inline
void SDDSWriter::writeRow(const SDDSDataRow& row) {
    row.write(os_m);
}


inline
void SDDSWriter::addColumn(const std::string& name,
                           const std::string& type,
                           const std::string& unit,
                           const std::string& desc)
{
    cols_m.push(std::make_tuple(name, type, unit, desc));
    columnNames_m.push_back(name);
}

template<typename T>
void SDDSDataRow::addColumn(const std::string& name,
                            const std::string& type,
                            const std::string& unit,
                            const std::string& desc,
                            const T& value) {
    if (columnValues_m.find(name) != columnValues_m.end()) return;

    variant_t var = value;
    columnValues_m.insert(std::make_pair(name, var));
    
    if ( header_m )
        columnInfo_m.insert(std::make_pair(name, std::make_tuple(type, unit, desc)));
}

inline
std::ostream& SDDSDataRow::write(std::ostream& os) const {
    if (!columnInfo_m.empty()) {
        int column = 1;
        for (const auto name: columnNames_m) {
            auto it = columnInfo_m.find(name);
            if (it == columnInfo_m.end()) {
                throw OpalException("SDDSDataRow::write",
                                    "info for column '" + name + "' not provided");
            }
            os << "&column\n"
               << indent_m << "name=" << name << ",\n"
               << indent_m << "type=" << std::get<0>(it->second) << ",\n"
               << indent_m << "units=" << std::get<1>(it->second) << ",\n"
               << indent_m << "description=\"" << column++ << " " << std::get<2>(it->second) << "\"\n"
               << "&end\n";
        }
        os << "&data\n"
           << indent_m << "mode=ascii" << ",\n"
           << indent_m << "no_row_counts=1" << "\n"
           << "&end";
    }

    for (const auto name: columnNames_m) {
        auto it = columnValues_m.find(name);
        if (it == columnValues_m.end()) {
            throw OpalException("SDDSDataRow::write",
                                "value for column '" + name + "' not provided");
        }

        os << it->second << std::setw(10) << "\t";
    }

    return os;
}

inline
std::ostream& operator<<(std::ostream& os, const SDDSDataRow& row) {
    return row.write(os);
}

#endif