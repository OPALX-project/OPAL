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
#include <sstream>

#include <boost/filesystem.hpp>
#include <boost/variant.hpp>

#include "Structure/SDDSColumn.h"
#include "Structure/SDDSColumnSet.h"

template <class T, unsigned Dim>
class PartBunchBase;

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

    virtual ~SDDSWriter() { };

    virtual void write(PartBunchBase<double, 3>* /*beam*/) { };

    /** \brief
     *  delete the last 'numberOfLines' lines of the file 'fileName'
     */
    void rewindLines(size_t numberOfLines);

    void replaceVersionString();

    double getLastValue(const std::string& column);

    bool exists() const;

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

    void addInfo(const std::string& mode,
                 const size_t& no_row_counts);

    void writeRow();

    void open();

    void close();

    /** \brief Write SDDS header.
     *
     * Writes the appropriate SDDS format header information, The SDDS tools can be used
     * for plotting data.
     */
    void writeHeader();

    template<typename T>
    std::string toString(const T& val);

    std::string fname_m;

    /** \brief First write to the statistics output file.
     *
     * Initially set to std::ios::out so that SDDS format header information is written to file
     * during the first write call to the statistics output file. Variable is then
     * reset to std::ios::app so that header information is only written once.
     */
    std::ios_base::openmode mode_m;

    SDDSColumnSet columns_m;

    /* The columns of the SDDS file need to be
     * added only once.
     */
    bool hasColumns() const;

private:

    void writeDescription();

    void writeParameters();

    void writeColumns();

    void writeInfo();

    std::ofstream os_m;

    std::string indent_m;

    desc_t desc_m;
    std::queue<param_t> params_m;
    std::queue<std::string> paramValues_m;
    data_t info_m;

    static constexpr
    unsigned int precision_m = 15;
};


inline
bool SDDSWriter::exists() const {
    return boost::filesystem::exists(fname_m);
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
void SDDSWriter::addInfo(const std::string& mode,
                         const size_t& no_row_counts) {
    info_m = std::make_pair(mode, no_row_counts);
}


inline
void SDDSWriter::writeRow() {
    columns_m.writeRow(os_m);
}


template<typename T>
std::string SDDSWriter::toString(const T& val) {
    std::ostringstream ss;
    ss.precision(precision_m);
    ss << val;
    return ss.str();
}


inline
bool SDDSWriter::hasColumns() const {
    return columns_m.hasColumns();
}

#endif