#include "ElementPositionWriter.h"

#include "Ippl.h"
#include "AbstractObjects/OpalData.h"

ElementPositionWriter::ElementPositionWriter(const std::string& fname)
    : SDDSWriter(fname, false)
{ }


void ElementPositionWriter::fillHeader() {

    if (this->hasColumns()) {
        return;
    }

    this->addDescription("Element positions '" + OpalData::getInstance()->getInputFn() + "'",
                         "element positions");
    this->addDefaultParameters();

    columns_m.addColumn("s",
                        "double",
                        "m",
                        "longitudinal position",
                        std::ios_base::fixed,
                        10);
    columns_m.addColumn("dipole",
                        "float",
                        "1",
                        "dipole field present",
                        std::ios_base::fixed,
                        4);
    columns_m.addColumn("quadrupole",
                        "float",
                        "1",
                        "quadrupole field present",
                        std::ios_base::fixed,
                        4);
    columns_m.addColumn("sextupole",
                        "float",
                        "1",
                        "sextupole field present",
                        std::ios_base::fixed,
                        4);
    columns_m.addColumn("octupole",
                        "float",
                        "1",
                        "octupole field present",
                        std::ios_base::fixed,
                        4);
    columns_m.addColumn("decapole",
                        "float",
                        "1",
                        "decapole field present",
                        std::ios_base::fixed,
                        4);
    columns_m.addColumn("multipole",
                        "float",
                        "1",
                        "higher multipole field present",
                        std::ios_base::fixed,
                        4);
    columns_m.addColumn("solenoid",
                        "float",
                        "1",
                        "solenoid field present",
                        std::ios_base::fixed,
                        4);
    columns_m.addColumn("rfcavity",
                        "float",
                        "1",
                        "RF field present",
                        std::ios_base::fixed,
                        4);
    columns_m.addColumn("monitor",
                        "float",
                        "1",
                        "monitor present",
                        std::ios_base::fixed,
                        4);
    columns_m.addColumn("other",
                        "float",
                        "1",
                        "other element present",
                        std::ios_base::fixed,
                        4);
    columns_m.addColumn("element_names",
                        "string",
                        "",
                        "names of elements");

    this->addInfo("ascii", 1);
}


void ElementPositionWriter::addRow(double spos,
                                   const std::vector<double> &row,
                                   const std::string &elements) {

    if ( Ippl::myNode() != 0 )
        return;


    this->fillHeader();

    this->open();

    this->writeHeader();

    static const std::vector<double> typeMultipliers = {3.3333e-1,
                                                        1.0,
                                                        0.5,
                                                        0.25,
                                                        1.0,
                                                        1.0,
                                                        1.0,
                                                        1.0,
                                                        1.0,
                                                        1.0};

    static const std::vector<std::string> columnNames = {"dipole",
                                                         "quadrupole",
                                                         "sextupole",
                                                         "octupole",
                                                         "decapole",
                                                         "multipole",
                                                         "solenoid",
                                                         "rfcavity",
                                                         "monitor",
                                                         "other"};

    columns_m.addColumnValue("s", spos);
    for (unsigned int i = 0; i < columnNames.size(); ++ i) {
        columns_m.addColumnValue(columnNames[i], row[i] * typeMultipliers[i]);
    }
    columns_m.addColumnValue("element_names", elements);

    this->writeRow();

    this->close();
}