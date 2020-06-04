//
// Class FieldWriter
//   This class writes the bunch internal fields on the grid to
//   file. It supports single core execution only.
//
// Copyright (c) 2020, Paul Scherrer Institut, Villigen PSI, Switzerland
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
#include <iomanip>
#include <fstream>

#include <boost/filesystem.hpp>
#include <boost/format.hpp>

#include "Utilities/Util.h"
#include "AbstractObjects/OpalData.h"
#include "Algorithms/PBunchDefs.h"

template<typename FieldType>
void FieldWriter::dumpField(FieldType& field, std::string name,
                            std::string unit, long long step,
                            FieldType* image)
{
    if (Ippl::getNodes() > 1) {
        return;
    }

    constexpr bool isVectorField = std::is_same<VField_t, FieldType>::value;
    std::string type = (isVectorField) ? "field" : "scalar";

    INFOMSG("*** START DUMPING " + Util::toUpper(name) + " FIELD ***" << endl);

    /* Save the files in the output directory of the simulation. The file
     * name of vector fields is
     *
     * 'basename'-'name'_field-'******'.dat
     *
     * and of scalar fields
     *
     * 'basename'-'name'_scalar-'******'.dat
     *
     * with
     *   'basename': OPAL input file name (*.in)
     *   'name':     field name (input argument of function)
     *   '******':   step padded with zeros to 6 digits
     */
    std::string dirname = OpalData::getInstance()->getAuxiliaryOutputDirectory();
    boost::filesystem::path file(dirname);
    boost::format filename("%1%-%2%-%|3$06|.dat");
    std::string basename = OpalData::getInstance()->getInputBasename();
    filename % basename % (name + std::string("_") + type) % step;
    file /= filename.str();

    std::ofstream fout(file.string(), std::ios::out);
    fout.precision(9);

    fout << "# " << name << " " << type << " data on grid" << std::endl
         << "#"
         << std::setw(4)  << "i"
         << std::setw(5)  << "j"
         << std::setw(5)  << "k"
         << std::setw(17) << "x [m]"
         << std::setw(17) << "y [m]"
         << std::setw(17) << "z [m]";
    if (isVectorField) {
        fout << std::setw(10) << name << "x [" << unit << "]"
             << std::setw(10) << name << "y [" << unit << "]"
             << std::setw(10) << name << "z [" << unit << "]";
    } else {
        fout << std::setw(13) << name << " [" << unit << "]";
    }

    if (image) {
        fout << std::setw(13) << name << " image [" << unit << "]";
    }

    fout << std::endl;

    Vector_t origin = field.get_mesh().get_origin();
    Vector_t spacing(field.get_mesh().get_meshSpacing(0),
                     field.get_mesh().get_meshSpacing(1),
                     field.get_mesh().get_meshSpacing(2));

    NDIndex<3> localIdx = field.getLayout().getLocalNDIndex();
    for (int x = localIdx[0].first(); x <= localIdx[0].last(); x++) {
        for (int y = localIdx[1].first(); y <= localIdx[1].last(); y++) {
            for (int z = localIdx[2].first(); z <= localIdx[2].last(); z++) {
                NDIndex<3> idx(Index(x, x), Index(y, y), Index(z, z));
                fout << std::setw(5) << x + 1
                     << std::setw(5) << y + 1
                     << std::setw(5) << z + 1
                     << std::setw(17) << origin(0) + x * spacing(0)
                     << std::setw(17) << origin(1) + y * spacing(1)
                     << std::setw(17) << origin(2) + z * spacing(2);
                if (isVectorField) {
                    Vector_t vfield = field.localElement(idx);
                    fout << std::setw(17) << vfield[0]
                         << std::setw(17) << vfield[1]
                         << std::setw(17) << vfield[2];
                } else {
                    fout << std::setw(17) << field.localElement(idx);
                }

                if (image) {
                    fout << std::setw(17) << image->localElement(idx);
                }
                fout << std::endl;
            }
        }
    }
    fout.close();
    INFOMSG("*** FINISHED DUMPING " + Util::toUpper(name) + " FIELD ***" << endl);
}