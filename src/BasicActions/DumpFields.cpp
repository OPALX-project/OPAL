//
// Class DumpFields
//   DumpFields dumps the static magnetic field of a Ring in a user-defined grid
//
// Copyright (c) 2016, Chris Rogers
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
#include "BasicActions/DumpFields.h"

#include "AbstractObjects/OpalData.h"
#include "AbsBeamline/Component.h"
#include "Attributes/Attributes.h"
#include "Fields/Interpolation/ThreeDGrid.h"
#include "Utilities/OpalException.h"
#include "Utilities/Util.h"

#include <boost/filesystem.hpp>

#include <fstream>

extern Inform* gmsg;

std::unordered_set<DumpFields*> DumpFields::dumpsSet_m;

DumpFields::DumpFields() :
    Action(SIZE, "DUMPFIELDS",
           "The \"DUMPFIELDS\" statement dumps a field map to a user-defined "
           "field file, for checking that fields are read in correctly "
           "from disk. The fields are written out on a Cartesian grid.") {
    // would be nice if "steps" could be integer
    itsAttr[FILE_NAME] = Attributes::makeString
        ("FILE_NAME", "Name of the file to which field data is dumped");

    itsAttr[X_START] = Attributes::makeReal
        ("X_START", "Start point in the grid in x [m]");

    itsAttr[DX] = Attributes::makeReal
        ("DX", "Grid step size in x [m]");

    itsAttr[X_STEPS] = Attributes::makeReal
        ("X_STEPS", "Number of steps in x");

    itsAttr[Y_START] = Attributes::makeReal
        ("Y_START", "Start point in the grid in y [m]");

    itsAttr[DY] = Attributes::makeReal
        ("DY", "Grid step size in y [m]");

    itsAttr[Y_STEPS] = Attributes::makeReal
        ("Y_STEPS", "Number of steps in y");

    itsAttr[Z_START] = Attributes::makeReal
        ("Z_START", "Start point in the grid in z [m]");

    itsAttr[DZ] = Attributes::makeReal
        ("DZ", "Grid step size in z [m]");

    itsAttr[Z_STEPS] = Attributes::makeReal
        ("Z_STEPS", "Number of steps in z");

    registerOwnership(AttributeHandler::STATEMENT);
}

DumpFields::DumpFields(const std::string& name, DumpFields* parent):
    Action(name, parent)
{}

DumpFields::~DumpFields() {
    delete grid_m;
    dumpsSet_m.erase(this);
}

DumpFields* DumpFields::clone(const std::string& name) {
    DumpFields* dumper = new DumpFields(name, this);
    if (grid_m != nullptr) {
        dumper->grid_m = grid_m->clone();
    }
    dumper->filename_m = filename_m;
    if (dumpsSet_m.find(this) != dumpsSet_m.end()) {
        dumpsSet_m.insert(dumper);
    }
    return dumper;
}

void DumpFields::execute() {
    buildGrid();
    // the routine for action (OpalParser/OpalParser) calls execute and then
    // deletes 'this'; so we must build a copy that lasts until the field maps
    // are constructed and we are ready for tracking (which is when the field
    // maps are written). Hence the clone call below.
    dumpsSet_m.insert(this->clone(""));
}

void DumpFields::buildGrid() {
    double x0 = Attributes::getReal(itsAttr[X_START]);
    double dx = Attributes::getReal(itsAttr[DX]);
    double nx = Attributes::getReal(itsAttr[X_STEPS]);

    double y0 = Attributes::getReal(itsAttr[Y_START]);
    double dy = Attributes::getReal(itsAttr[DY]);
    double ny = Attributes::getReal(itsAttr[Y_STEPS]);

    double z0 = Attributes::getReal(itsAttr[Z_START]);
    double dz = Attributes::getReal(itsAttr[DZ]);
    double nz = Attributes::getReal(itsAttr[Z_STEPS]);

    checkInt(nx, "X_STEPS");
    checkInt(ny, "Y_STEPS");
    checkInt(nz, "Z_STEPS");
    delete grid_m;

    grid_m = new interpolation::ThreeDGrid(dx, dy, dz,
                                           x0, y0, z0,
                                           nx, ny, nz);

    filename_m = Attributes::getString(itsAttr[FILE_NAME]);
}

void DumpFields::writeFields(Component* field) {
    typedef std::unordered_set<DumpFields*>::iterator dump_iter;
    for (dump_iter it = dumpsSet_m.begin(); it != dumpsSet_m.end(); ++it) {
        (*it)->writeFieldThis(field);
    }
}

void DumpFields::checkInt(double real, std::string name, double tolerance) {
    if (std::abs(std::floor(real) - real) > tolerance) {
        throw OpalException("DumpFields::checkInt",
                            "Value for " + name +
                            " should be an integer but a real value was found");
    }
    if (std::floor(real) < 0.5) {
        throw OpalException("DumpFields::checkInt",
                            "Value for " + name + " should be 1 or more");
    }
}

void DumpFields::writeFieldThis(Component* field) {
    if (grid_m == nullptr) {
        throw OpalException("DumpFields::writeFieldThis",
                            "The grid was nullptr; there was a problem with the DumpFields initialisation.");
    }
    if (field == nullptr) {
        throw OpalException("DumpFields::writeFieldThis",
                            "The field to be written was nullptr.");
    }

    *gmsg << *this << endl;

    std::string fname;
    if (boost::filesystem::path(filename_m).is_absolute() == true) {
        fname = filename_m;
    } else {
        fname = Util::combineFilePath({
            OpalData::getInstance()->getAuxiliaryOutputDirectory(),
            filename_m
        });
    }

    double time = 0.;
    Vector_t point(0., 0., 0.);
    Vector_t centroid(0., 0., 0.);
    std::ofstream fout(fname.c_str(), std::ofstream::out);
    if (!fout.good()) {
        throw OpalException("DumpFields::writeFieldThis",
                            "Failed to open DumpFields file " + filename_m);
    }
    // set precision
    fout << grid_m->end().toInteger() << "\n";
    fout << 1 << " x [m]\n";
    fout << 2 << " y [m]\n";
    fout << 3 << " z [m]\n";
    fout << 4 << " Bx [kGauss]\n";
    fout << 5 << " By [kGauss]\n";
    fout << 6 << " Bz [kGauss]\n";
    fout << 0 << std::endl;
    for (interpolation::Mesh::Iterator it = grid_m->begin();
         it < grid_m->end();
         ++it) {
        Vector_t E(0., 0., 0.);
        Vector_t B(0., 0., 0.);
        it.getPosition(&point[0]);
        field->apply(point, centroid, time, E, B);
        fout << point[0] << " " << point[1] << " " << point[2] << " ";
        fout << B[0] << " " << B[1] << " " << B[2] << "\n";
    }
    if (!fout.good()) {
        throw OpalException("DumpFields::writeFieldThis",
                            "Something went wrong during writing " + filename_m);
    }
    fout.close();
}

void DumpFields::print(std::ostream& os) const {
    os << "* ************* D U M P  F I E L D S *********************************************** " << std::endl;
    os << "* File name: '" << filename_m << "'\n"
       << "* X_START = "  << Attributes::getReal(itsAttr[X_START]) << " [m]\n"
       << "* DX      = "  << Attributes::getReal(itsAttr[DX])      << " [m]\n"
       << "* X_STEPS = "  << Attributes::getReal(itsAttr[X_STEPS]) << '\n'
       << "* Y_START = "  << Attributes::getReal(itsAttr[Y_START]) << " [m]\n"
       << "* DY      = "  << Attributes::getReal(itsAttr[DY])      << " [m]\n"
       << "* Y_STEPS = "  << Attributes::getReal(itsAttr[Y_STEPS]) << '\n'
       << "* Z_START = "  << Attributes::getReal(itsAttr[Z_START]) << " [m]\n"
       << "* DZ      = "  << Attributes::getReal(itsAttr[DZ])      << " [m]\n"
       << "* Z_STEPS = "  << Attributes::getReal(itsAttr[Z_STEPS]) << '\n';
    os << "* ********************************************************************************** " << std::endl;
}
