//
// Class DumpEMFields
//   DumpEMFields dumps the dynamically changing fields of a Ring in a user-
//   defined grid.
//
// Copyright (c) 2017, Chris Rogers
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
#include "BasicActions/DumpEMFields.h"

#include "AbstractObjects/OpalData.h"
#include "AbsBeamline/Component.h"
#include "Attributes/Attributes.h"
#include "Fields/Interpolation/NDGrid.h"
#include "Physics/Units.h"
#include "Utilities/OpalException.h"
#include "Utilities/Util.h"

#include <boost/filesystem.hpp>

#include <cmath>
#include <fstream>
#include <map>

extern Inform* gmsg;

std::unordered_set<DumpEMFields*> DumpEMFields::dumpsSet_m;

DumpEMFields::DumpEMFields() :
    Action(SIZE, "DUMPEMFIELDS",
           "The \"DUMPEMFIELDS\" statement dumps a field map to a user-defined "
           "field file, for checking that fields are generated correctly. "
           "The fields are written out on a grid in space and time."),
    grid_m(nullptr),
    filename_m("") {

    // would be nice if "steps" could be integer
    itsAttr[FILE_NAME] = Attributes::makeString
        ("FILE_NAME", "Name of the file to which field data is dumped");

    itsAttr[COORDINATE_SYSTEM] = Attributes::makePredefinedString
        ("COORDINATE_SYSTEM", "Choose to use CARTESIAN or CYLINDRICAL coordinates", {"CARTESIAN", "CYLINDRICAL"}, "CARTESIAN");

    itsAttr[X_START] = Attributes::makeReal
        ("X_START", "(Cartesian) Start point in the grid in x [m]");

    itsAttr[DX] = Attributes::makeReal
        ("DX", "(Cartesian) Grid step size in x [m]");

    itsAttr[X_STEPS] = Attributes::makeReal
        ("X_STEPS", "(Cartesian) Number of steps in x");

    itsAttr[Y_START] = Attributes::makeReal
        ("Y_START", "(Cartesian) Start point in the grid in y [m]");

    itsAttr[DY] = Attributes::makeReal
        ("DY", "(Cartesian) Grid step size in y [m]");

    itsAttr[Y_STEPS] = Attributes::makeReal
        ("Y_STEPS", "(Cartesian) Number of steps in y");

    itsAttr[Z_START] = Attributes::makeReal
        ("Z_START", "Start point in the grid in z [m]");

    itsAttr[DZ] = Attributes::makeReal
        ("DZ", "Grid step size in z [m]");

    itsAttr[Z_STEPS] = Attributes::makeReal
        ("Z_STEPS", "Number of steps in z");

    itsAttr[T_START] = Attributes::makeReal
        ("T_START", "Start point in the grid in time [ns]");

    itsAttr[DT] = Attributes::makeReal
        ("DT", "Grid step size in time [ns]");

    itsAttr[T_STEPS] = Attributes::makeReal
        ("T_STEPS", "Number of steps in time");

    itsAttr[R_START] = Attributes::makeReal
        ("R_START", "(Cylindrical) Start point in the grid in radius [m]");

    itsAttr[DR] = Attributes::makeReal
        ("DR", "(Cylindrical) Grid step size in radius [m]");

    itsAttr[R_STEPS] = Attributes::makeReal
        ("R_STEPS", "(Cylindrical) Number of steps in radius");

    itsAttr[PHI_START] = Attributes::makeReal
        ("PHI_START", "(Cylindrical) Start point in the grid in phi [rad]");

    itsAttr[DPHI] = Attributes::makeReal
        ("DPHI", "(Cylindrical) Grid step size in phi [rad]");

    itsAttr[PHI_STEPS] = Attributes::makeReal
        ("PHI_STEPS", "(Cylindrical) Number of steps in phi");

    registerOwnership(AttributeHandler::STATEMENT);
}

DumpEMFields::DumpEMFields(const std::string& name, DumpEMFields* parent):
    Action(name, parent), grid_m(nullptr)
{}

DumpEMFields::~DumpEMFields() {
    delete grid_m;
    dumpsSet_m.erase(this);
}

DumpEMFields* DumpEMFields::clone(const std::string& name) {
    DumpEMFields* dumper = new DumpEMFields(name, this);
    if (grid_m != nullptr) {
        dumper->grid_m = grid_m->clone();
    }
    dumper->filename_m = filename_m;
    dumper->coordinates_m = coordinates_m;
    if (dumpsSet_m.find(this) != dumpsSet_m.end()) {
        dumpsSet_m.insert(dumper);
    }
    return dumper;
}

void DumpEMFields::parseCoordinateSystem() {
    static const std::map<std::string, CoordinateSystem> stringCoordinateSystem_s = {
        {"CARTESIAN",   CoordinateSystem::CARTESIAN},
        {"CYLINDRICAL", CoordinateSystem::CYLINDRICAL}
    };
    coordinates_m = stringCoordinateSystem_s.at(Attributes::getString(itsAttr[COORDINATE_SYSTEM]));
}

void DumpEMFields::execute() {
    buildGrid();
    // the routine for action (OpalParser/OpalParser) calls execute and then
    // deletes 'this'; so we must build a copy that lasts until the field maps
    // are constructed and we are ready for tracking (which is when the field
    // maps are written). Hence the clone call below.
    dumpsSet_m.insert(this->clone(""));
}

void DumpEMFields::buildGrid() {
    std::vector<double> spacing(4);
    std::vector<double> origin(4);
    std::vector<int> gridSize(4);
    parseCoordinateSystem();

    switch (coordinates_m) {
        case CoordinateSystem::CARTESIAN: {
            origin[0] = Attributes::getReal(itsAttr[X_START]);
            spacing[0] = Attributes::getReal(itsAttr[DX]);
            double nx = Attributes::getReal(itsAttr[X_STEPS]);
            checkInt(nx, "X_STEPS");
            gridSize[0] = nx;

            origin[1] = Attributes::getReal(itsAttr[Y_START]);
            spacing[1] = Attributes::getReal(itsAttr[DY]);
            double ny = Attributes::getReal(itsAttr[Y_STEPS]);
            checkInt(ny, "Y_STEPS");
            gridSize[1] = ny;

            break;
        }
        case CoordinateSystem::CYLINDRICAL: {
            origin[0] = Attributes::getReal(itsAttr[R_START]);
            spacing[0] = Attributes::getReal(itsAttr[DR]);
            double nr = Attributes::getReal(itsAttr[R_STEPS]);
            checkInt(nr, "R_STEPS");
            gridSize[0] = nr;

            origin[1] = Attributes::getReal(itsAttr[PHI_START]);
            spacing[1] = Attributes::getReal(itsAttr[DPHI]);
            double nphi = Attributes::getReal(itsAttr[PHI_STEPS]);
            checkInt(nphi, "PHI_STEPS");
            gridSize[1] = nphi;

            break;
        }
    }

    origin[2] = Attributes::getReal(itsAttr[Z_START]);
    spacing[2] = Attributes::getReal(itsAttr[DZ]);
    double nz = Attributes::getReal(itsAttr[Z_STEPS]);
    checkInt(nz, "Z_STEPS");
    gridSize[2] = nz;

    origin[3] = Attributes::getReal(itsAttr[T_START]);
    spacing[3] = Attributes::getReal(itsAttr[DT]);
    double nt = Attributes::getReal(itsAttr[T_STEPS]);
    checkInt(nt, "T_STEPS");
    gridSize[3] = nt;

    if (grid_m != nullptr) {
        delete grid_m;
    }

    grid_m = new interpolation::NDGrid(4, &gridSize[0], &spacing[0], &origin[0]);

    filename_m = Attributes::getString(itsAttr[FILE_NAME]);
}

void DumpEMFields::writeFields(Component* field) {
    typedef std::unordered_set<DumpEMFields*>::iterator dump_iter;
    for (dump_iter it = dumpsSet_m.begin(); it != dumpsSet_m.end(); ++it) {
        (*it)->writeFieldThis(field);
    }
}

void DumpEMFields::checkInt(double real, std::string name, double tolerance) {
    real += tolerance; // prevent rounding error
    if (std::abs(std::floor(real) - real) > 2*tolerance) {
        throw OpalException("DumpEMFields::checkInt",
                            "Value for " + name +
                            " should be an integer but a real value was found");
    }
    if (std::floor(real) < 0.5) {
        throw OpalException("DumpEMFields::checkInt",
                            "Value for " + name + " should be 1 or more");
    }
}

void DumpEMFields::writeHeader(std::ofstream& fout) const {
    fout << grid_m->end().toInteger() << "\n";
    switch (coordinates_m) {
        case CoordinateSystem::CARTESIAN: {
            fout << 1 << "  x [m]\n";
            fout << 2 << "  y [m]\n";
            fout << 3 << "  z [m]\n";
            fout << 4 << "  t [ns]\n";
            fout << 5 << "  Bx [kGauss]\n";
            fout << 6 << "  By [kGauss]\n";
            fout << 7 << "  Bz [kGauss]\n";
            fout << 8 << "  Ex [MV/m]\n";
            fout << 9 << "  Ey [MV/m]\n";
            fout << 10 << " Ez [MV/m]\n";
            break;
        }
        case CoordinateSystem::CYLINDRICAL: {
            fout << 1 << "  r [m]\n";
            fout << 2 << "  phi [deg]\n";
            fout << 3 << "  z [m]\n";
            fout << 4 << "  t [ns]\n";
            fout << 5 << "  Br   [kGauss]\n";
            fout << 6 << "  Bphi [kGauss]\n";
            fout << 7 << "  Bz   [kGauss]\n";
            fout << 8 << "  Er   [MV/m]\n";
            fout << 9 << "  Ephi [MV/m]\n";
            fout << 10 << " Ez   [MV/m]\n";
            break;
        }
    }
    fout << 0 << std::endl;
}

void DumpEMFields::writeFieldLine(Component* field,
                                  const Vector_t& pointIn,
                                  const double& time,
                                  std::ofstream& fout) const {
    Vector_t centroid(0., 0., 0.);
    Vector_t E(0., 0., 0.);
    Vector_t B(0., 0., 0.);
    Vector_t point = pointIn;
    if (coordinates_m == CoordinateSystem::CYLINDRICAL) {
        // pointIn is r, phi, z
        point[0] = std::cos(pointIn[1])*pointIn[0];
        point[1] = std::sin(pointIn[1])*pointIn[0];
    }

    field->apply(point, centroid, time, E, B);
    Vector_t Bout = B;
    Vector_t Eout = E;
    if (coordinates_m == CoordinateSystem::CYLINDRICAL) {
        // pointIn is r, phi, z
        Bout[0] =  B[0]*std::cos(pointIn[1])+B[1]*std::sin(pointIn[1]);
        Bout[1] = -B[0]*std::sin(pointIn[1])+B[1]*std::cos(pointIn[1]);
        Eout[0] =  E[0]*std::cos(pointIn[1])+E[1]*std::sin(pointIn[1]);
        Eout[1] = -E[0]*std::sin(pointIn[1])+E[1]*std::cos(pointIn[1]);
        fout << pointIn[0] << " " << pointIn[1]*Units::rad2deg << " " << pointIn[2] << " " << time << " ";
    } else {
        fout << pointIn[0] << " " << pointIn[1] << " " << pointIn[2] << " " << time << " ";
    }

    fout << Bout[0] << " " << Bout[1] << " " << Bout[2] << " ";
    fout << Eout[0] << " " << Eout[1] << " " << Eout[2] << "\n";
}

void DumpEMFields::writeFieldThis(Component* field) {
    if (grid_m == nullptr) {
        throw OpalException("DumpEMFields::writeFieldThis",
                            "The grid was nullptr; there was a problem with the DumpEMFields initialisation.");
    }
    if (field == nullptr) {
        throw OpalException("DumpEMFields::writeFieldThis",
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

    std::vector<double> point_std(4);
    Vector_t point(0., 0., 0.);
    std::ofstream fout;
    try {
        fout.open(fname.c_str(), std::ofstream::out);
    } catch (std::exception& exc) {
        throw OpalException("DumpEMFields::writeFieldThis",
                            "Failed to open DumpEMFields file " + filename_m);
    }
    if (!fout.good()) {
        throw OpalException("DumpEMFields::writeFieldThis",
                            "Failed to open DumpEMFields file " + filename_m);
    }
    // set precision
    writeHeader(fout);
    for (interpolation::Mesh::Iterator it = grid_m->begin();
         it < grid_m->end();
         ++it) {
        it.getPosition(&point_std[0]);
        for (size_t i = 0; i < 3; ++i) {
            point[i] = point_std[i];
        }
        double time = point_std[3];
        writeFieldLine(field, point, time, fout);
    }
    if (!fout.good()) {
        throw OpalException("DumpEMFields::writeFieldThis",
                            "Something went wrong during writing " + filename_m);
    }
    fout.close();
}

void DumpEMFields::print(std::ostream& os) const {
    os << "* ************* D U M P  E M  F I E L D S ****************************************** " << std::endl;
    os << "* File name: '" << filename_m << "'\n";
    if (coordinates_m == CoordinateSystem::CARTESIAN) {
        os << "* Coordinate system: " << Attributes::getString(itsAttr[COORDINATE_SYSTEM]) << '\n'
           << "* X_START   = " << Attributes::getReal(itsAttr[X_START]) << " [m]\n"
           << "* DX        = " << Attributes::getReal(itsAttr[DX])      << " [m]\n"
           << "* X_STEPS   = " << Attributes::getReal(itsAttr[X_STEPS]) << '\n'
           << "* Y_START   = " << Attributes::getReal(itsAttr[Y_START]) << " [m]\n"
           << "* DY        = " << Attributes::getReal(itsAttr[DY])      << " [m]\n"
           << "* Y_STEPS   = " << Attributes::getReal(itsAttr[Y_STEPS]) << '\n';
    } else if (coordinates_m == CoordinateSystem::CYLINDRICAL) {
        os << "* Coordinate system: " << Attributes::getString(itsAttr[COORDINATE_SYSTEM]) << '\n'
           << "* R_START   = " << Attributes::getReal(itsAttr[R_START])   << " [m]\n"
           << "* DR        = " << Attributes::getReal(itsAttr[DR])        << " [m]\n"
           << "* R_STEPS   = " << Attributes::getReal(itsAttr[R_STEPS])   << '\n'
           << "* PHI_START = " << Attributes::getReal(itsAttr[PHI_START]) << " [rad]\n"
           << "* DPHI      = " << Attributes::getReal(itsAttr[DPHI])      << " [rad]\n"
           << "* PHI_STEPS = " << Attributes::getReal(itsAttr[PHI_STEPS]) << '\n';
    }
    os << "* Z_START   = " << Attributes::getReal(itsAttr[Z_START]) << " [m]\n"
       << "* DZ        = " << Attributes::getReal(itsAttr[DZ])      << " [m]\n"
       << "* Z_STEPS   = " << Attributes::getReal(itsAttr[Z_STEPS]) << '\n'
       << "* T_START   = " << Attributes::getReal(itsAttr[T_START]) << " [ns]\n"
       << "* DT        = " << Attributes::getReal(itsAttr[DT])      << " [ns]\n"
       << "* T_STEPS   = " << Attributes::getReal(itsAttr[T_STEPS]) << '\n';
    os << "* ********************************************************************************** " << std::endl;
}
