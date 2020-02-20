/*
 *  Copyright (c) 2016, Chris Rogers
 *  All rights reserved.
 *  Redistribution and use in source and binary forms, with or without
 *  modification, are permitted provided that the following conditions are met:
 *  1. Redistributions of source code must retain the above copyright notice,
 *     this list of conditions and the following disclaimer.
 *  2. Redistributions in binary form must reproduce the above copyright notice,
 *     this list of conditions and the following disclaimer in the documentation
 *     and/or other materials provided with the distribution.
 *  3. Neither the name of STFC nor the names of its contributors may be used to
 *     endorse or promote products derived from this software without specific
 *     prior written permission.
 *
 *  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 *  AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 *  IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 *  ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
 *  LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 *  CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 *  SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 *  INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 *  CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 *  ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 *  POSSIBILITY OF SUCH DAMAGE.
 */


#include <fstream>

#include "Fields/Interpolation/NDGrid.h"  // classic
#include "AbsBeamline/Component.h"  // classic
#include "Utilities/OpalException.h"
#include "Attributes/Attributes.h"
#include "BasicActions/DumpEMFields.h"

// extern Inform *gmsg;

std::unordered_set<DumpEMFields*> DumpEMFields::dumpsSet_m;

std::string DumpEMFields::dumpemfields_docstring =
std::string("The \"DUMPEMFIELDS\" statement dumps a field map to a user-defined")+
std::string(" field file, for checking that fields are generated correctly.")+
std::string(" The fields are written out on a grid in space and time.");

DumpEMFields::DumpEMFields() :
                    Action(20, "DUMPEMFIELDS", dumpemfields_docstring.c_str()),
                    grid_m(NULL), filename_m("") {
    // would be nice if "steps" could be integer
    itsAttr[0] = Attributes::makeReal
                 ("X_START", "(Cartesian) Start point in the grid in x [m]");
    itsAttr[1] = Attributes::makeReal
                 ("DX", "(Cartesian) Grid step size in x [m]");
    itsAttr[2] = Attributes::makeReal
                 ("X_STEPS", "(Cartesian) Number of steps in x");
    itsAttr[3] = Attributes::makeReal
                 ("Y_START", "(Cartesian) Start point in the grid in y [m]");
    itsAttr[4] = Attributes::makeReal
                 ("DY", "(Cartesian) Grid step size in y [m]");
    itsAttr[5] = Attributes::makeReal
                 ("Y_STEPS", "(Cartesian) Number of steps in y");
    itsAttr[6] = Attributes::makeReal
                 ("Z_START", "Start point in the grid in z [m]");
    itsAttr[7] = Attributes::makeReal
                 ("DZ", "Grid step size in z [m]");
    itsAttr[8] = Attributes::makeReal
                 ("Z_STEPS", "Number of steps in z");
    itsAttr[9] = Attributes::makeReal
                 ("T_START", "Start point in the grid in time [ns]");
    itsAttr[10] = Attributes::makeReal
                 ("DT", "Grid step size in time [ns]");
    itsAttr[11] = Attributes::makeReal
                 ("T_STEPS", "Number of steps in time");
    itsAttr[12] = Attributes::makeString
                 ("FILE_NAME", "Name of the file to which field data is dumped");
    itsAttr[13] = Attributes::makeString("COORDINATE_SYSTEM",
                "Choose to use 'Cartesian' or 'Cylindrical' coordinates");
    itsAttr[14] = Attributes::makeReal
                 ("R_START", "(Cylindrical) Start point in the grid in radius [m]");
    itsAttr[15] = Attributes::makeReal
                 ("DR", "(Cylindrical) Grid step size in radius [m]");
    itsAttr[16] = Attributes::makeReal
                 ("R_STEPS", "(Cylindrical) Number of steps in radius");
    itsAttr[17] = Attributes::makeReal
                 ("PHI_START", "(Cylindrical) Start point in the grid in phi [rad]");
    itsAttr[18] = Attributes::makeReal
                 ("DPHI", "(Cylindrical) Grid step size in phi [rad]");
    itsAttr[19] = Attributes::makeReal
                 ("PHI_STEPS", "(Cylindrical) Number of steps in phi");
}

DumpEMFields::~DumpEMFields() {
    delete grid_m;
    dumpsSet_m.erase(this);
}

DumpEMFields* DumpEMFields::clone(const std::string &/*name*/) {
    DumpEMFields* dumper = new DumpEMFields();
    if (grid_m != NULL) {
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
    if (!itsAttr[13]) {
        coordinates_m = CARTESIAN;
        return;
    }
    std::string coordStr = Attributes::getString(itsAttr[13]);
    for (size_t i = 0; i < coordStr.size(); ++i) {
        coordStr[i] = tolower(coordStr[i]);
    }
    if (coordStr == "cylindrical") {
        coordinates_m = CYLINDRICAL;
    } else if (coordStr == "cartesian") {
        coordinates_m = CARTESIAN;
    } else {
        throw OpalException("DumpEMFields::parseCoordinateSystem",
          std::string("Failed to parse '")+Attributes::getString(itsAttr[13])+
          std::string("'. OPAL expected either 'cylindrical' or 'cartesian'."));
    }
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

    if (coordinates_m == CYLINDRICAL) {
        origin[0] = Attributes::getReal(itsAttr[14]);
        spacing[0] = Attributes::getReal(itsAttr[15]);
        double nr = Attributes::getReal(itsAttr[16]);
        checkInt(nr, "R_STEPS");
        gridSize[0] = nr;

        origin[1] = Attributes::getReal(itsAttr[17]);
        spacing[1] = Attributes::getReal(itsAttr[18]);
        double nphi = Attributes::getReal(itsAttr[19]);
        checkInt(nphi, "PHI_STEPS");
        gridSize[1] = nphi;
    } else {
        origin[0] = Attributes::getReal(itsAttr[0]);
        spacing[0] = Attributes::getReal(itsAttr[1]);
        double nx = Attributes::getReal(itsAttr[2]);
        checkInt(nx, "X_STEPS");
        gridSize[0] = nx;

        origin[1] = Attributes::getReal(itsAttr[3]);
        spacing[1] = Attributes::getReal(itsAttr[4]);
        double ny = Attributes::getReal(itsAttr[5]);
        checkInt(ny, "Y_STEPS");
        gridSize[1] = ny;
    }
    origin[2] = Attributes::getReal(itsAttr[6]);
    spacing[2] = Attributes::getReal(itsAttr[7]);
    double nz = Attributes::getReal(itsAttr[8]);
    checkInt(nz, "Z_STEPS");
    gridSize[2] = nz;

    origin[3] = Attributes::getReal(itsAttr[9]);
    spacing[3] = Attributes::getReal(itsAttr[10]);
    double nt = Attributes::getReal(itsAttr[11]);
    checkInt(nt, "T_STEPS");
    gridSize[3] = nt;

    if (grid_m != NULL) {
        delete grid_m;
    }

    grid_m = new interpolation::NDGrid(4, &gridSize[0], &spacing[0], &origin[0]);

    filename_m = Attributes::getString(itsAttr[12]);
}

void DumpEMFields::writeFields(Component* field) {
    typedef std::unordered_set<DumpEMFields*>::iterator dump_iter;
    for (dump_iter it = dumpsSet_m.begin(); it != dumpsSet_m.end(); ++it) {
        (*it)->writeFieldThis(field);
    }
}

void DumpEMFields::checkInt(double real, std::string name, double tolerance) {
    real += tolerance; // prevent rounding error
    if (fabs(floor(real) - real) > 2*tolerance) {
        throw OpalException("DumpEMFields::checkInt",
                            "Value for "+name+
                            " should be an integer but a real value was found");
    }
    if (floor(real) < 0.5) {
        throw OpalException("DumpEMFields::checkInt",
                            "Value for "+name+" should be 1 or more");
    }
}

void DumpEMFields::writeHeader(std::ofstream& fout) const {
    fout << grid_m->end().toInteger() << "\n";
    if (coordinates_m == CYLINDRICAL) {
        fout << 1 << "  r [mm]\n";
        fout << 2 << "  phi [degree]\n";
    } else {
        fout << 1 << "  x [mm]\n";
        fout << 2 << "  y [mm]\n";
    }
    fout << 3 << "  z [mm]\n";
    fout << 4 << "  t [ns]\n";
    if (coordinates_m == CYLINDRICAL) {
        fout << 5 << "  Br [kGauss]\n";
        fout << 6 << "  Bphi [kGauss]\n";
        fout << 7 << "  Bz [kGauss]\n";
        fout << 8 << "  Er [MV/m]\n";
        fout << 9 << "  Ephi [MV/m]\n";
        fout << 10 << " Ez [MV/m]\n";
    } else {
        fout << 5 << "  Bx [kGauss]\n";
        fout << 6 << "  By [kGauss]\n";
        fout << 7 << "  Bz [kGauss]\n";
        fout << 8 << "  Ex [MV/m]\n";
        fout << 9 << "  Ey [MV/m]\n";
        fout << 10 << " Ez [MV/m]\n";
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
    if (coordinates_m == CYLINDRICAL) {
        // pointIn is r, phi, z 
        point[0] = cos(pointIn[1])*pointIn[0]; 
        point[1] = sin(pointIn[1])*pointIn[0];
    }

    field->apply(point, centroid, time, E, B);
    Vector_t Bout = B;
    Vector_t Eout = E;
    if (coordinates_m == CYLINDRICAL) {
        // pointIn is r, phi, z 
        Bout[0] = B[0]*cos(pointIn[1])+B[1]*sin(pointIn[1]); 
        Bout[1] = -B[0]*sin(pointIn[1])+B[1]*cos(pointIn[1]); 
        Eout[0] = E[0]*cos(pointIn[1])+E[1]*sin(pointIn[1]); 
        Eout[1] = -E[0]*sin(pointIn[1])+E[1]*cos(pointIn[1]);
        fout << pointIn[0] << " " << pointIn[1]*DEGREE << " " << pointIn[2] << " " << time << " ";
    } else {
        fout << pointIn[0] << " " << pointIn[1] << " " << pointIn[2] << " " << time << " ";
    }

    fout << Bout[0] << " " << Bout[1] << " " << Bout[2] << " ";
    fout << Eout[0] << " " << Eout[1] << " " << Eout[2] << "\n";
}

void DumpEMFields::writeFieldThis(Component* field) {
    if (grid_m == NULL) {
        throw OpalException("DumpEMFields::writeFieldThis",
                            "The grid was NULL; there was a problem with the DumpEMFields initialisation.");
    }

    if (field == NULL) {
        throw OpalException("DumpEMFields::writeFieldThis",
                            "The field to be written was NULL.");
    }
    std::vector<double> point_std(4);
    Vector_t point(0., 0., 0.);
    std::ofstream fout;
    try {
        fout.open(filename_m.c_str(), std::ofstream::out);
    } catch (std::exception& exc) {
        throw OpalException("DumpEMFields::writeFieldThis",
                            "Failed to open DumpEMFields file "+filename_m);
    }
    if (!fout.good()) {
        throw OpalException("DumpEMFields::writeFieldThis",
                            "Failed to open DumpEMFields file "+filename_m);
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
                            "Something went wrong during writing "+filename_m);
    }
    fout.close();
}

