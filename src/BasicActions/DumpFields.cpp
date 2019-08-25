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

#include "Fields/Interpolation/ThreeDGrid.h"  // classic
#include "AbsBeamline/Component.h"  // classic
#include "Utilities/OpalException.h"
#include "Attributes/Attributes.h"
#include "BasicActions/DumpFields.h"

std::unordered_set<DumpFields*> DumpFields::dumpsSet_m;

std::string DumpFields::dumpfields_docstring =
std::string("The \"DUMPFIELDS\" statement dumps a field map to a user-defined")+
std::string(" field file, for checking that fields are read in correctly")+
std::string(" from disk. The fields are written out on a Cartesian grid.");

DumpFields::DumpFields() :
                    Action(10, "DUMPFIELDS", dumpfields_docstring.c_str()) {
    // would be nice if "steps" could be integer
    itsAttr[0] = Attributes::makeReal
                 ("X_START", "Start point in the grid in x [m]");
    itsAttr[1] = Attributes::makeReal
                 ("DX", "Grid step size in x [m]");
    itsAttr[2] = Attributes::makeReal
                 ("X_STEPS", "Number of steps in x");
    itsAttr[3] = Attributes::makeReal
                 ("Y_START", "Start point in the grid in y [m]");
    itsAttr[4] = Attributes::makeReal
                 ("DY", "Grid step size in y [m]");
    itsAttr[5] = Attributes::makeReal
                 ("Y_STEPS", "Number of steps in y");
    itsAttr[6] = Attributes::makeReal
                 ("Z_START", "Start point in the grid in z [m]");
    itsAttr[7] = Attributes::makeReal
                 ("DZ", "Grid step size in z [m]");
    itsAttr[8] = Attributes::makeReal
                 ("Z_STEPS", "Number of steps in z");
    itsAttr[9] = Attributes::makeString
                 ("FILE_NAME", "Name of the file to which field data is dumped");

    registerOwnership(AttributeHandler::STATEMENT);
}

DumpFields::DumpFields(const std::string &name, DumpFields *parent):
    Action(name, parent)
{}

DumpFields::~DumpFields() {
    delete grid_m;
    dumpsSet_m.erase(this);
}

DumpFields* DumpFields::clone(const std::string &name) {
    DumpFields* dumper = new DumpFields(name, this);
    if (grid_m != NULL) {
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
    double x0 = Attributes::getReal(itsAttr[0]);
    double dx = Attributes::getReal(itsAttr[1]);
    double nx = Attributes::getReal(itsAttr[2]);

    double y0 = Attributes::getReal(itsAttr[3]);
    double dy = Attributes::getReal(itsAttr[4]);
    double ny = Attributes::getReal(itsAttr[5]);

    double z0 = Attributes::getReal(itsAttr[6]);
    double dz = Attributes::getReal(itsAttr[7]);
    double nz = Attributes::getReal(itsAttr[8]);

    checkInt(nx, "X_STEPS");
    checkInt(ny, "Y_STEPS");
    checkInt(nz, "Z_STEPS");
    delete grid_m;

    grid_m = new interpolation::ThreeDGrid(dx, dy, dz,
                                           x0, y0, z0,
                                           nx, ny, nz);

    filename_m = Attributes::getString(itsAttr[9]);
}

void DumpFields::writeFields(Component* field) {
    typedef std::unordered_set<DumpFields*>::iterator dump_iter;
    for (dump_iter it = dumpsSet_m.begin(); it != dumpsSet_m.end(); ++it) {
        (*it)->writeFieldThis(field);
    }
}

void DumpFields::checkInt(double real, std::string name, double tolerance) {
    if (fabs(floor(real) - real) > tolerance) {
        throw OpalException("DumpFields::checkInt",
                            "Value for "+name+
                            " should be an integer but a real value was found");
    }
    if (floor(real) < 0.5) {
        throw OpalException("DumpFields::checkInt",
                            "Value for "+name+" should be 1 or more");
    }
}

void DumpFields::writeFieldThis(Component* field) {
    if (grid_m == NULL) {
        throw OpalException("DumpFields::writeFieldThis",
                            "The grid was NULL; there was a problem with the DumpFields initialisation.");
    }
    if (field == NULL) {
        throw OpalException("DumpFields::writeFieldThis",
                            "The field to be written was NULL.");
    }
    double time = 0.;
    Vector_t point(0., 0., 0.);
    Vector_t centroid(0., 0., 0.);
    std::ofstream fout(filename_m.c_str(), std::ofstream::out);
    if (!fout.good()) {
        throw OpalException("DumpFields::writeFieldThis",
                            "Failed to open DumpFields file "+filename_m);
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
                            "Something went wrong during writing "+filename_m);
    }
    fout.close();
}
