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
#ifndef OPAL_BASICACTIONS_DUMPEMFIELDS_HH
#define OPAL_BASICACTIONS_DUMPEMFIELDS_HH

#include "AbsBeamline/Component.h"
#include "AbstractObjects/Action.h"

#include <string>
#include <unordered_set>
#include <vector>

namespace interpolation {
  class NDGrid;
}
class Component;

/** DumpEMFields dumps the dynamically changing fields of a Ring in a user-
 *  defined grid.
 *
 *  The idea is to print out the field map across a 4D grid in space-time for
 *  debugging purposes. The problem is to manage the DumpEMFields object
 *  through three phases of program execution; initial construction, parsing
 *  and then actual field map writing (where we need to somehow let DumpFields
 *  know what the field maps are). So for each DumpEMFields object created, we
 *  store in a set. When the execute() method is called, DumpEMFields builds a
 *  grid using the parsed information.
 *
 *  When the ParallelCyclotronTracker is about to start tracking, it calls
 *  writeFields method which loops over the static set of DumpEMFields and
 *  writes each one. It is not the cleanest implementation, but I can't see a
 *  better way.
 *
 *  The DumpEMFields themselves operate by iterating over a NDGrid object
 *  and looking up the field/writing it out on each grid point.
 *
 */
class DumpEMFields: public Action {

public:
    /// The common attributes of DumpEMFields.
    enum {
        FILE_NAME,
        COORDINATE_SYSTEM,
        X_START,
        DX,
        X_STEPS,
        Y_START,
        DY,
        Y_STEPS,
        Z_START,
        DZ,
        Z_STEPS,
        T_START,
        DT,
        T_STEPS,
        R_START,
        DR,
        R_STEPS,
        PHI_START,
        DPHI,
        PHI_STEPS,
        SIZE
    };

    /** Constructor */
    DumpEMFields();

    /** Constructor */
    DumpEMFields(const std::string& name, DumpEMFields* parent);

    /** Destructor deletes grid_m and if in the dumps set, take it out */
    virtual ~DumpEMFields();

    /** Make a clone (overloadable copy-constructor).
     *    @param name not used
     *  If this is in the dumpsSet_m, so will the clone. Not sure how the
     *  itsAttr stuff works, so this may not get properly copied?
     */
    virtual DumpEMFields* clone(const std::string& name);

    /** Builds the grid but does not write the field map
     *
     *  Builds a grid of points in x-y-z space using the NDGrid algorithm.
     *  Checks that X_STEPS, Y_STEPS, Z_STEPS are integers or throws
     *  OpalException.
     */
    virtual void execute();

    /** Write the fields for all defined DumpEMFields objects
     *    @param field borrowed reference to the Component object that holds the
     *    field map; caller owns the memory.
     *  Iterates over the DumpEMFields in the dumpsSet_m and calls writeFieldThis
     *  on each DumpEMFields. This writes each field map in turn. Format is:
     *  <number of rows>
     *  <column 1> <units>
     *  <column 2> <units>
     *  <column 3> <units>
     *  <column 4> <units>
     *  <column 5> <units>
     *  <column 6> <units>
     *  0
     *  <field map data>
     */
    static void writeFields(Component* field);

    /** Print the attributes of DumpEMFields to standard out */
    void print(std::ostream& os) const;

private:

    enum class CoordinateSystem: unsigned short {
        CARTESIAN,
        CYLINDRICAL
    };

    virtual void writeFieldThis(Component* field);
    virtual void buildGrid();
    void parseCoordinateSystem();
    static void checkInt(double value, std::string name, double tolerance = 1e-9);
    void writeHeader(std::ofstream& fout) const;
    void writeFieldLine(Component* field,
                        const Vector_t& point,
                        const double& time,
                        std::ofstream& fout) const;

    interpolation::NDGrid* grid_m;
    std::string filename_m;

    CoordinateSystem coordinates_m = CoordinateSystem::CARTESIAN;

    static std::unordered_set<DumpEMFields*> dumpsSet_m;

    DumpEMFields(const DumpEMFields& dump);  // disabled
    DumpEMFields& operator=(const DumpEMFields& dump);  // disabled
};

inline std::ostream& operator<<(std::ostream& os, const DumpEMFields& b) {
    b.print(os);
    return os;
}

#endif // ifdef OPAL_DUMPFIELDS_HH
