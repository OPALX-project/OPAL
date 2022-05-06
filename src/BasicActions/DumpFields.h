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
#ifndef OPAL_DUMPFIELDS_HH
#define OPAL_DUMPFIELDS_HH

#include "AbsBeamline/Component.h"
#include "AbstractObjects/Action.h"

#include <string>
#include <unordered_set>

namespace interpolation {
  class ThreeDGrid;
}

class Component;

/** DumpFields dumps the static magnetic field of a Ring in a user-defined grid
 *
 *  The idea is to print out the field map across a grid in space for
 *  debugging purposes. The problem is to manage the DumpFields object through
 *  three phases of program execution; initial construction, parsing and then
 *  actual field map writing (where we need to somehow let DumpFields know what
 *  the field maps are). So for each DumpFields object created, we store in a
 *  set. When the execute() method is called, DumpFields builds a grid using
 *  the parsed information.
 *
 *  When the ParallelCyclotronTracker is about to start tracking, it calls
 *  writeFields method which loops over the static set of DumpFields and writes
 *  each one. It is not the cleanest implementation, but I can't see a better
 *  way.
 *
 *  The DumpFields themselves operate by iterating over a ThreeDGrid object
 *  and looking up the field/writing it out on each grid point.
 *
 *  In order to dump time dependent fields, for example RF, see the
 *  DumpEMFields action.
 */
class DumpFields: public Action {

public:
    /// The common attributes of DumpFields.
    enum {
        FILE_NAME,
        X_START,
        DX,
        X_STEPS,
        Y_START,
        DY,
        Y_STEPS,
        Z_START,
        DZ,
        Z_STEPS,
        SIZE
    };

    /** Constructor */
    DumpFields();

    /** Constructor */
    DumpFields(const std::string& name, DumpFields* parent);

    /** Destructor deletes grid_m and if in the dumps set, take it out */
    virtual ~DumpFields();

    /** Make a clone (overloadable copy-constructor).
     *    @param name not used
     *  If this is in the dumpsSet_m, so will the clone. Not sure how the
     *  itsAttr stuff works, so this may not get properly copied?
     */
    virtual DumpFields* clone(const std::string& name);

    /** Builds the grid but does not write the field map
     *
     *  Builds a grid of points in x-y-z space using the ThreeDGrid algorithm.
     *  Checks that X_STEPS, Y_STEPS, Z_STEPS are integers or throws
     *  OpalException.
     */
    virtual void execute();

    /** Write the fields for all defined DumpFields objects
     *    @param field borrowed reference to the Component object that holds the
     *    field map; caller owns the memory.
     *  Iterates over the DumpFields in the dumpsSet_m and calls writeFieldThis
     *  on each DumpFields. This writes each field map in turn. Format is:
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

    /** Print the attributes of DumpFields to standard out */
    void print(std::ostream& os) const;

private:
    virtual void writeFieldThis(Component* field);
    virtual void buildGrid();
    static void checkInt(double value, std::string name, double tolerance = 1e-9);

    interpolation::ThreeDGrid* grid_m = nullptr;

    std::string filename_m;

    static std::unordered_set<DumpFields*> dumpsSet_m;

    DumpFields(const DumpFields& dump);  // disabled
    DumpFields& operator=(const DumpFields& dump);  // disabled
};

inline std::ostream& operator<<(std::ostream& os, const DumpFields& b) {
    b.print(os);
    return os;
}

#endif // ifdef OPAL_DUMPFIELDS_HH
