//
// Class OpalParallelPlate
//   The ParallelPlate element.
//
// Copyright (c) 200x - 2020, Paul Scherrer Institut, Villigen PSI, Switzerland
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
#ifndef OPAL_OpalParallelPlate_HH
#define OPAL_OpalParallelPlate_HH

#include "Elements/OpalElement.h"

class BoundaryGeometry;

class OpalParallelPlate: public OpalElement {

public:

    /// The attributes of class OpalParallelPlate.
    enum {
        VOLT = COMMON,  // The peak voltage.
        GEOMETRY,       // geometry of boundary
        FREQ,           // The RF frequency.
        LAG,            // The phase lag.

        PLENGTH,           //distance between two plates or length in s direction
        SIZE
    };

    /// Exemplar constructor.
    OpalParallelPlate();

    virtual ~OpalParallelPlate();

    /// Make clone.
    virtual OpalParallelPlate *clone(const std::string &name);

    /// Fill in all registered attributes.
    virtual void fillRegisteredAttributes(const ElementBase &);

    /// Update the embedded CLASSIC cavity.
    virtual void update();

private:

    // Not implemented.
    OpalParallelPlate(const OpalParallelPlate &);
    void operator=(const OpalParallelPlate &);

    // Clone constructor.
    OpalParallelPlate(const std::string &name, OpalParallelPlate *parent);



    BoundaryGeometry *obgeo_m;


};

#endif // OPAL_OpalParallelPlate_HH
