//
// Class UndulatorRep
//   Defines a concrete undulator/wiggler representation.
//
// Copyright (c) 2020, Arnau Albà, Paul Scherrer Institut, Villigen PSI, Switzerland
// All rights reserved.
//
// Implemented as part of the MSc thesis
// "Start-to-End Modelling of the AWA Micro-Bunched Electron Cooling POP-Experiment"
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
#ifndef CLASSIC_UndulatorRep_HH
#define CLASSIC_UndulatorRep_HH

#include "AbsBeamline/Undulator.h"

#include "BeamlineGeometry/StraightGeometry.h"
#include "Fields/NullField.h"

class UndulatorRep : public Undulator {
public:
    /// Constructor with given name.
    explicit UndulatorRep(const std::string& name);

    UndulatorRep();
    UndulatorRep(const UndulatorRep&);
    virtual ~UndulatorRep();

    /// Return clone.
    //  Return an identical deep copy of the element.
    virtual ElementBase* clone() const;

    /// Construct a read/write channel.
    //  This method constructs a Channel permitting read/write access to
    //  the attribute [b]aKey[/b] and returns it.
    //  If the attribute does not exist, it returns nullptr.
    virtual Channel* getChannel(const std::string& aKey, bool = false);

    /// Get field.
    //  Version for non-constant object.
    virtual NullField& getField();

    /// Get field.
    //  Version for constant object.
    virtual const NullField& getField() const;

    /// Get geometry.
    //  Version for non-constant object.
    virtual StraightGeometry& getGeometry();

    /// Get geometry.
    //  Version for constant object.
    virtual const StraightGeometry& getGeometry() const;

private:
    // Not implemented.
    void operator=(const UndulatorRep&);

    /// The zero magnetic field.
    NullField field;

    /// The geometry.
    StraightGeometry geometry;
};

#endif  // CLASSIC_UndulatorRep_HH
