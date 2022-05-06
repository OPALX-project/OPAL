//
// Class TravelingWaveRep
//   Representation for a traveling wave.
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
#ifndef CLASSIC_TravelingWaveRep_HH
#define CLASSIC_TravelingWaveRep_HH

#include "AbsBeamline/TravelingWave.h"
#include "BeamlineGeometry/StraightGeometry.h"
#include "Fields/AcceleratingField.h"


class TravelingWaveRep: public TravelingWave {

public:

    /// Constructor with given name.
    explicit TravelingWaveRep(const std::string &name);

    TravelingWaveRep();
    TravelingWaveRep(const TravelingWaveRep &);
    virtual ~TravelingWaveRep();

    /// Return clone.
    //  Return an identical deep copy of the element.
    virtual ElementBase *clone() const;

    /// Construct a read/write channel.
    //  This method constructs a Channel permitting read/write access to
    //  the attribute [b]aKey[/b] and returns it.
    //  If the attribute does not exist, it returns nullptr.
    virtual Channel *getChannel(const std::string &aKey, bool = false);

    /// Get field.
    //  Version for non-constant object.
    virtual AcceleratingField &getField();

    /// Get field.
    //  Version for constant object.
    virtual const AcceleratingField &getField() const;

    /// Get geometry.
    //  Return the element geometry.
    //  Version for non-constant object.
    virtual StraightGeometry &getGeometry();

    /// Get geometry.
    //  Return the element geometry
    //  Version for constant object.
    virtual const StraightGeometry &getGeometry() const;

    /// Get amplitude.
    //  Return the RF amplitude in Volts.
    virtual double getAmplitude() const;

    /// Get frequency.
    //  Return the RF frequency in Hertz.
    virtual double getFrequency() const;

    /// Get phase.
    //  Return the RF phase in radians.
    virtual double getPhase() const;

    /// Set amplitude.
    //  Assign the RF amplitude in Volts.
    virtual void setAmplitude(double V);

    /// Set frequency.
    //  Assign the RF frequency in Hertz.
    virtual void setFrequency(double f);

    /// Set phase.
    //  Assing the RF phase in radians.
    virtual void setPhase(double phi);

    /// Set ignore switch.
    //  If [b]ignore[/b] is true, cavities are ignored for static calculations.
    static void setIgnore(bool ignore = false);

private:

    // Not implemented.
    void operator=(const TravelingWaveRep &);

    /// The cavity's geometry.
    StraightGeometry geometry;

    /// The cavity's field.
    AcceleratingField field;

    /// Cavities are ignored (amplitude = 0) when this switch is set.
    static bool ignoreCavities;
};

#endif // CLASSIC_TravelingWaveRep_HH
