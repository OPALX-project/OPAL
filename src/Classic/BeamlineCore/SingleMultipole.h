//
// Class SingleMultipole
//   Representation for single multipoles.
//   Template for representation of single multipoles.
//   Represents all the basic (design) multipole magnets found in an
//   accelerator.  A single multipole has only one multipole component,
//   the pole-number of which cannot be changed (once a quadrupole, always
//   a quadrupole).  This differs from a MultipoleRep object which can
//   have an arbitrary number of multipole components.
//   [P]
//   This template class can be used to instantiate classes like:
//   [UL]
//   [LI]Quadrupole (order = 2),
//   [LI]Sextupole (order = 3),
//   [LI]Octupole (order = 4),
//   [LI]SkewQuadrupole (order = -2),
//   [LI]SkewSextupole (order = -3),
//   [LI]SkewOctupole (order = -4).
//   [/UL]
//   The order and the skew flag are encoded in the template parameter.
//   A positive [b]order[/b] implies a normal multipole,
//   A negative [b]order[/b] implies a skew multipole.
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
#ifndef CLASSIC_SingleMultipole_HH
#define CLASSIC_SingleMultipole_HH

#include "AbsBeamline/Multipole.h"
#include "BeamlineGeometry/StraightGeometry.h"
#include "Channels/IndirectChannel.h"
#include "Fields/BSingleMultipoleField.h"


template <int order>
class SingleMultipole: public Multipole {

public:

    /// Constructor with given name.
    explicit SingleMultipole(const std::string &name);

    SingleMultipole();
    SingleMultipole(const SingleMultipole &);
    virtual ~SingleMultipole();

    /// Get field.
    //  Version for non-constant object.
    virtual BMultipoleField &getField() ;

    /// Get field.
    //  Version for constant object.
    virtual const BMultipoleField &getField() const;

    /// Get geometry.
    //  Version for non-constant object.
    virtual StraightGeometry &getGeometry();

    /// Get geometry.
    //  Version for constant object.
    virtual const StraightGeometry &getGeometry() const;

    /// Get component.
    //  Return the single multipole component in T/m**(n-1).
    virtual double getComponent() const;

    /// Set component.
    //  Assign the single multipole component in T/m**(n-1).
    virtual void setComponent(double Bn);

    /// Return clone.
    //  Return an identical deep copy of the element.
    virtual ElementBase *clone() const;

    /// Construct a read/write channel.
    //  This method constructs a Channel permitting read/write access to
    //  the attribute [b]aKey[/b] and returns it.
    //  If the attribute does not exist, it returns nullptr.
    virtual Channel *getChannel(const std::string &aKey, bool = false);

private:

    // Not implemented.
    void operator=(const SingleMultipole &);

    // A temporary for magnetic field conversion.
    mutable BMultipoleField tempField;

    /// Multipole geometry.
    StraightGeometry geometry;

    /// The single multipole component.
    BSingleMultipoleField<order> field;

    // The type string returned.
    static const std::string type;

    // Attribute access table.
    struct Entry {
        const char *name;
        double(SingleMultipole<order>::*get)() const;
        void (SingleMultipole<order>::*set)(double);
    };

    // The table of attributes.
    static const Entry entries[];
};


// Implementation of template class SingleMultipole
// ------------------------------------------------------------------------

template <int order>
SingleMultipole<order>::SingleMultipole():
    Multipole(),
    geometry(),
    field()
{}


template <int order>
SingleMultipole<order>::SingleMultipole
(const SingleMultipole &multipole):
    Multipole(multipole),
    geometry(multipole.geometry),
    field(multipole.field)
{}


template <int order>
SingleMultipole<order>::SingleMultipole(const std::string &name):
    Multipole(name),
    geometry(),
    field()
{}


template <int order>
SingleMultipole<order>::~SingleMultipole()
{}


template <int order> inline
double SingleMultipole<order>::getComponent() const {
    return field.getComponent();
}

template <int order> inline
void SingleMultipole<order>::setComponent(double value) {
    field.setComponent(value);
}


template <int order>
BMultipoleField &SingleMultipole<order>::getField() {
    tempField = BMultipoleField(field);
    return tempField;
}


template <int order>
const BMultipoleField &SingleMultipole<order>::getField() const {
    tempField = BMultipoleField(field);
    return tempField;
}


template <int order> inline
StraightGeometry &SingleMultipole<order>::getGeometry() {
    return geometry;
}

template <int order> inline
const StraightGeometry &SingleMultipole<order>::getGeometry() const {
    return geometry;
}


template <int order> inline
ElementBase *SingleMultipole<order>::clone() const {
    return new SingleMultipole<order>(*this);
}


template <int order> inline
Channel *SingleMultipole<order>::getChannel(const std::string &aKey, bool) {
    for(const Entry *entry = entries; entry->name != 0; ++entry) {
        if(aKey == entry->name) {
            return new IndirectChannel<SingleMultipole<order> >
                   (*this, entry->get, entry->set);
        }
    }

    return ElementBase::getChannel(aKey);
}

#endif // __SingleMultipole_HH
