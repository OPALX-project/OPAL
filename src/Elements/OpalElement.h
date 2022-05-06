//
// Class OpalElement
//   Base class for all beam line elements.
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
#ifndef OPAL_OpalElement_HH
#define OPAL_OpalElement_HH

#include "AbstractObjects/Element.h"
#include "MemoryManagement/OwnPtr.h"
#include <map>
#include <string>

class Statement;

class OpalElement: public Element {

public:

    /// The common attributes for all elements.
    enum {
        TYPE,           // The design type.
        APERT,          // The aperture data.
        LENGTH,         // The element length.
        ELEMEDGE,       // The position of the element (in path length)
        WAKEF,          // The wake function to be used
        PARTICLEMATTERINTERACTION, // The particle mater interaction handler to be used
        ORIGIN,         // The location of the element in floor coordinates
        ORIENTATION,    // The orientation of the element (Tait Bryan angles)
        X,              // The x-coordinate of the location of the element in floor coordinates
        Y,              // The y-coordinate of the location of the element in floor coordinates
        Z,              // The z-coordinate of the location of the element in floor coordinates
        THETA,          // The rotation about the y-axis
        PHI,            // The rotation about the x-axis
        PSI,            // The rotation about the z-axis
        DX,             // Misalignment in x (local coordinate system)
        DY,             // Misalignment in y (local coordinate system)
        DZ,             // Misalignment in z (local coordinate system)
        DTHETA,         // The rotation around y axis in rad.
        DPHI,           // The rotation around x axis in rad.
        DPSI,           // The rotation around s axis in rad.
        OUTFN,          // Output filename
        DELETEONTRANSVERSEEXIT, // Flag whether particles should be deleted if exit transversally
        COMMON
    };

    virtual ~OpalElement();

    /// Return element length.
    virtual double getLength() const;

    /// Return the element's type name.
    const std::string getTypeName() const;

    //return the element aperture vector
    std::pair<ApertureType, std::vector<double> > getApert() const;

    /// Return the element's type name.
    const std::string getWakeF() const;

    const std::string getParticleMatterInteraction() const;

    const std::string getWMaterial() const;

    const std::string getWakeGeom() const;

    std::vector<double> getWakeParam() const;

    const std::string getWakeConductivity() const;

    /// Parse the element.
    //  This special version for elements handles unknown attributes by
    //  appending them to the attribute list.
    virtual void parse(Statement&);

    /// Print the object.
    //  This special version handles special printing in OPAL-8 format.
    virtual void print(std::ostream&) const;

    /// Update the embedded CLASSIC element.
    virtual void update();

    /// Transmit the ``unknown'' (not known to OPAL) attributes to CLASSIC.
    virtual void updateUnknown(ElementBase*);

protected:

    /// Exemplar constructor.
    OpalElement(int size, const char* name, const char* help);

    /// Clone constructor.
    OpalElement(const std::string& name, OpalElement* parent);

    /// Print multipole components in OPAL-8 format.
    //  This function is accessible to all multipole-like elements
    //  (RBend, SBend, Quadrupole, Sextupole, Octupole, Multipole).
    static void printMultipoleStrength(std::ostream& os,
                                       int order,
                                       int& len,
                                       const std::string& sName,
                                       const std::string& tName,
                                       const Attribute& length,
                                       const Attribute& vNorm,
                                       const Attribute& vSkew);

    /// Print an attribute with a OPAL-8 name (as an expression).
    static void printAttribute(std::ostream& os,
                               const std::string& name,
                               const std::string& image,
                               int& len);

    /// Print an attribute with a OPAL-8 name (as a constant).
    static void printAttribute(std::ostream& os,
                               const std::string& name,
                               double value,
                               int& len);

    void registerOwnership() const;

private:

    // Not implemented.
    OpalElement();
    void operator=(const OpalElement&);

    // The original size of the attribute list.
    int itsSize;
};

#endif // OPAL_OpalElement_HH
