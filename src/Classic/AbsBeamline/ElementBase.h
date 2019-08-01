#ifndef CLASSIC_ElementBase_HH
#define CLASSIC_ElementBase_HH 1

// ------------------------------------------------------------------------
// $RCSfile: ElementBase.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.2 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: ElementBase
//   The very base class for beam line representation objects.  A beam line
//   is modelled as a composite structure having a single root object
//   (the top level beam line), which contains both ``single'' leaf-type
//   elements (Components), as well as sub-lines (composites).
//
// ------------------------------------------------------------------------
// Class category: AbsBeamline
// ------------------------------------------------------------------------
//
// $Date: 2000/12/16 16:26:43 $
// $Author: mad $
//
// ------------------------------------------------------------------------

#include "AbsBeamline/AttributeSet.h"
#include "BeamlineGeometry/Geometry.h"
#include "MemoryManagement/RCObject.h"
#include "Algorithms/Vektor.h"
#include "Algorithms/Quaternion.h"
#include "Algorithms/CoordinateSystemTrafo.h"
#include "Utilities/GeneralClassicException.h"

#include <string>
#include <queue>

class BeamlineVisitor;
class Channel;
class ConstChannel;
class ElementImage;

enum ElemType {
    isDrift,
    isSolenoid,
    isSource,
    isRF,
    isDipole,
    isMultipole,
    isOther
};



class WakeFunction;
class ParticleMatterInteractionHandler;
class BoundaryGeometry;

// Class ElementBase
// ------------------------------------------------------------------------
/// Interface for basic beam line object.
//  This class defines the abstract interface for all objects which can be
//  contained in a beam line. ElementBase forms the base class for two distinct
//  but related hierarchies of objects:
//  [OL]
//  [LI]
//  A set of concrete accelerator element classes, which compose the standard
//  accelerator component library (SACL).
//  [LI]
//  A second hierarchy which parallels the SACL, acting as container or
//  wrapper-like objects.  These latter classes are used to construct
//  beam-lines composed of referenced ``design'' components, together with
//  beam-line position dependent extrinsic state (e. g. errors). These wrapper
//  objects are by default unique.
//  [/OL]
//  Also derived from ElementBase there is a class AlignWrapper, which can
//  be used to store misalignment errors or deliberate displacements.
//  Any element can have only one AlignWrapper and one FieldWrapper.  To be
//  processed in the correct order, the AlignWrapper must be the outermost.
//  [pre]
//    AlignWrapper --> FieldWrapper --> Element
//  [/pre]
//  To ensure this structure, wrapper elements cannot be constructed directly,
//  one should rather call one of:
//  [dl]
//  [dt]makeAlignWrapper()  [dd]make an align wrapper.
//  [dt]makeFieldWrapper()  [dd]make a field wrapper.
//  [dt]makeWrappers()      [dd]make both wrappers.
//  [/dl]
//  An existing wrapper can be removed by
//  [dl]
//  [dt]removeAlignWrapper()[dd]remove align wrapper.
//  [dt]removeFieldWrapper()[dd]remove field wrapper.
//  [/dl]
//  Instances of the concrete classes for single elements are by default
//  sharable.  Instances of beam lines, wrappers and integrators are by
//  default non-sharable, but they may be made sharable by a call to
//  [b]makeSharable()[/b].
//  [p]
//  An ElementBase object can return two lengths, which may be different:
//  [OL]
//  [LI]
//  The arc length along the geometry.
//  [LI]
//  The design length, often measured along a straight line.
//  [/OL]
//  Class ElementBase contains a map of name versus value for user-defined
//  attributes (see file AbsBeamline/AttributeSet.hh).  The map is primarily
//  intended for processes that require algorithm-specific data in the
//  accelerator model.
//  [P]
//  The class ElementBase has as base class the abstract class RCObject.
//  Virtual derivation is used to make multiple inheritance possible for
//  derived concrete classes. ElementBase implements three copy modes:
//  [OL]
//  [LI]
//  Copy by reference: Call RCObject::addReference() and use [b]this[/b].
//  [LI]
//  Copy structure: use ElementBase::copyStructure().
//  During copying of the structure, all sharable items are re-used, while
//  all non-sharable ones are cloned.
//  [LI]
//  Copy by cloning: use ElementBase::clone().
//  This returns a full deep copy.
//  [/OL]


class ElementBase: public RCObject {

public:

    /// Constructor with given name.
    explicit ElementBase(const std::string &name);

    ElementBase();
    ElementBase(const ElementBase &);
    virtual ~ElementBase();

    /// Get element name.
    virtual const std::string &getName() const;

    /// Set element name.
    virtual void setName(const std::string &name);

    enum ApertureType {RECTANGULAR
                     , ELLIPTICAL
                     , CONIC_RECTANGULAR
                     , CONIC_ELLIPTICAL
    };

    enum ElementType {ALIGNWRAPPER
                    , BEAMBEAM
                    , BEAMBEAM3D
                    , BEAMLINE
                    , CCOLLIMATOR
                    , CORRECTOR
                    , CORRECTORWRAPPER
                    , CYCLOTRON
                    , CYCLOTRONWRAPPER
                    , CYCLOTRONVALLEY
                    , DEGRADER
                    , DIAGNOSTIC
                    , DRIFT
                    , FLEXIBLECOLLIMATOR
                    , INTEGRATOR
                    , LAMBERTSON
                    , MARKER
                    , MONITOR
                    , MPSPLITINTEGRATOR
                    , MULTIPOLE
                    , MULTIPOLET
                    , MULTIPOLEWRAPPER
                    , OFFSET
                    , PARALLELPLATE
                    , PATCH
                    , PROBE
                    , RBEND
                    , RBEND3D
                    , RBENDWRAPPER
                    , RFCAVITY
                    , RFQUADRUPOLE
                    , RING
                    , SBEND3D
                    , SBEND
                    , SBENDWRAPPER
                    , SEPARATOR
                    , SEPTUM
                    , SOLENOID
                    , SOURCE
                    , STRIPPER
                    , TRAVELINGWAVE
                    , VARIABLERFCAVITY
                    , ANY};

    /// Get element type std::string.
    virtual ElementType getType() const = 0;

    std::string getTypeString() const;
    static std::string getTypeString(ElementType type);

    /// Get geometry.
    //  Return the element geometry.
    //  Version for non-constant object.
    virtual BGeometryBase  &getGeometry() = 0;

    /// Get geometry.
    //  Return the element geometry
    //  Version for constant object.
    virtual const BGeometryBase  &getGeometry() const = 0;

    /// Get arc length.
    //  Return the entire arc length measured along the design orbit
    virtual double getArcLength() const;

    /// Get design length.
    //  Return the design length defined by the geometry.
    //  This may be the arc length or the straight length.
    virtual double getElementLength() const;

    /// Set design length.
    //  Set the design length defined by the geometry.
    //  This may be the arc length or the straight length.
    virtual void setElementLength(double length);

    virtual void getElementDimensions(double &begin,
                                      double &end) const {
        begin = 0.0;
        end = getElementLength();
    }

    /// Get origin position.
    //  Return the arc length from the entrance to the origin of the element
    //  (origin >= 0)
    virtual double getOrigin() const;

    /// Get entrance position.
    //  Return the arc length from the origin to the entrance of the element
    //  (entrance <= 0)
    virtual double getEntrance() const;

    /// Get exit position.
    //  Return the arc length from the origin to the exit of the element
    //  (exit >= 0)
    virtual double getExit() const;

    /// Get transform.
    //  Return the transform of the local coordinate system from the
    //  position [b]fromS[/b] to the position [b]toS[/b].
    virtual Euclid3D getTransform(double fromS, double toS) const;

    /// Get transform.
    //  Equivalent to getTransform(0.0, s).
    //  Return the transform of the local coordinate system from the
    //  origin and [b]s[/b].
    virtual Euclid3D getTransform(double s) const;

    /// Get transform.
    //  Equivalent to getTransform(getEntrance(), getExit()).
    //  Return the transform of the local coordinate system from the
    //  entrance to the exit of the element.
    virtual Euclid3D getTotalTransform() const;

    /// Get transform.
    //  Equivalent to getTransform(0.0, getEntrance()).
    //  Return the transform of the local coordinate system from the
    //  origin to the entrance of the element.
    virtual Euclid3D getEntranceFrame() const;

    /// Get transform.
    //  Equivalent to getTransform(0.0, getExit()).
    //  Return the transform of the local coordinate system from the
    //  origin to the exit of the element.
    virtual Euclid3D getExitFrame() const;

    /// Get patch.
    //  Returns the entrance patch (transformation) which is used to transform
    //  the global geometry to the local geometry for a misaligned element
    //  at its entrance. The default behaviour returns identity transformation.
    //  This function should be overridden by derived concrete classes which
    //  model complex geometries.
    virtual Euclid3D getEntrancePatch() const;

    /// Get patch.
    //  Returns the entrance patch (transformation) which is used to transform
    //  the local geometry to the global geometry for a misaligned element
    //  at its exit. The default behaviour returns identity transformation.
    //  This function should be overridden by derived concrete classes which
    //  model complex geometries.
    virtual Euclid3D getExitPatch() const;

    /// Get attribute value.
    //  If the attribute does not exist, return zero.
    virtual double getAttribute(const std::string &aKey) const;

    /// Test for existence of an attribute.
    //  If the attribute exists, return true, otherwise false.
    virtual bool hasAttribute(const std::string &aKey) const;

    /// Remove an existing attribute.
    virtual void removeAttribute(const std::string &aKey);

    /// Set value of an attribute.
    virtual void setAttribute(const std::string &aKey, double val);

    /// Construct a read/write channel.
    //  This method constructs a Channel permitting read/write access to
    //  the attribute [b]aKey[/b] and returns it.
    //  If the attribute does not exist, it returns NULL.
    virtual Channel *getChannel(const std::string &aKey, bool create = false);

    /// Construct a read-only channel.
    //  This method constructs a Channel permitting read-only access to
    //  the attribute [b]aKey[/b] and returns it.
    //  If the attribute does not exist, it returns NULL.
    virtual const ConstChannel *getConstChannel(const std::string &aKey) const;

    /// Construct an image.
    //  Return the image of the element, containing the name and type std::string
    //  of the element, and a copy of the user-defined attributes.
    virtual ElementImage *getImage() const;

    /// Apply visitor.
    //  This method must be overridden by derived classes. It should call the
    //  method of the visitor corresponding to the element class.
    //  If any error occurs, this method throws an exception.
    virtual void accept(BeamlineVisitor &visitor) const = 0;


    /// Return clone.
    //  Return an identical deep copy of the element.
    virtual ElementBase *clone() const = 0;

    /// Make a structural copy.
    //  Return a fresh copy of any beam line structure is made,
    //  but sharable elements remain shared.
    virtual ElementBase *copyStructure();

    /// Test if the element can be shared.
    bool isSharable() const;

    /// Set sharable flag.
    //  The whole structure depending on [b]this[/b] is marked as sharable.
    //  After this call a [b]copyStructure()[/b] call reuses the element.
    virtual void makeSharable();

    /// Allow misalignment.
    //  Build an AlignWrapper pointing to the element and return a pointer to
    //  that wrapper.  If the element cannot be misaligned, or already has an
    //  AlignWrapper, return a pointer to the element.
    //  Wrappers are non-sharable, unless otherwise defined.
    virtual ElementBase *makeAlignWrapper();

    /// Allow field errors.
    //  Build a FieldWrapper pointing to the element and return a pointer to
    //  that wrapper.  If the element cannot have field errors, or already has
    //  a FieldWrapper, return a pointer to the element.
    //  Wrappers are non-sharable, unless otherwise defined.
    virtual ElementBase *makeFieldWrapper();

    /// Allow errors.
    //  Equivalent to the calls
    //  [pre]
    //    makeFieldWrapper()->makeAlignWrapper()
    //  [/pre].
    //  Wrappers are non-sharable, unless otherwise defined.
    virtual ElementBase *makeWrappers();

    /// Remove align wrapper.
    //  Remove the align wrapper.
    virtual ElementBase *removeAlignWrapper();

    /// Remove align wrapper.
    //  Remove the align wrapper.
    virtual const ElementBase *removeAlignWrapper() const ;

    /// Remove field wrapper.
    //  Remove the field wrapper.
    virtual ElementBase *removeFieldWrapper();

    /// Remove field wrapper.
    //  Remove the field wrapper for constant object.
    virtual const ElementBase *removeFieldWrapper() const;

    /// Return the design element.
    //  Return [b]this[/b], if the element is not wrapped.
    virtual ElementBase *removeWrappers();

    /// Return the design element.
    //  Return [b]this[/b], if the element is not wrapped.
    virtual const  ElementBase *removeWrappers() const ;

    /// Update element.
    //  This method stores all attributes contained in the AttributeSet to
    //  "*this".  The return value [b]true[/b] indicates success.
    bool update(const AttributeSet &);

    ///@{ Access to ELEMEDGE attribute
    void setElementPosition(double elemedge);
    double getElementPosition() const;
    bool isElementPositionSet() const;
    ///@}
    /// attach a boundary geometry field to the element
    virtual void setBoundaryGeometry(BoundaryGeometry *geo);

    /// return the attached boundary geometrt object if there is any
    virtual BoundaryGeometry *getBoundaryGeometry() const;

    virtual bool hasBoundaryGeometry() const;


    /// attach a wake field to the element
    virtual void setWake(WakeFunction *wf);

    /// return the attached wake object if there is any
    virtual WakeFunction *getWake() const;

    virtual bool hasWake() const;

    virtual void setParticleMatterInteraction(ParticleMatterInteractionHandler *spys);

    virtual ParticleMatterInteractionHandler *getParticleMatterInteraction() const;

    virtual bool hasParticleMatterInteraction() const;

    /// returns element type as enumeration needed in the envelope tracker
    ElemType getElType() const;

    /// set the element type as enumeration needed in the envelope tracker
    void setElType(ElemType elt);

    void setCSTrafoGlobal2Local(const CoordinateSystemTrafo &ori);
    CoordinateSystemTrafo getCSTrafoGlobal2Local() const;
    void releasePosition();
    void fixPosition();
    bool isPositioned() const;

    virtual CoordinateSystemTrafo getEdgeToBegin() const;
    virtual CoordinateSystemTrafo getEdgeToEnd() const;

    void setAperture(const ApertureType& type, const std::vector<double> &args);
    std::pair<ElementBase::ApertureType, std::vector<double> > getAperture() const;

    virtual bool isInside(const Vector_t &r) const;

    void setMisalignment(double x, double y, double s);
    void setMisalignment(const CoordinateSystemTrafo &cst);

    void getMisalignment(double &x, double &y, double &s) const;
    CoordinateSystemTrafo getMisalignment() const;

    void setActionRange(const std::queue<std::pair<double, double> > &range);
    void setCurrentSCoordinate(double s);

    /// Set rotation about z axis in bend frame.
    void setRotationAboutZ(double rotation);
    double getRotationAboutZ() const;

    void checkLengthResolution(double dt,
                               const Vector_t &P) const;
protected:
    bool isInsideTransverse(const Vector_t &r, double f = 1) const;

    // Sharable flag.
    // If this flag is true, the element is always shared.
    mutable bool shareFlag;

    CoordinateSystemTrafo csTrafoGlobal2Local_m;
    CoordinateSystemTrafo misalignment_m;

    std::pair<ApertureType, std::vector<double> > aperture_m;

    double elementEdge_m;

    double rotationZAxis_m;

private:

    // Not implemented.
    void operator=(const ElementBase &);

    // The element's name
    std::string elementID;

    // The user-defined set of attributes.
    AttributeSet userAttribs;

    WakeFunction *wake_m;

    BoundaryGeometry *bgeometry_m;

    ParticleMatterInteractionHandler *parmatint_m;

    ElemType elType_m;

    bool positionIsFixed;
    ///@{ ELEMEDGE attribute
    double elementPosition_m;
    bool elemedgeSet_m;
    ///@}
    std::queue<std::pair<double, double> > actionRange_m;
};


// Inline functions.
// ------------------------------------------------------------------------

inline
double ElementBase::getArcLength() const
{ return getGeometry().getArcLength(); }

inline
double ElementBase::getElementLength() const
{ return getGeometry().getElementLength(); }

inline
void ElementBase::setElementLength(double length)
{ getGeometry().setElementLength(length); }

inline
double ElementBase::getOrigin() const
{ return getGeometry().getOrigin(); }

inline
double ElementBase::getEntrance() const
{ return getGeometry().getEntrance(); }

inline
double ElementBase::getExit() const
{ return getGeometry().getExit(); }

inline
Euclid3D ElementBase::getTransform(double fromS, double toS) const
{ return getGeometry().getTransform(fromS, toS); }

inline
Euclid3D ElementBase::getTotalTransform() const
{ return getGeometry().getTotalTransform(); }

inline
Euclid3D ElementBase::getTransform(double s) const
{ return getGeometry().getTransform(s); }

inline
Euclid3D ElementBase::getEntranceFrame() const
{ return getGeometry().getEntranceFrame(); }

inline
Euclid3D ElementBase::getExitFrame() const
{ return getGeometry().getExitFrame(); }

inline
Euclid3D ElementBase::getEntrancePatch() const
{ return getGeometry().getEntrancePatch(); }

inline
Euclid3D ElementBase::getExitPatch() const
{ return getGeometry().getExitPatch(); }

inline
bool ElementBase::isSharable() const
{ return shareFlag; }

inline
WakeFunction *ElementBase::getWake() const
{ return wake_m; }

inline
bool ElementBase::hasWake() const
{ return wake_m != NULL; }

inline
BoundaryGeometry *ElementBase::getBoundaryGeometry() const
{ return bgeometry_m; }

inline
bool ElementBase::hasBoundaryGeometry() const
{ return bgeometry_m != NULL; }

inline
ParticleMatterInteractionHandler *ElementBase::getParticleMatterInteraction() const
{ return parmatint_m; }

inline
bool ElementBase::hasParticleMatterInteraction() const
{ return parmatint_m != NULL; }

inline
ElemType ElementBase::getElType() const
{ return elType_m;}

inline
void ElementBase::setElType(ElemType elt)
{ elType_m = elt;}

inline
void ElementBase::setCSTrafoGlobal2Local(const CoordinateSystemTrafo &trafo)
{
    if (positionIsFixed) return;

    csTrafoGlobal2Local_m = trafo;
}

inline
CoordinateSystemTrafo ElementBase::getCSTrafoGlobal2Local() const
{ return csTrafoGlobal2Local_m; }

inline
CoordinateSystemTrafo ElementBase::getEdgeToBegin() const
{
    CoordinateSystemTrafo ret(Vector_t(0, 0, 0),
                              Quaternion(1, 0, 0, 0));

    return ret;
}

inline
CoordinateSystemTrafo ElementBase::getEdgeToEnd() const
{
    CoordinateSystemTrafo ret(Vector_t(0, 0, getElementLength()),
                              Quaternion(1, 0, 0, 0));

    return ret;
}

inline
void ElementBase::setAperture(const ApertureType& type, const std::vector<double> &args)
{
    aperture_m.first = type;
    aperture_m.second = args;
}

inline
std::pair<ElementBase::ApertureType, std::vector<double> > ElementBase::getAperture() const
{
    return aperture_m;
}

inline
bool ElementBase::isInside(const Vector_t &r) const
{
    const double length = getElementLength();
    return isInsideTransverse(r, r(2) / length * aperture_m.second[2]) && r(2) >= 0.0 && r(2) < length;
}

inline
bool ElementBase::isInsideTransverse(const Vector_t &r, double f) const
{
    switch(aperture_m.first) {
    case RECTANGULAR:
        return (std::abs(r[0]) < aperture_m.second[0] && std::abs(r[1]) < aperture_m.second[1]);
    case ELLIPTICAL:
        return (std::pow(r[0] / aperture_m.second[0], 2) + std::pow(r[1] / aperture_m.second[1], 2) < 1.0);
    case CONIC_RECTANGULAR:
        return (std::abs(r[0]) < f * aperture_m.second[0] && std::abs(r[1]) < f * aperture_m.second[1]);
    case CONIC_ELLIPTICAL:
        return (std::pow(r[0] / (f * aperture_m.second[0]), 2) + std::pow(r[1] / (f * aperture_m.second[1]), 2) < 1.0);
    default:
        return false;
    }
}

inline
void ElementBase::setMisalignment(const CoordinateSystemTrafo &cst) {
    misalignment_m = cst;
}

inline
CoordinateSystemTrafo ElementBase::getMisalignment() const {
    return misalignment_m;
}

inline
void ElementBase::releasePosition() {
    positionIsFixed = false;
}

inline
void ElementBase::fixPosition() {
    positionIsFixed = true;
}

inline
bool ElementBase::isPositioned() const {
    return positionIsFixed;
}

inline
void ElementBase::setActionRange(const std::queue<std::pair<double, double> > &range) {
    actionRange_m = range;

    if (actionRange_m.size() > 0)
        elementEdge_m = actionRange_m.front().first;
}

inline
void ElementBase::setRotationAboutZ(double rotation) {
    rotationZAxis_m = rotation;
}

inline
double ElementBase::getRotationAboutZ() const {
    return rotationZAxis_m;
}

inline
std::string ElementBase::getTypeString() const
{ return getTypeString(getType());}

inline
void ElementBase::setElementPosition(double elemedge) {
    elementPosition_m = elemedge;
    elemedgeSet_m = true;
}

inline
double ElementBase::getElementPosition() const {
    if (elemedgeSet_m)
        return elementPosition_m;

    throw GeneralClassicException("ElementBase::getElementPosition()",
                                  std::string("ELEMEDGE for \"") + getName() + "\" not set");
}

inline
bool ElementBase::isElementPositionSet() const
{ return elemedgeSet_m; }

#endif // CLASSIC_ElementBase_HH