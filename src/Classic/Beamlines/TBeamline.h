#ifndef CLASSIC_TBeamline_HH
#define CLASSIC_TBeamline_HH

// ------------------------------------------------------------------------
// $RCSfile: TBeamline.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.3.2.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Template class: TBeamline<class T>
//
// ------------------------------------------------------------------------
// Class category: Beamlines
// ------------------------------------------------------------------------
//
// $Date: 2004/11/18 22:18:06 $
// $Author: jsberg $
//
// ------------------------------------------------------------------------

#include "Beamlines/Beamline.h"
#include "AbsBeamline/BeamlineVisitor.h"
#include "Beamlines/BeamlineGeometry.h"
#include "Utilities/CLRangeError.h"
#include "Algorithms/Vektor.h"
#include "Algorithms/Quaternion.h"
#include <algorithm>
#include <list>


// Template class TBeamline
// ------------------------------------------------------------------------
/// Template class for beam lines.
//  Instantiation with different T types allows attachment of additional
//  data to each position in the line.

template <class T>
class TBeamline: public Beamline, public std::list<T> {

public:

    /// Default constructor.
    TBeamline();

    /// Constructor with given name.
    explicit TBeamline(const std::string &name);

    TBeamline(const TBeamline<T> &right);
    virtual ~TBeamline();

    /// Apply BeamlineVisitor to this line.
    virtual void accept(BeamlineVisitor &) const;

    /// Apply visitor to all elements of the line.
    //  If the parameter [b]r2l[/b] is true, the line is traversed from
    //  right (s=C) to left (s=0).
    //  If any error occurs, this method may throw an exception.
    virtual void iterate(BeamlineVisitor &, bool r2l) const;

    /// Make clone.
    virtual TBeamline<T> *clone() const;

    /// Make structure copy.
    virtual TBeamline<T> *copyStructure();

    /// Set sharable flag.
    //  The whole beamline and the elements depending on [b]this[/b] are
    //  marked as sharable.  After this call a [b]copyStructure()[/b] call
    //  reuses the element.
    virtual void makeSharable();

    /// Get geometry.
    //  Version for non-constant object.
    virtual BeamlineGeometry &getGeometry();

    /// Get geometry.
    //  Version for constant object.
    virtual const BeamlineGeometry &getGeometry() const;

    /// Get arc length.
    //  Return the length of the geometry, measured along the design orbit.
    virtual double getArcLength() const;

    /// Get design length.
    //  Return the length of the geometry, measured along the design polygone.
    virtual double getElementLength() const;

    /// Get origin position.
    //  Return the arc length from the entrance to the origin of the
    //  geometry (non-negative).
    virtual double getOrigin() const;

    /// Get entrance position.
    //  Return the arc length from the origin to the entrance of the
    //  geometry (non-positive).
    virtual double getEntrance() const;

    /// Get exit position.
    //  Return the arc length from the origin to the exit of the
    //  geometry (non-negative).
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
    virtual Euclid3D getExitFrame()     const;

    /// Construct an image.
    //  Return the image of the element, containing the name and type string
    //  of the element, and a copy of the user-defined attributes.
    virtual ElementImage *getImage() const;

    /// Get beamline type
    virtual ElementBase::ElementType getType() const;

    /// Append a T object.
    virtual void append(const T &);

    /// Prepend a T object.
    virtual void prepend(const T &);

    void setOrigin3D(const Vector_t& ori);
    Vector_t getOrigin3D() const;
    void setInitialDirection(const Quaternion& rot);
    Quaternion getInitialDirection() const;

    void setRelativeFlag(bool flag);
    bool getRelativeFlag() const;
protected:

    /// The beamline geometry.
    //  Exists to match the interface for ElementBase.
    BeamlineGeometry itsGeometry;

    Vector_t itsOrigin_m;
    Quaternion itsCoordTrafoTo_m;
    bool relativePositions_m;
};


// Implementation of template class TBeamline
// ------------------------------------------------------------------------

template <class T>
TBeamline<T>::TBeamline():
    Beamline(),
    std::list<T>(),
    itsGeometry(*this),
    itsOrigin_m(0),
    itsCoordTrafoTo_m(1.0, 0.0, 0.0, 0.0),
    relativePositions_m(false)
{}


template <class T>
TBeamline<T>::TBeamline(const std::string &name):
    Beamline(name),
    std::list<T>(),
    itsGeometry(*this),
    itsOrigin_m(0),
    itsCoordTrafoTo_m(1.0, 0.0, 0.0, 0.0),
    relativePositions_m(false)
{}


template <class T>
TBeamline<T>::TBeamline(const TBeamline<T> &rhs):
    Beamline(rhs),
    std::list<T>(rhs),
    itsGeometry(*this),
    itsOrigin_m(rhs.itsOrigin_m),
    itsCoordTrafoTo_m(rhs.itsCoordTrafoTo_m),
    relativePositions_m(rhs.relativePositions_m)
{}


template <class T>
TBeamline<T>::~TBeamline()
{}


template <class T>
void TBeamline<T>::accept(BeamlineVisitor &visitor) const {
    visitor.visitBeamline(*this);
}


template <class T>
void TBeamline<T>::iterate(BeamlineVisitor &visitor, bool r2l) const {
    if(r2l) {
        for(typename std::list<T>::const_reverse_iterator op = this->rbegin();
            op != this->rend(); ++op) {
            op->accept(visitor);
        }
    } else {
        for(typename std::list<T>::const_iterator op = this->begin();
            op != this->end(); ++op) {
            op->accept(visitor);
        }
    }
}


template <class T>
TBeamline<T> *TBeamline<T>::clone() const {
    TBeamline<T> *line = new TBeamline(getName());

    for(typename std::list<T>::const_iterator op = this->begin();
        op != this->end(); ++op) {
        // Make copy of the T object containing a deep copy of its child.
        T newObj(*op);
        newObj.setElement(op->getElement()->clone());
        line->append(newObj);
    }

    line->itsOrigin_m = itsOrigin_m;
    line->itsCoordTrafoTo_m = itsCoordTrafoTo_m;
    line->relativePositions_m = relativePositions_m;

    return line;
}

// 21-6-2000 ada change iterator to const_iterator
template <class T>
TBeamline<T> *TBeamline<T>::copyStructure() {
    if(isSharable()) {
        return this;
    } else {
        TBeamline<T> *line = new TBeamline(getName());

        for(typename std::list<T>::const_iterator iter = this->begin();
            iter != this->end(); ++iter) {
            // The copy constructor ensures proper transmission of data.
            T newObj(*iter);
            newObj.setElement(iter->getElement()->copyStructure());
            line->append(newObj);
        }

        line->itsOrigin_m = itsOrigin_m;
        line->itsCoordTrafoTo_m = itsCoordTrafoTo_m;
        line->relativePositions_m = relativePositions_m;

        return line;
    }
}

// 21-6-2000 ada change iterator to const_iterator
template <class T> inline
void TBeamline<T>::makeSharable() {
    shareFlag = true;

    for(typename std::list<T>::const_iterator iter = this->begin();
        iter != this->end(); ++iter) {
        iter->getElement()->makeSharable();
    }
}


template <class T> inline
BeamlineGeometry &TBeamline<T>::getGeometry() {
    return itsGeometry;
}


template <class T> inline
const BeamlineGeometry &TBeamline<T>::getGeometry() const {
    return itsGeometry;
}


template <class T>
double TBeamline<T>::getArcLength() const {
    double length = 0.0;

    for(typename std::list<T>::const_iterator iter = this->begin();
        iter != this->end(); ++iter) {
        length += iter->getElement()->getArcLength();
    }

    return length;
}


template <class T>
double TBeamline<T>::getElementLength() const {
    double length = 0.0;

    for(typename std::list<T>::const_iterator iter = this->begin();
        iter != this->end(); ++iter) {
        length += iter->getElement()->getElementLength();
    }

    return length;
}


template <class T>
double TBeamline<T>::getOrigin() const {
    return (getArcLength() / 2.0);
}


template <class T>
double TBeamline<T>::getEntrance() const {
    return (- getOrigin());
}


template <class T>
double TBeamline<T>::getExit() const {
    return (getArcLength() / 2.0);
}


template <class T>
Euclid3D TBeamline<T>::getTransform(double fromS, double toS) const {
    Euclid3D transform;

    if(fromS < toS) {
        double s1 = getEntrance();
        typename std::list<T>::const_iterator iter = this->begin();

        while(iter != this->end()  &&  s1 <= toS) {
            const ElementBase &element = *iter->getElement();
            double l = element.getArcLength();
            double s2 = s1 + l;

            if(s2 > fromS) {
                double s0 = (s1 + s2) / 2.0;
                double arc1 = std::max(s1, fromS) - s0;
                double arc2 = std::min(s2, toS)   - s0;
                transform *= element.getTransform(arc1, arc2);
            }

            s1 = s2;
            ++iter;
        }
    } else {
        double s1 = getExit();
        typename std::list<T>::const_reverse_iterator iter = this->rbegin();

        while(iter != this->rend()  &&  s1 >= toS) {
            const ElementBase &element = *iter->getElement();
            double l = element.getArcLength();
            double s2 = s1 - l;

            if(s2 < fromS) {
                double s0 = (s1 + s2) / 2.0;
                double arc1 = std::min(s1, fromS) - s0;
                double arc2 = std::max(s2, toS)   - s0;
                transform *= element.getTransform(arc1, arc2);
            }

            s1 = s2;
            ++iter;
        }
    }

    return transform;
}


template <class T>
Euclid3D TBeamline<T>::getTotalTransform() const {
    Euclid3D transform;

    for(typename std::list<T>::const_iterator iter = this->begin();
        iter != this->end(); ++iter) {
        transform.dotBy(iter->getElement()->getTotalTransform());
    }

    return transform;
}


template <class T>
Euclid3D TBeamline<T>::getTransform(double s) const {
    return getTransform(0.0, s);
}


template <class T>
Euclid3D TBeamline<T>::getEntranceFrame() const {
    return getTransform(0.0, getEntrance());
}


template <class T>
Euclid3D TBeamline<T>::getExitFrame() const {
    return getTransform(0.0, getExit());
}


template <class T> inline
ElementImage *TBeamline<T>::getImage() const {
    // this code needs revision when the class ElementImage is defined
    return NULL;
}


template <class T> inline
ElementBase::ElementType TBeamline<T>::getType() const {
    return BEAMLINE;
}


template <class T> inline
void TBeamline<T>::append(const T &obj) {
    this->push_back(obj);
}


template <class T> inline
void TBeamline<T>::prepend(const T &obj) {
    this->push_front(obj);
}

template <class T> inline
void TBeamline<T>::setOrigin3D(const Vector_t& ori) {
    itsOrigin_m = ori;
}

template <class T> inline
Vector_t TBeamline<T>::getOrigin3D() const {
    return itsOrigin_m;
}

template <class T> inline
void TBeamline<T>::setInitialDirection(const Quaternion& trafoTo) {
    itsCoordTrafoTo_m = trafoTo;
}

template <class T> inline
Quaternion TBeamline<T>::getInitialDirection() const {
    return itsCoordTrafoTo_m;
}

template <class T> inline
void TBeamline<T>::setRelativeFlag(bool flag) {
    relativePositions_m = flag;
}

template <class T> inline
bool TBeamline<T>::getRelativeFlag() const {
    return relativePositions_m;
}
#endif // CLASSIC_TBeamline_HH