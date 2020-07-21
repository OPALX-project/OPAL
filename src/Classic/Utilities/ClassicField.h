#ifndef CLASSIC_FIELD_H
#define CLASSIC_FIELD_H

#include <vector>
#include <list>
#include <memory>
#include "AbsBeamline/Component.h"
#include "Algorithms/Quaternion.h"
#include "Algorithms/CoordinateSystemTrafo.h"

class ClassicField {
public:
    ClassicField(std::shared_ptr<Component>, const double &, const double &);
    ~ClassicField();
    std::shared_ptr<Component> getElement();
    std::shared_ptr<const Component> getElement() const;
    double getLength() const;
    const double &getStart() const;
    const double &getEnd() const;
    void setStart(const double & z);
    void setEnd(const double & z);
    const bool &isOn() const;
    void setOn(const double &kinematicEnergy);
    void setOff();

    static bool SortAsc(const ClassicField &fle1, const ClassicField &fle2) {
        return (fle1.start_m < fle2.start_m);
    }

    static bool ZeroLength(const ClassicField &fle) {
        return (fle.getLength() < 1.e-6);
    }

    ElementBase::BoundingBox getBoundingBoxInLabCoords() const;

    CoordinateSystemTrafo getCoordTransformationTo() const ;
    void setCoordTransformationTo(const CoordinateSystemTrafo &trafo);
    bool isPositioned() const;
    void fixPosition();
    unsigned int order_m;
private:
    std::shared_ptr<Component> element_m;
    double start_m;
    double end_m;
    bool is_on_m;
};

typedef std::list<ClassicField> FieldList;

inline std::shared_ptr<Component> ClassicField::getElement() {
    return element_m;
}

inline std::shared_ptr<const Component> ClassicField::getElement() const {
    return element_m;
}

inline double ClassicField::getLength() const {
    return end_m - start_m;
}

inline const double &ClassicField::getStart() const {
    return start_m;
}

inline const double &ClassicField::getEnd() const {
    return end_m;
}

inline const bool &ClassicField::isOn() const {
    return is_on_m;
}

inline void ClassicField::setStart(const double & z) {
    start_m = z;
}

inline void ClassicField::setEnd(const double & z) {
    end_m = z;
}

inline
CoordinateSystemTrafo ClassicField::getCoordTransformationTo() const {
    return element_m->getCSTrafoGlobal2Local();
}

inline
void ClassicField::setCoordTransformationTo(const CoordinateSystemTrafo &trafo) {
    element_m->setCSTrafoGlobal2Local(trafo);
}

inline
bool ClassicField::isPositioned() const {
    return element_m->isPositioned();
}

inline
void ClassicField::fixPosition() {
    element_m->fixPosition();
}

inline
ElementBase::BoundingBox ClassicField::getBoundingBoxInLabCoords() const {
    return element_m->getBoundingBoxInLabCoords();
}
#endif // CLASSIC_FIELD_H
