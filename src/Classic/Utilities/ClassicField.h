#ifndef CLASSIC_FIELD_H
#define CLASSIC_FIELD_H

#include <list>
#include <memory>
#include "AbsBeamline/Component.h"
#include "Structure/BoundingBox.h"

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
        return (fle1.start_m < fle2.start_m
                || (fle1.start_m == fle2.start_m && fle1.element_m->getName() < fle2.element_m->getName()));
    }

    static bool ZeroLength(const ClassicField &fle) {
        return (fle.getLength() < 1.e-6);
    }

    BoundingBox getBoundingBoxInLabCoords() const;

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
BoundingBox ClassicField::getBoundingBoxInLabCoords() const {
    return element_m->getBoundingBoxInLabCoords();
}
#endif // CLASSIC_FIELD_H