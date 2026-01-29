#ifndef MSLANG_RECTANGLE_H
#define MSLANG_RECTANGLE_H

#include "Utilities/MSLang.h"

namespace mslang {
    struct Rectangle: public Base {
        double width_m;
        double height_m;

        Rectangle():
            Base(),
            width_m(0.0),
            height_m(0.0)
        { }

        Rectangle(const Rectangle &right):
            Base(right),
            width_m(right.width_m),
            height_m(right.height_m)
        { }

        virtual ~Rectangle() { }

        virtual void print(int indentwidth);
        virtual void computeBoundingBox();
        virtual bool isInside(const Vector_t &R) const;
        virtual void writeGnuplot(std::ofstream &out) const;
        virtual void apply(std::vector<std::shared_ptr<Base> > &bfuncs);
        virtual std::shared_ptr<Base> clone() const;
        static bool parse_detail(iterator &it, const iterator &end, Function* fun);
    };
}

#endif