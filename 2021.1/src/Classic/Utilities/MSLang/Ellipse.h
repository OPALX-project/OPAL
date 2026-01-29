#ifndef MSLANG_ELLIPSE_H
#define MSLANG_ELLIPSE_H

#include "Utilities/MSLang.h"

namespace mslang {
    struct Ellipse: public Base {
        double width_m;
        double height_m;

        Ellipse():
            Base(),
            width_m(0.0),
            height_m(0.0)
        { }

        Ellipse(const Ellipse &right):
            Base(right),
            width_m(right.width_m),
            height_m(right.height_m)
        { }

        virtual ~Ellipse() { }

        virtual void print(int indentwidth);
        virtual void writeGnuplot(std::ofstream &out) const;
        virtual void apply(std::vector<std::shared_ptr<Base> > &bfuncs);
        virtual std::shared_ptr<Base> clone() const;
        virtual void computeBoundingBox();
        virtual bool isInside(const Vector_t &R) const;
        static bool parse_detail(iterator &it, const iterator &end, Function* fun);
    };
}

#endif