#ifndef MSLANG_TRIANGLE_H
#define MSLANG_TRIANGLE_H

#include "Utilities/MSLang.h"

namespace mslang {
    struct Triangle: public Base {
        std::vector<Vector_t> nodes_m;
        Triangle():
            Base(),
            nodes_m(std::vector<Vector_t>(3, Vector_t(0, 0, 1)))
        { }

        Triangle(const Triangle &right):
            Base(right),
            nodes_m(right.nodes_m)
        { }

        virtual ~Triangle()
        { }

        virtual void print(int indentwidth);
        virtual void apply(std::vector<std::shared_ptr<Base> > &bfuncs);
        virtual std::shared_ptr<Base> clone() const;
        virtual void writeGnuplot(std::ofstream &out) const;
        virtual void computeBoundingBox();
        double crossProduct(const Vector_t &pt, unsigned int nodeNum) const;
        virtual bool isInside(const Vector_t &R) const;
        void orientNodesCCW();
    };
}

#endif