#ifndef MSLANG_H
#define MSLANG_H

#include "Utilities/MSLang/BoundingBox.h"
#include "Utilities/MSLang/AffineTransformation.h"
#include "Algorithms/Vektor.h"
#include "AppTypes/Tenzor.h"

#include <string>
#include <iostream>
#include <fstream>
#include <list>
#include <memory>

namespace mslang {
    typedef std::string::iterator iterator;

    inline
    double euclidean_norm2D(Vector_t v) {
        v[2] = 0.0;
        return euclidean_norm(v);
    }

    struct Base;

    struct Function {
        virtual ~Function() {};

        virtual void print(int indent) = 0;
        virtual void apply(std::vector<std::shared_ptr<Base>> &bfuncs) = 0;

        static bool parse(iterator &it, const iterator &end, Function* &fun);

        static const std::string UDouble;
        static const std::string Double;
        static const std::string UInt;
        static const std::string FCall;
    };

    struct Base: public Function {
        AffineTransformation trafo_m;
        BoundingBox bb_m;
        std::vector<std::shared_ptr<Base> > divisor_m;

        Base():
            trafo_m()
        { }

        Base(const Base &right):
            trafo_m(right.trafo_m),
            bb_m(right.bb_m)
        { }

        virtual ~Base() {
            // for (auto item: divisor_m) {
            //     item.reset();
            // }
            divisor_m.clear();
        }

        virtual std::shared_ptr<Base> clone() const = 0;
        virtual void writeGnuplot(std::ofstream &out) const = 0;
        virtual void computeBoundingBox() = 0;
        virtual bool isInside(const Vector_t &R) const = 0;
        virtual void divideBy(std::vector<std::shared_ptr<Base> > &divisors) {
            for (auto item: divisors) {
                if (bb_m.doesIntersect(item->bb_m)) {
                    divisor_m.emplace_back(item->clone());
                }
            }
        }
    };

    bool parse(std::string str, Function* &fun);
}

#endif
