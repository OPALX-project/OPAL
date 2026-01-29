#ifndef MSLANG_DIFFERENCE_H
#define MSLANG_DIFFERENCE_H

#include "Utilities/MSLang.h"

namespace mslang {
    struct Difference: public Function {
        Function* dividend_m;
        Function* divisor_m;

        Difference()
        { }

        Difference(const Difference &right):
            dividend_m(right.dividend_m),
            divisor_m(right.divisor_m)
        { }

        virtual ~Difference() {
            delete dividend_m;
            delete divisor_m;
        }

        virtual void print(int indentwidth);
        virtual void apply(std::vector<std::shared_ptr<Base> > &bfuncs);
        static bool parse_detail(iterator &it, const iterator &end, Function* &fun);
    };
}

#endif