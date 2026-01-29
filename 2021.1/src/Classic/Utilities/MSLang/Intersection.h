#ifndef MSLANG_INTERSECTION_H
#define MSLANG_INTERSECTION_H

#include "Utilities/MSLang.h"

namespace mslang {
    struct Intersection: public Function {
        Function *firstOperand_m;
        Function *secondOperand_m;

        Intersection()
        { }

        Intersection(const Intersection &right):
            firstOperand_m(right.firstOperand_m),
            secondOperand_m(right.secondOperand_m)
        { }

        virtual ~Intersection() {
            delete firstOperand_m;
            delete secondOperand_m;
        }

        virtual void print(int indentwidth);
        virtual void apply(std::vector<std::shared_ptr<Base> > &bfuncs);
        static bool parse_detail(iterator &it, const iterator &end, Function* &fun);
    };
}

#endif