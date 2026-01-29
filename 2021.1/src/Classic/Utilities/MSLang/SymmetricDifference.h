#ifndef MSLANG_SYMMETRICDIFFERENCE_H
#define MSLANG_SYMMETRICDIFFERENCE_H

#include "Utilities/MSLang.h"

namespace mslang {
    struct SymmetricDifference: public Function {
        Function *firstOperand_m;
        Function *secondOperand_m;

        SymmetricDifference()
        { }

        SymmetricDifference(const SymmetricDifference &right):
            firstOperand_m(right.firstOperand_m),
            secondOperand_m(right.secondOperand_m)
        { }

        virtual ~SymmetricDifference() {
            delete firstOperand_m;
            delete secondOperand_m;
        }

        virtual void print(int indentwidth);
        virtual void apply(std::vector<std::shared_ptr<Base> > &bfuncs);
        static bool parse_detail(iterator &it, const iterator &end, Function* &fun);
    };
}

#endif