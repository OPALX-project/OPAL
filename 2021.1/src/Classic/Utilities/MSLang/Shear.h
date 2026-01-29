#ifndef MSLANG_SHEAR_H
#define MSLANG_SHEAR_H

#include "Utilities/MSLang.h"

namespace mslang {
    struct Shear: public Function {
        Function* func_m;
        double angleX_m;
        double angleY_m;

        virtual ~Shear() {
            delete func_m;
        }

        virtual void print(int indentwidth);
        void applyShear(std::vector<std::shared_ptr<Base> > &bfuncs);
        virtual void apply(std::vector<std::shared_ptr<Base> > &bfuncs);
        static bool parse_detail(iterator &it, const iterator &end, Function* &fun);
    };
}

#endif