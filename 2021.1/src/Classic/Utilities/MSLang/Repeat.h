#ifndef MSLANG_REPEAT_H
#define MSLANG_REPEAT_H

#include "Utilities/MSLang.h"

namespace mslang {
    struct Repeat: public Function {
        Function* func_m;
        unsigned int N_m;
        double shiftx_m;
        double shifty_m;
        double rot_m;

        virtual ~Repeat() {
            delete func_m;
        }

        virtual void print(int indentwidth);
        virtual void apply(std::vector<std::shared_ptr<Base> > &bfuncs);
        static bool parse_detail(iterator &it, const iterator &end, Function* &fun);
    };
}

#endif