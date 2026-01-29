#ifndef MSLANG_UNION_H
#define MSLANG_UNION_H

#include "Utilities/MSLang.h"

namespace mslang {
    struct Union: public Function {
        std::vector<Function*> funcs_m;

        virtual ~Union () {
            for (Function* func: funcs_m) {
                delete func;
            }
        }

        virtual void print(int indentwidth);
        virtual void apply(std::vector<std::shared_ptr<Base> > &bfuncs);
        static bool parse_detail(iterator &it, const iterator &end, Function* &fun);
    };
}

#endif