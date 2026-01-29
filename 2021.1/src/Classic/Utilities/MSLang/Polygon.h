#ifndef MSLANG_POLYGON_H
#define MSLANG_POLYGON_H

#include "Utilities/MSLang.h"
#include "Utilities/MSLang/Triangle.h"

namespace mslang {
    struct Polygon: public Function {
        std::vector<Triangle> triangles_m;

        void triangulize(std::vector<Vector_t> &nodes);
        static bool parse_detail(iterator &it, const iterator &end, Function* &fun);
        virtual void print(int ident);
        virtual void apply(std::vector<std::shared_ptr<Base> > &bfuncs);
    };
}

#endif