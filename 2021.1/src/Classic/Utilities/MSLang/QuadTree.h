#ifndef MSLANG_QUADTREE_H
#define MSLANG_QUADTREE_H

#include "Utilities/MSLang.h"

#include <memory>

namespace mslang {
    struct QuadTree {
        int level_m;
        std::list<std::shared_ptr<Base> > objects_m;
        BoundingBox bb_m;
        std::vector<std::shared_ptr<QuadTree> > nodes_m;

        QuadTree():
            level_m(0),
            bb_m(),
            nodes_m()
        { }

        QuadTree(int l, const BoundingBox &b):
            level_m(l),
            bb_m(b),
            nodes_m()
        { }

        QuadTree(const QuadTree &right);

        ~QuadTree();

        void reset();

        void operator=(const QuadTree &right);

        void transferIfInside(std::list<std::shared_ptr<Base> > &objs);
        void buildUp();

        void writeGnuplot(std::ostream &out) const;

        bool isInside(const Vector_t &R) const;

        void getDepth(unsigned int &d) const;
    };
}

#endif