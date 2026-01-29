#include "Utilities/MSLang/QuadTree.h"
#include <utility>

namespace mslang {
    QuadTree::QuadTree(const QuadTree &right):
        level_m(right.level_m),
        objects_m(right.objects_m.begin(),
                  right.objects_m.end()),
        bb_m(right.bb_m),
        nodes_m()
    {
        if (right.nodes_m.size() != 0) {
            for (unsigned int i = 0; i < 4u; ++ i) {
                nodes_m.emplace_back(new QuadTree(*right.nodes_m[i]));
            }
        }
    }

    QuadTree::~QuadTree() {
        // for (std::shared_ptr<Base> &obj: objects_m)
        //     obj.reset();

        objects_m.clear();
        nodes_m.clear();
    }

    void QuadTree::reset() {
        objects_m.clear();
        nodes_m.clear();
    }

    void QuadTree::operator=(const QuadTree &right) {
        level_m = right.level_m;
        objects_m.insert(objects_m.end(),
                         right.objects_m.begin(),
                         right.objects_m.end());
        bb_m = right.bb_m;

        if (nodes_m.size() != 0) nodes_m.clear();

        if (right.nodes_m.size() != 0) {
            for (unsigned int i = 0; i < 4u; ++ i) {
                nodes_m.emplace_back(new QuadTree(*right.nodes_m[i]));
            }
        }
    }

    void QuadTree::transferIfInside(std::list<std::shared_ptr<Base> > &objs) {
        for (std::shared_ptr<Base> &obj: objs) {
            if (bb_m.isInside(obj->bb_m)) {
                objects_m.emplace_back(std::move(obj));
            }
        }

        objs.remove_if([](const std::shared_ptr<Base> obj) { return !obj; });
    }

    void QuadTree::buildUp() {
        double X[] = {bb_m.center_m[0] - 0.5 * bb_m.width_m,
                      bb_m.center_m[0],
                      bb_m.center_m[0] + 0.5 * bb_m.width_m};
        double Y[] = {bb_m.center_m[1] - 0.5 * bb_m.height_m,
                      bb_m.center_m[1],
                      bb_m.center_m[1] + 0.5 * bb_m.height_m};

        bool allEmpty = true;

        nodes_m.reserve(4);
        for (unsigned int i = 0; i < 4u; ++ i) {
            nodes_m.emplace_back(new QuadTree(level_m + 1,
                                              BoundingBox(Vector_t(X[i / 2], Y[i % 2], 0.0),
                                                          Vector_t(X[i / 2 + 1], Y[i % 2 + 1], 0.0))));

            nodes_m.back()->transferIfInside(objects_m);

            if (nodes_m.back()->objects_m.size() != 0) {
                allEmpty = false;
            }
        }

        if (allEmpty) {
            nodes_m.clear();
            return;
        }

        for (unsigned int i = 0; i < 4u; ++ i) {
            nodes_m[i]->buildUp();
        }
    }

    void QuadTree::writeGnuplot(std::ostream &out) const {
        out << "# level: " << level_m << ", size: " << objects_m.size() << std::endl;
        bb_m.writeGnuplot(out);
        out << "# num holes: " << objects_m.size() << std::endl;
        out << std::endl;

        if (nodes_m.size() != 0) {
            for (unsigned int i = 0; i < 4u; ++ i) {
                nodes_m[i]->writeGnuplot(out);
            }
        }
    }

    bool QuadTree::isInside(const Vector_t &R) const {
        if (nodes_m.size() != 0) {
            Vector_t X = R - bb_m.center_m;
            unsigned int idx = (X[1] < 0.0 ? 0: 1);
            idx += (X[0] < 0.0 ? 0: 2);

            if (nodes_m[idx]->isInside(R)) {
                return true;
            }
        }

        for (const std::shared_ptr<Base> & obj: objects_m) {
            if (obj->isInside(R)) {
                return true;
            }
        }

        return false;
    }

    void QuadTree::getDepth(unsigned int &d) const {
        if (nodes_m.size() > 0) {
            for (const auto & node: nodes_m) {
                node->getDepth(d);
            }
        } else {
            if ((unsigned int)level_m > d) d = level_m;
        }
    }
}