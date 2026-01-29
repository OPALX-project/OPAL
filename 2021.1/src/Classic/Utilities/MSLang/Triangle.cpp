#include "Utilities/MSLang/Triangle.h"
#include "Physics/Physics.h"

namespace mslang {
    void Triangle::print(int indentwidth) {
        std::string indent(indentwidth, ' ');
        std::string indentbase(4, ' ');
        Vector_t origin = trafo_m.getOrigin();
        double angle = trafo_m.getAngle() * Physics::rad2deg;

        std::cout << indent << "triangle, \n";

        for (unsigned int i = 0; i < 3u; ++ i) {
            std::cout << indent << indentbase << "node " << i << ": " << nodes_m[i] << "\n";
        }
        std::cout << indent << indentbase << "origin: " << origin[0] << ", " << origin[1] << ",\n"
                  << indent << indentbase << "angle: " << angle << "\n"
                  << indent << indentbase << trafo_m(0, 0) << "\t" << trafo_m(0, 1) << "\t" << trafo_m(0, 2) << "\n"
                  << indent << indentbase << trafo_m(1, 0) << "\t" << trafo_m(1, 1) << "\t" << trafo_m(1, 2) << "\n"
                  << indent << indentbase << trafo_m(2, 0) << "\t" << trafo_m(2, 1) << "\t" << trafo_m(2, 2) << std::endl;
    }
    void Triangle::apply(std::vector<std::shared_ptr<Base> > &bfuncs) {
        bfuncs.emplace_back(this->clone());
    }

    std::shared_ptr<Base>  Triangle::clone() const {
        std::shared_ptr<Triangle> tri(new Triangle(*this));
        tri->trafo_m = trafo_m;
        tri->bb_m = bb_m;

        for (auto item: divisor_m) {
            tri->divisor_m.emplace_back(item->clone());
        }

        return std::static_pointer_cast<Base>(tri);
    }

    void Triangle::writeGnuplot(std::ofstream &out) const {
        unsigned int width = out.precision() + 8;

        for (unsigned int i = 0; i < 4u; ++ i) {
            Vector_t corner = trafo_m.transformFrom(nodes_m[i % 3u]);

            out << std::setw(width) << corner[0]
                << std::setw(width) << corner[1]
                << std::endl;
        }
        out << std::endl;

        for (auto item: divisor_m) {
            item->writeGnuplot(out);
        }

        // bb_m.writeGnuplot(out);
    }

    void Triangle::computeBoundingBox() {
        std::vector<Vector_t> corners;
        for (unsigned int i = 0; i < 3u; ++ i) {
            corners.push_back(trafo_m.transformFrom(nodes_m[i]));
        }

        Vector_t llc = corners[0], urc = corners[0];
        for (unsigned int i = 1u; i < 3u; ++ i) {
            if (corners[i][0] < llc[0]) llc[0] = corners[i][0];
            else if (corners[i][0] > urc[0]) urc[0] = corners[i][0];

            if (corners[i][1] < llc[1]) llc[1] = corners[i][1];
            else if (corners[i][1] > urc[1]) urc[1] = corners[i][1];
        }

        bb_m = BoundingBox(llc, urc);

        for (auto item: divisor_m) {
            item->computeBoundingBox();
        }
    }

    double Triangle::crossProduct(const Vector_t &pt, unsigned int nodeNum) const {
        nodeNum = nodeNum % 3u;
        unsigned int nextNode = (nodeNum + 1) % 3u;
        Vector_t nodeToPt = pt - nodes_m[nodeNum];
        Vector_t nodeToNext = nodes_m[nextNode] - nodes_m[nodeNum];

        return nodeToPt[0] * nodeToNext[1] - nodeToPt[1] * nodeToNext[0];

    }

    bool Triangle::isInside(const Vector_t &R) const {
        Vector_t X = trafo_m.transformTo(R);

        bool test0 = (crossProduct(X, 0) <= 0.0);
        bool test1 = (crossProduct(X, 1) <= 0.0);
        bool test2 = (crossProduct(X, 2) <= 0.0);

        if (!(test0 && test1 && test2)) return false;

        for (auto item: divisor_m)
            if (item->isInside(R)) return false;

        return true;
    }

    void Triangle::orientNodesCCW() {
        if (crossProduct(nodes_m[0], 1) > 0.0) {
            std::swap(nodes_m[1], nodes_m[2]);
        }
    }
}
