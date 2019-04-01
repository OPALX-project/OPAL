#include "Utilities/MSLang/Rectangle.h"
#include "Utilities/MSLang/ArgumentExtractor.h"
#include "Utilities/MSLang/matheval.h"
#include "Physics/Physics.h"

#include <boost/regex.hpp>

namespace mslang {
    void Rectangle::print(int indentwidth) {
        std::string indent(indentwidth, ' ');
        std::string indent2(indentwidth + 8, ' ');
        Vector_t origin = trafo_m.getOrigin();
        double angle = trafo_m.getAngle() * Physics::rad2deg;
        std::cout << indent << "rectangle, \n"
                  << indent2 << "w: " << width_m << ", \n"
                  << indent2 << "h: " << height_m << ", \n"
                  << indent2 << "origin: " << origin[0] << ", " << origin[1] << ",\n"
                  << indent2 << "angle: " << angle << "\n"
                  << indent2 << trafo_m(0, 0) << "\t" << trafo_m(0, 1) << "\t" << trafo_m(0, 2) << "\n"
                  << indent2 << trafo_m(1, 0) << "\t" << trafo_m(1, 1) << "\t" << trafo_m(1, 2) << "\n"
                  << indent2 << trafo_m(2, 0) << "\t" << trafo_m(2, 1) << "\t" << trafo_m(2, 2) << std::endl;
    }

    void Rectangle::computeBoundingBox() {
        std::vector<Vector_t> corners({Vector_t(0.5 * width_m, 0.5 * height_m, 0),
                    Vector_t(-0.5 * width_m, 0.5 * height_m, 0),
                    Vector_t(-0.5 * width_m, -0.5 * height_m, 0),
                    Vector_t(0.5 * width_m, -0.5 * height_m, 0)});

        for (Vector_t &v: corners) {
            v = trafo_m.transformFrom(v);
        }

        Vector_t llc = corners[0], urc = corners[0];
        for (unsigned int i = 1; i < 4; ++ i) {
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

    bool Rectangle::isInside(const Vector_t &R) const {
        if (!bb_m.isInside(R)) return false;

        Vector_t X = trafo_m.transformTo(R);
        if (2 * std::abs(X[0]) <= width_m &&
            2 * std::abs(X[1]) <= height_m) {
            for (auto item: divisor_m) {
                if (item->isInside(R))
                    return false;
            }
            return true;
        }

        return false;
    }

    void Rectangle::writeGnuplot(std::ofstream &out) const {
        std::vector<Vector_t> pts({Vector_t(0.5 * width_m, 0.5 * height_m, 0),
                    Vector_t(-0.5 * width_m, 0.5 * height_m, 0),
                    Vector_t(-0.5 * width_m, -0.5 * height_m, 0),
                    Vector_t(0.5 * width_m, -0.5 * height_m, 0)});
        unsigned int width = out.precision() + 8;
        for (unsigned int i = 0; i < 5; ++ i) {
            Vector_t pt = pts[i % 4];
            pt = trafo_m.transformFrom(pt);

            out << std::setw(width) << pt[0]
                << std::setw(width) << pt[1]
                << std::endl;
        }
        out << std::endl;

        for (auto item: divisor_m) {
            item->writeGnuplot(out);
        }

        // bb_m.writeGnuplot(out);
    }

    void Rectangle::apply(std::vector<std::shared_ptr<Base> > &bfuncs) {
        bfuncs.emplace_back(this->clone());
    }

    std::shared_ptr<Base> Rectangle::clone() const {
        std::shared_ptr<Rectangle> rect(new Rectangle);
        rect->width_m = width_m;
        rect->height_m = height_m;
        rect->trafo_m = trafo_m;
        rect->bb_m = bb_m;

        for (auto item: divisor_m) {
            rect->divisor_m.emplace_back(item->clone());
        }

        return std::static_pointer_cast<Base>(rect);
    }

    bool Rectangle::parse_detail(iterator &it, const iterator &end, Function* fun) {
        std::string str(it, end);
        boost::regex argumentList(UDouble + "," + UDouble + "(\\).*)");
        boost::smatch what;

        Rectangle *rect = static_cast<Rectangle*>(fun);
        ArgumentExtractor arguments(str);
        try {
            rect->width_m = parseMathExpression(arguments.get(0));
            rect->height_m = parseMathExpression(arguments.get(1));
        } catch (std::runtime_error &e) {
            std::cout << e.what() << std::endl;
            return false;
        }

        if (rect->width_m < 0.0) {
            std::cout << "Rectangle: a negative width provided '"
                      << arguments.get(0) << " = " << rect->width_m << "'"
                      << std::endl;
            return false;
        }
        if (rect->height_m < 0.0) {
            std::cout << "Rectangle: a negative height provided '"
                      << arguments.get(1) << " = " << rect->height_m << "'"
                      << std::endl;
            return false;
        }

        it += (arguments.getLengthConsumed() + 1);

        return true;
    }
}
