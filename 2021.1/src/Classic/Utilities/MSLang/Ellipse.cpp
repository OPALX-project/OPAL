#include "Utilities/MSLang/Ellipse.h"
#include "Utilities/MSLang/ArgumentExtractor.h"
#include "Utilities/MSLang/matheval.h"
#include "Physics/Physics.h"

#include <boost/regex.hpp>

namespace mslang {
    void Ellipse::print(int indentwidth) {
        std::string indent(indentwidth, ' ');
        std::string indent2(indentwidth + 8, ' ');
        Vector_t origin = trafo_m.getOrigin();
        double angle = trafo_m.getAngle() * Physics::rad2deg;
        std::cout << indent << "ellipse, \n"
                  << indent2 << "w: " << width_m << ", \n"
                  << indent2 << "h: " << height_m << ", \n"
                  << indent2 << "origin: " << origin[0] << ", " << origin[1] << ",\n"
                  << indent2 << "angle: " << angle << "\n"
                  << indent2 << std::setw(14) << trafo_m(0, 0) << std::setw(14) << trafo_m(0, 1) << std::setw(14) << trafo_m(0, 2) << "\n"
                  << indent2 << std::setw(14) << trafo_m(1, 0) << std::setw(14) << trafo_m(1, 1) << std::setw(14) << trafo_m(1, 2) << "\n"
                  << indent2 << std::setw(14) << trafo_m(2, 0) << std::setw(14) << trafo_m(2, 1) << std::setw(14) << trafo_m(2, 2)
                  << std::endl;
    }

    void Ellipse::writeGnuplot(std::ofstream &out) const {
        const unsigned int N = 101;
        const double dp = Physics::two_pi / (N - 1);
        const unsigned int colwidth = out.precision() + 8;

        double phi = 0;
        for (unsigned int i = 0; i < N; ++ i, phi += dp) {
            Vector_t pt(0.0);
            pt[0] = std::copysign(sqrt(std::pow(height_m * width_m * 0.25, 2) /
                                       (std::pow(height_m * 0.5, 2) +
                                        std::pow(width_m * 0.5 * tan(phi), 2))),
                                  cos(phi));
            pt[1] = pt[0] * tan(phi);
            pt = trafo_m.transformFrom(pt);

            out << std::setw(colwidth) << pt[0]
                << std::setw(colwidth) << pt[1]
                << std::endl;
        }
        out << std::endl;

        for (auto item: divisor_m) {
            item->writeGnuplot(out);
        }

        // bb_m.writeGnuplot(out);
    }

    void Ellipse::apply(std::vector<std::shared_ptr<Base> > &bfuncs) {
        bfuncs.emplace_back(this->clone());
    }

    std::shared_ptr<Base> Ellipse::clone() const{
        std::shared_ptr<Ellipse> elps(new Ellipse);
        elps->width_m = width_m;
        elps->height_m = height_m;
        elps->trafo_m = trafo_m;
        elps->bb_m = bb_m;

        for (auto item: divisor_m) {
            elps->divisor_m.emplace_back(item->clone());
        }

        return std::static_pointer_cast<Base>(elps);
    }

    void Ellipse::computeBoundingBox() {
        Vector_t llc(0.0), urc(0.0);
        const Vector_t e_x(1.0, 0.0, 0.0), e_y(0.0, 1.0, 0.0);
        const Vector_t center = trafo_m.transformFrom(Vector_t(0.0));
        const Vector_t e_xp = trafo_m.transformFrom(e_x) - center;
        const Vector_t e_yp = trafo_m.transformFrom(e_y) - center;
        const double &M11 = e_xp[0];
        const double &M12 = e_yp[0];
        const double &M21 = e_xp[1];
        const double &M22 = e_yp[1];

        double t = atan2(height_m * M12, width_m * M11);
        double halfwidth = 0.5 * (M11 * width_m * cos(t) +
                                  M12 * height_m * sin(t));
        llc[0] = center[0] - std::abs(halfwidth);
        urc[0] = center[0] + std::abs(halfwidth);

        t = atan2(height_m * M22, width_m * M21);

        double halfheight = 0.5 * (M21 * width_m * cos(t) +
                                   M22 * height_m * sin(t));

        llc[1] = center[1] - std::abs(halfheight);
        urc[1] = center[1] + std::abs(halfheight);

        bb_m = BoundingBox(llc, urc);

        for (auto item: divisor_m) {
            item->computeBoundingBox();
        }
    }

    bool Ellipse::isInside(const Vector_t &R) const {
        if (!bb_m.isInside(R)) return false;

        Vector_t X = trafo_m.transformTo(R);
        if (4 * (std::pow(X[0] / width_m, 2) + std::pow(X[1] / height_m, 2)) <= 1) {

            for (auto item: divisor_m) {
                if (item->isInside(R)) return false;
            }

            return true;
        }

        return false;
    }

    bool Ellipse::parse_detail(iterator &it, const iterator &end, Function* fun) {

        ArgumentExtractor arguments(std::string(it, end));
        double width, height;
        try {
            width = parseMathExpression(arguments.get(0));
            height = parseMathExpression(arguments.get(1));
        } catch (std::runtime_error &e) {
            std::cout << e.what() << std::endl;
            return false;
        }

        Ellipse *elps = static_cast<Ellipse*>(fun);
        elps->width_m = width;
        elps->height_m = height;

        if (elps->width_m < 0.0) {
            std::cout << "Ellipse: a negative width provided '"
                      << arguments.get(0) << " = " << elps->width_m << "'"
                      << std::endl;
            return false;
        }
        if (elps->height_m < 0.0) {
            std::cout << "Ellipse: a negative height provided '"
                      << arguments.get(1) << " = " << elps->height_m << "'"
                      << std::endl;
            return false;
        }

        it += (arguments.getLengthConsumed() + 1);

        return true;
    }
}
