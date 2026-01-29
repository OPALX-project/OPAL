#include "Utilities/MSLang/Shear.h"
#include "Utilities/MSLang/ArgumentExtractor.h"
#include "Utilities/MSLang/matheval.h"

#include <boost/regex.hpp>

namespace mslang {
    void Shear::print(int indentwidth) {
        std::string indent(indentwidth, ' ');
        std::string indent2(indentwidth + 8, ' ');
        std::cout << indent << "shear, " << std::endl;
        func_m->print(indentwidth + 8);
        if (std::abs(angleX_m) > 0.0) {
            std::cout << ",\n"
                      << indent2 << "angle X: " << angleX_m;
        } else {
            std::cout << ",\n"
                      << indent2 << "angle Y: " << angleY_m;
        }
    }

    void Shear::applyShear(std::vector<std::shared_ptr<Base> > &bfuncs) {
        AffineTransformation shear(Vector_t(1.0, tan(angleX_m), 0.0),
                                   Vector_t(-tan(angleY_m), 1.0, 0.0));

        const unsigned int size = bfuncs.size();

        for (unsigned int j = 0; j < size; ++ j) {
            std::shared_ptr<Base> &obj = bfuncs[j];
            obj->trafo_m = obj->trafo_m.mult(shear);

            if (obj->divisor_m.size() > 0)
                applyShear(obj->divisor_m);
        }
    }

    void Shear::apply(std::vector<std::shared_ptr<Base> > &bfuncs) {
        func_m->apply(bfuncs);
        applyShear(bfuncs);
    }

    bool Shear::parse_detail(iterator &it, const iterator &end, Function* &fun) {
        Shear *shr = static_cast<Shear*>(fun);
        if (!parse(it, end, shr->func_m)) return false;

        ArgumentExtractor arguments(std::string(++ it, end));
        try {
            shr->angleX_m = parseMathExpression(arguments.get(0));
            shr->angleY_m = parseMathExpression(arguments.get(1));
        } catch (std::runtime_error &e) {
            std::cout << e.what() << std::endl;
            return false;
        }

        it += (arguments.getLengthConsumed() + 1);

        return true;
    }
}