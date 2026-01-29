#include "Utilities/MSLang/Translation.h"
#include "Utilities/MSLang/ArgumentExtractor.h"
#include "Utilities/MSLang/matheval.h"
#include "Physics/Physics.h"

#include <boost/regex.hpp>
#include <string>

namespace mslang {
    void Translation::print(int indentwidth) {
        std::string indent(indentwidth, ' ');
        std::string indent2(indentwidth + 8, ' ');
        std::cout << indent << "translate, " << std::endl;
        func_m->print(indentwidth + 8);
        std::cout << ",\n"
                  << indent2 << "dx: " << shiftx_m << ", \n"
                  << indent2 << "dy: " << shifty_m;
    }

    void Translation::applyTranslation(std::vector<std::shared_ptr<Base> > &bfuncs) {
        AffineTransformation shift(Vector_t(1.0, 0.0, -shiftx_m),
                                   Vector_t(0.0, 1.0, -shifty_m));

        const unsigned int size = bfuncs.size();
        for (unsigned int j = 0; j < size; ++ j) {
            std::shared_ptr<Base> &obj = bfuncs[j];
            obj->trafo_m = obj->trafo_m.mult(shift);

            if (obj->divisor_m.size() > 0)
                applyTranslation(obj->divisor_m);
        }
    }

    void Translation::apply(std::vector<std::shared_ptr<Base> > &bfuncs) {
        func_m->apply(bfuncs);
        applyTranslation(bfuncs);
    }

    bool Translation::parse_detail(iterator &it, const iterator &end, Function* &fun) {
        Translation *trans = static_cast<Translation*>(fun);
        if (!parse(it, end, trans->func_m)) return false;

        ArgumentExtractor arguments(std::string(++ it, end));
        try {
            trans->shiftx_m = parseMathExpression(arguments.get(0));
            trans->shifty_m = parseMathExpression(arguments.get(1));
        } catch (std::runtime_error &e) {
            std::cout << e.what() << std::endl;
            return false;
        }

        it += (arguments.getLengthConsumed() + 1);

        return true;
    }
}