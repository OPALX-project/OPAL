#include "Utilities/MSLang/Repeat.h"
#include "Utilities/MSLang/ArgumentExtractor.h"
#include "Utilities/MSLang/matheval.h"

namespace mslang {
    void Repeat::print(int indentwidth) {
        std::string indent(indentwidth, ' ');
        std::string indent2(indentwidth + 8, ' ');
        std::cout << indent << "repeat, " << std::endl;
        func_m->print(indentwidth + 8);
        std::cout << ",\n"
                  << indent2 << "N: " << N_m << ", \n"
                  << indent2 << "dx: " << shiftx_m << ", \n"
                  << indent2 << "dy: " << shifty_m;
    }

    void Repeat::apply(std::vector<std::shared_ptr<Base> > &bfuncs) {
        AffineTransformation trafo(Vector_t(cos(rot_m), sin(rot_m), -shiftx_m),
                                   Vector_t(-sin(rot_m), cos(rot_m), -shifty_m));

        func_m->apply(bfuncs);
        const unsigned int size = bfuncs.size();

        AffineTransformation current_trafo = trafo;
        for (unsigned int i = 0; i < N_m; ++ i) {
            for (unsigned int j = 0; j < size; ++ j) {
                std::shared_ptr<Base> obj(bfuncs[j]->clone());
                obj->trafo_m = obj->trafo_m.mult(current_trafo);
                bfuncs.emplace_back(std::move(obj));
            }

            current_trafo = current_trafo.mult(trafo);
        }
    }

    bool Repeat::parse_detail(iterator &it, const iterator &end, Function* &fun) {
        Repeat *rep = static_cast<Repeat*>(fun);
        if (!parse(it, end, rep->func_m)) return false;

        int numRepetitions = 0;
        ArgumentExtractor arguments(std::string(++ it, end));
        try {
            if (arguments.getNumArguments() == 3) {
                numRepetitions = parseMathExpression(arguments.get(0));
                rep->shiftx_m = parseMathExpression(arguments.get(1));
                rep->shifty_m = parseMathExpression(arguments.get(2));
                rep->rot_m = 0.0;

            } else if (arguments.getNumArguments() == 2) {
                numRepetitions = parseMathExpression(arguments.get(0));
                rep->shiftx_m = 0.0;
                rep->shifty_m = 0.0;
                rep->rot_m = parseMathExpression(arguments.get(1));
            } else {
                std::cout << "Repeat: number of arguments not supported" << std::endl;
                return false;
            }
        } catch (std::runtime_error &e) {
            std::cout << e.what() << std::endl;
            return false;
        }

        if (numRepetitions < 0) {
            std::cout << "Repeat: a negative number of repetitions provided '"
                      << arguments.get(0) << " = " << rep->N_m << "'"
                      << std::endl;
            return false;
        }

        rep->N_m = numRepetitions;

        it += (arguments.getLengthConsumed() + 1);

        return true;
    }
}