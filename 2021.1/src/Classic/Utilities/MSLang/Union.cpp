#include "Utilities/MSLang/Union.h"

#include <boost/regex.hpp>

namespace mslang {
    void Union::print(int indentwidth) {
        std::string indent(indentwidth, ' ');
        std::string indent2(indentwidth + 8, ' ');
        std::string indent3(indentwidth + 16, ' ');
        std::cout << indent << "union, " << std::endl;
        std::cout << indent2 << "funcs: {\n";
        funcs_m.front()->print(indentwidth + 16);
        for (unsigned int i = 1; i < funcs_m.size(); ++ i) {
            std::cout << "\n"
                      << indent3 << "," << std::endl;
            funcs_m[i]->print(indentwidth + 16);
        }
        std::cout << "\n"
                  << indent2 << "} ";
    }

    void Union::apply(std::vector<std::shared_ptr<Base> > &bfuncs) {
        for (unsigned int i = 0; i < funcs_m.size(); ++ i) {
            std::vector<std::shared_ptr<Base> > children;
            Function *func = funcs_m[i];
            func->apply(children);
            bfuncs.insert(bfuncs.end(), children.begin(), children.end());
        }
    }

    bool Union::parse_detail(iterator &it, const iterator &end, Function* &fun) {
        Union *unin = static_cast<Union*>(fun);
        unin->funcs_m.push_back(NULL);
        if (!parse(it, end, unin->funcs_m.back())) return false;

        boost::regex argumentList("(,[a-z]+\\(.*)");
        boost::regex endParenthesis("\\)(.*)");
        boost::smatch what;

        std::string str(it, end);
        while (boost::regex_match(str, what, argumentList)) {
            iterator it2 = it + 1;
            unin->funcs_m.push_back(NULL);

            if (!parse(it2, end, unin->funcs_m.back())) return false;

            it = it2;
            str = std::string(it, end);
        }

        str = std::string(it, end);
        if (!boost::regex_match(str, what, endParenthesis)) return false;

        std::string fullMatch = what[0];
        std::string rest = what[1];

        it += (fullMatch.size() - rest.size());

        return true;
    }
}