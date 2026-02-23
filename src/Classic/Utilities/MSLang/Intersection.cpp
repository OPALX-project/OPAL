#include "Utilities/MSLang/Intersection.h"

#include <iostream>
#include <memory>
#include <regex>
#include <string>
#include <vector>

namespace mslang {
    void Intersection::print(int indentwidth) {
        std::string indent(' ', indentwidth);
        std::cout << indent << "Intersection\n"
                  << indent << "    first operand\n";
        firstOperand_m->print(indentwidth + 8);

        std::cout << indent << "    second operand\n";
        secondOperand_m->print(indentwidth + 8);
    }

    void Intersection::apply(std::vector<std::shared_ptr<Base> > &bfuncs) {
        std::vector<std::shared_ptr<Base> > first, firstrep, second;

        firstOperand_m->apply(first);
        firstOperand_m->apply(firstrep);
        secondOperand_m->apply(second);

        for (auto item: firstrep) {
            item->divideBy(second);
        }

        for (auto item: first) {
            item->divideBy(firstrep);
            bfuncs.emplace_back(item->clone());
        }

        for (auto item: first)
            item.reset();
        for (auto item: firstrep)
            item.reset();
        for (auto item: second)
            item.reset();
    }

    bool Intersection::parse_detail(iterator &it, const iterator &end, Function* &fun) {
        Intersection *inter = static_cast<Intersection*>(fun);
        if (!parse(it, end, inter->firstOperand_m)) return false;

        std::regex argumentList("(,[a-z]+\\(.*)");
        std::regex endParenthesis("\\)(.*)");
        std::smatch what;

        std::string str(it, end);
        if (!std::regex_match(str, what, argumentList)) return false;

        iterator it2 = it + 1;
        if (!parse(it2, end, inter->secondOperand_m)) return false;

        it = it2;
        str = std::string(it, end);
        if (!std::regex_match(str, what, endParenthesis)) return false;

        std::string fullMatch = what[0];
        std::string rest = what[1];

        it += (fullMatch.size() - rest.size());

        return true;
    }
}
