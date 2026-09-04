#include "Utilities/MSLang/SymmetricDifference.h"

#include <memory>
#include <regex>
#include <string>
#include <vector>
#include <iostream>

namespace mslang {
    void SymmetricDifference::print(int indentwidth) {
        std::string indent(' ', indentwidth);
        std::cout << indent << "Symmetric division\n"
                  << indent << "    first operand\n";
        firstOperand_m->print(indentwidth + 8);

        std::cout << indent << "    second operand\n";
        secondOperand_m->print(indentwidth + 8);
    }

    void SymmetricDifference::apply(std::vector<std::shared_ptr<Base> > &bfuncs) {
        std::vector<std::shared_ptr<Base> > first, second;

        firstOperand_m->apply(first);
        secondOperand_m->apply(second);
        for (auto item: first) {
            item->divideBy(second);
            bfuncs.emplace_back(item->clone());
        }

        for (auto item: first)
            item.reset();

        first.clear();

        firstOperand_m->apply(first);
        for (auto item: second) {
            item->divideBy(first);
            bfuncs.emplace_back(item->clone());
        }

        for (auto item: first)
            item.reset();
        for (auto item: second)
            item.reset();
    }

    bool SymmetricDifference::parse_detail(iterator &it, const iterator &end, Function* &fun) {
        SymmetricDifference *dif = static_cast<SymmetricDifference*>(fun);
        if (!parse(it, end, dif->firstOperand_m)) return false;

        std::regex argumentList("(,[a-z]+\\(.*)");
        std::regex endParenthesis("\\)(.*)");
        std::smatch what;

        std::string str(it, end);
        if (!std::regex_match(str, what, argumentList)) return false;

        iterator it2 = it + 1;
        if (!parse(it2, end, dif->secondOperand_m)) return false;

        it = it2;
        str = std::string(it, end);
        if (!std::regex_match(str, what, endParenthesis)) return false;

        std::string fullMatch = what[0];
        std::string rest = what[1];

        it += (fullMatch.size() - rest.size());

        return true;
    }
}
