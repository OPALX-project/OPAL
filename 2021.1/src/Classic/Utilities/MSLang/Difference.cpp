#include "Utilities/MSLang/Difference.h"

#include <boost/regex.hpp>

namespace mslang {
    void Difference::print(int indentwidth) {
        std::string indent(' ', indentwidth);
        std::cout << indent << "Difference\n"
                  << indent << "    nominators\n";
        dividend_m->print(indentwidth + 8);

        std::cout << indent << "    denominators\n";
        divisor_m->print(indentwidth + 8);
    }

    void Difference::apply(std::vector<std::shared_ptr<Base> > &bfuncs) {
        std::vector<std::shared_ptr<Base> > nom, denom;

        dividend_m->apply(nom);
        divisor_m->apply(denom);
        for (auto item: nom) {
            item->divideBy(denom);
            bfuncs.emplace_back(item->clone());
        }

        for (auto item: nom)
            item.reset();
        for (auto item: denom)
            item.reset();
    }

    bool Difference::parse_detail(iterator &it, const iterator &end, Function* &fun) {
        Difference *dif = static_cast<Difference*>(fun);
        if (!parse(it, end, dif->dividend_m)) return false;

        boost::regex argumentList("(,[a-z]+\\(.*)");
        boost::regex endParenthesis("\\)(.*)");
        boost::smatch what;

        std::string str(it, end);
        if (!boost::regex_match(str, what, argumentList)) return false;

        iterator it2 = it + 1;
        if (!parse(it2, end, dif->divisor_m)) return false;

        it = it2;
        str = std::string(it, end);
        if (!boost::regex_match(str, what, endParenthesis)) return false;

        std::string fullMatch = what[0];
        std::string rest = what[1];

        it += (fullMatch.size() - rest.size());

        return true;
    }
}
