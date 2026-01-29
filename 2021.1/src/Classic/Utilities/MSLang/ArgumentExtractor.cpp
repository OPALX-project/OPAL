#include "Utilities/MSLang/ArgumentExtractor.h"

#include <stdexcept>

namespace mslang {
    ArgumentExtractor::ArgumentExtractor(const std::string &input):
        inputArguments_m(input)
    {
        const unsigned int length = input.length();
        unsigned int pos = 0;
        unsigned short level = 1;
        unsigned int start = 0;

        for (; pos < length; ++ pos) {
            if (input[pos] == ')') {
                -- level;
                if (level == 0) {
                    argumentBoundaries_m.push_back(std::make_pair(start, pos - start));
                    break;
                }
            }
            else if (input[pos] == '(')
                ++ level;
            else if ((input[pos] == ','
                      || input[pos] == ';')
                     && level == 1) {
                argumentBoundaries_m.push_back(std::make_pair(start, pos - start));
                start = pos + 1;
            }
        }

        inputArguments_m = inputArguments_m.substr(0, pos);
        lengthConsumed_m = pos;
    }

    std::string ArgumentExtractor::get(unsigned int i) const {
        if (i >= argumentBoundaries_m.size())
            throw std::out_of_range("function only has " + std::to_string(argumentBoundaries_m.size()) + " arguments");

        auto boundaries = argumentBoundaries_m[i];

        return inputArguments_m.substr(boundaries.first, boundaries.second);
    }

    unsigned int ArgumentExtractor::getLengthConsumed() const {
        return lengthConsumed_m;
    }

    unsigned int ArgumentExtractor::getNumArguments() const {
        return argumentBoundaries_m.size();
    }

}