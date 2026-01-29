#include "Utilities/MSLang/matheval.hpp"

namespace mslang {
    double parseMathExpression(const std::string &str) {
        return matheval::parse<double>(str, {});
    }
}