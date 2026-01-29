#ifndef MSLANG_ARGUMENTEXTRACTOR_H
#define MSLANG_ARGUMENTEXTRACTOR_H

#include <string>
#include <vector>
#include <utility>

namespace mslang {
    struct ArgumentExtractor {
        ArgumentExtractor(const std::string &input);

        std::string get(unsigned int i) const;
        unsigned int getLengthConsumed() const;
        unsigned int getNumArguments() const;

        std::vector<std::pair<unsigned int, unsigned int>> argumentBoundaries_m;
        std::string inputArguments_m;
        unsigned int lengthConsumed_m;
    };
}
#endif