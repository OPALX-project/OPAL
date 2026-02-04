#ifndef FUNCTION_HPP
#define FUNCTION_HPP

#include <functional>
#include <string>
#include <tuple>
#include <variant>
#include <vector>

namespace client {
    namespace function {
        typedef std::variant<double, bool, std::string> argument_t;

        typedef std::vector<argument_t> arguments_t;

        typedef std::function<std::tuple<double, bool> (arguments_t)> type;

        typedef std::pair<std::string, type> named_t;
    }
}

#endif
