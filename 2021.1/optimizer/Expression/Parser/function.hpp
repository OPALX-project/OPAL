#ifndef FUNCTION_HPP
#define FUNCTION_HPP

#include <vector>

#include "boost/function.hpp"
#include "boost/tuple/tuple.hpp"
#include "boost/variant/variant.hpp"

namespace client { namespace function {

    typedef boost::variant <
          double
        , bool
        , std::string
        >
    argument_t;

    typedef std::vector<argument_t> arguments_t;

    typedef boost::function<boost::tuple<double, bool> (arguments_t)> type;

    typedef std::pair<std::string, type> named_t;

}}

#endif
