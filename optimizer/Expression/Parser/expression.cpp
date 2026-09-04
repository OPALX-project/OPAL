#include <string>

#include "expression_def.hpp"

typedef std::string::const_iterator iterator_type;
template struct client::parser::expression<iterator_type>;