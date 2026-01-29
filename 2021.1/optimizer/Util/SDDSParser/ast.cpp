//
// Namespace ast
//
// Copyright (c) 2015, Christof Metzger-Kraus, Helmholtz-Zentrum Berlin
// All rights reserved
//
// This file is part of OPAL.
//
// OPAL is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// You should have received a copy of the GNU General Public License
// along with OPAL. If not, see <https://www.gnu.org/licenses/>.
//
#include "ast_def.hpp"
#include "skipper.hpp"

typedef std::string::const_iterator iterator_t;
typedef SDDS::parser::skipper<iterator_t> skipper_t;
template struct SDDS::parser::string<iterator_t, skipper_t>;
template struct SDDS::parser::qstring<iterator_t, skipper_t>;
