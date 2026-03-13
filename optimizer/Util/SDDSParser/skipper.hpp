//
// Struct skipper
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
#ifndef SKIPPER_HPP_
#define SKIPPER_HPP_

// Skipper functionality is now handled by SimpleParser::skipWhitespace()
// This file kept for compatibility but no longer needed
/// \todo Remove this compatibility header once all downstream users migrate.

namespace SDDS { namespace parser
{
    // Skipper is now integrated into SimpleParser

}}

#endif


