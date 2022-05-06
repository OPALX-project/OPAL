//
// Copyright (c) 2008 - 2020, Paul Scherrer Institut, Villigen PSI, Switzerland
//
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
#ifndef GSLERRORHANDLING_H
#define GSLERRORHANDLING_H

#include "Utilities/OpalException.h"

inline void handleGSLErrors(const char *reason,
                            const char *file,
                            int /*line*/,
                            int /*gsl_errno*/) {
    throw OpalException(file, reason);
}

#endif