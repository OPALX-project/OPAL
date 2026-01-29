//
// Class FieldWriter
//   This class writes the bunch internal fields on the grid to
//   file. It supports single core execution only.
//
// Copyright (c) 2020, Paul Scherrer Institut, Villigen PSI, Switzerland
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
#ifndef OPAL_FIELD_WRITER_H
#define OPAL_FIELD_WRITER_H

#include <string>

class FieldWriter
{
public:

    /// Dump a scalar or vector field to a file.
    /*
     * @param[in] field is the scalar or vector field on the grid
     * @param[in] name is the field name
     * @param[in] unit of the field
     * @param[in] step of the output
     * @param[in] image of the potential (optional)
     */
    template<typename FieldType>
    void dumpField(FieldType& field, std::string name,
                   std::string unit, long long step,
                   FieldType* image = nullptr);
};

#include "FieldWriter.hpp"

#endif