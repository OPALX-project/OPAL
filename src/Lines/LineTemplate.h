//
// Class LineTemplate
//
//   An ``archetype'' for a OPAL beam line with arguments.
//   The model is stored in form of a MacroStream.  A call to the macro line
//   is expanded by first replacing the arguments, and then parsing the
//   resulting stream as a LINE definition.
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

#ifndef OPAL_LineTemplate_HH
#define OPAL_LineTemplate_HH

#include "OpalParser/Macro.h"
#include "OpalParser/MacroStream.h"

class Line;
class Statement;
class TokenStream;

class LineTemplate: public Macro {

    friend class Line;

public:

    LineTemplate();
    virtual ~LineTemplate();

    /// Make clone.
    //  Throw OpalException, since the template cannot be cloned.
    virtual LineTemplate *clone(const std::string &name);

    /// Make line instance.
    //  The instance gets the name [b]name[/b], and its actual arguments
    //  are read from [b]stat[/b]. The parser is ignored.
    virtual Object *makeInstance
    (const std::string &name, Statement &stat, const Parser *);

    /// Make a line template.
    //  Return nullptr, since one cannot make a template from a template.
    virtual Object *makeTemplate(const std::string &, TokenStream &, Statement &);

    /// Parse the line template.
    //  Read the actual arguments from [b]stat[/b]. [b]is[/b] is not used.
    void parseTemplate(TokenStream &is, Statement &stat);

private:

    // Not implemented.
    LineTemplate(const LineTemplate &);
    void operator=(const LineTemplate &);

    // Clone constructor.
    LineTemplate(const std::string &name, Object *parent);

    // The contained beam line element list.
    MacroStream body;
};

#endif // OPAL_LineTemplate_HH
