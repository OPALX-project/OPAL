//
// Class SequenceTemplate
//
//   An ``archetype'' for a SEQUENCE with arguments.
//   The model is stored in form of a MacroStream.  A call to the macro
//   sequence is expanded by first replacing the arguments, and then parsing
//   the resulting stream as a SEQUENCE definition.
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

#ifndef OPAL_SequenceTemplate_HH
#define OPAL_SequenceTemplate_HH

#include "OpalParser/Macro.h"
#include "Parser/SimpleStatement.h"
#include "OpalParser/MacroStream.h"
#include <list>

class Sequence;
class Statement;
class TokenStream;

class SequenceTemplate: public Macro {

    friend class Sequence;

public:

    SequenceTemplate();
    virtual ~SequenceTemplate();

    /// Make clone.
    //  Throw OpalException, since the template cannot be cloned.
    virtual SequenceTemplate *clone(const std::string &name);

    /// Make line instance.
    //  The instance gets the name [b]name[/b], and its actual arguments
    //  are read from [b]stat[/b].  The parser is ignored.
    virtual Object *makeInstance
    (const std::string &name, Statement &, const Parser *);

    /// Make a sequence template.
    //  Return nullptr, since one cannot make a template from a template.
    virtual Object *makeTemplate(const std::string &name, TokenStream &, Statement &);

    /// Parse the sequence template.
    void parseTemplate(TokenStream &, Statement &);

private:

    // Not implemented.
    SequenceTemplate(const SequenceTemplate &);
    void operator=(const SequenceTemplate &);

    // Clone constructor.
    SequenceTemplate(const std::string &name, Object *parent);

    // The contained beam sequence element list.
    MacroStream body;
};

#endif // OPAL_SequenceTemplate_HH
