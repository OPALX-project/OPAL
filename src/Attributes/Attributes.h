#ifndef OPAL_Attributes_HH
#define OPAL_Attributes_HH

// ------------------------------------------------------------------------
// $RCSfile: Attributes.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Namespace Attributes
//
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:33:37 $
// $Author: Andreas Adelmann $
//
// ------------------------------------------------------------------------

#include "AbstractObjects/Attribute.h"
#include "Parser/Token.h"
#include <list>
#include <string>
#include <vector>

class PlaceRep;
class RangeRep;
class TableRowRep;

// Namespace Attributes.
// ------------------------------------------------------------------------
/// A collection of routines to construct object Attributes and retrieve
//  their values.
//  Each exemplar constructor must use these functions to fill in the
//  attribute handlers for all command attributes.
//  This namespace is part of the interface to the class Attribute,
//  used to store object attributes.

namespace Attributes {

    /// Make logical attribute.
    //  Initial value undefined.
    extern Attribute makeBool(const std::string &name, const std::string &help);

    /// Make logical attribute.
    //  Initial value defined.
    extern Attribute
    makeBool(const std::string &name, const std::string &help, bool initial);

    /// Return logical value.
    //  Evaluate any expression and return value.
    //  If the expression is non-deferred, the value is cached.
    //  Re-evaluation takes place only if a new definition has been read.
    bool getBool(const Attribute &);

    /// Set logical value.
    extern void setBool(Attribute &, bool);


    /// Create a logical array attribute.
    //  Initial value is empty array.
    extern Attribute makeBoolArray(const std::string &name, const std::string &help);

    /// Get logical array value.
    extern std::vector<bool> getBoolArray(const Attribute &);

    /// Set logical array value.
    extern void setBoolArray(Attribute &, const std::vector<bool> &);


    /// Create a place attribute.
    //  Initial value is "#S".
    extern Attribute makePlace(const std::string &name, const std::string &help);

    /// Get place value.
    extern PlaceRep getPlace(const Attribute &);

    /// Set place value.
    extern void setPlace(Attribute &, const PlaceRep &);


    /// Create a range attribute.
    //  Initial value is "FULL".
    extern Attribute makeRange(const std::string &name, const std::string &help);

    /// Get range value.
    extern RangeRep getRange(const Attribute &);

    /// Set range value.
    extern void setRange(Attribute &, const RangeRep &);


    /// Make real attribute.
    //  Initial value is undefined.
    extern Attribute
    makeReal(const std::string &name, const std::string &help);

    /// Make real attribute.
    //  Initial value is defined.
    extern Attribute
    makeReal(const std::string &name, const std::string &help, double initial);

    /// Return real value.
    extern double getReal(const Attribute &attr);

    /// Set real value.
    extern void setReal(Attribute &attr, double val);


    /// Create real array attribute.
    //  Initial value is empty array.
    extern Attribute makeRealArray(const std::string &name, const std::string &help);

    /// Get array value.
    extern std::vector<double> getRealArray(const Attribute &);

    /// Set array value.
    extern void setRealArray(Attribute &, const std::vector<double> &);


    /// Create a reference attribute.
    //  Initial value is undefined.
    extern Attribute makeReference(const std::string &name, const std::string &help);


    /// Make string attribute.
    //  Initial value undefined.
    extern Attribute makeString(const std::string &name, const std::string &help);

    /// Make string attribute.
    //  Initial value is defined.
    extern Attribute
    makeString(const std::string &name, const std::string &help, const std::string &initial);

    /// Get string value.
    extern std::string getString(const Attribute &attr);

    /// Set string value.
    extern void setString(Attribute &attr, const std::string &val);

    /// Make uppercase string attribute.
    //  Initial value undefined.
    extern Attribute makeUpperCaseString(const std::string &name, const std::string &help);

    /// Make uppercase string attribute.
    //  Initial value is defined.
    extern Attribute
    makeUpperCaseString(const std::string &name, const std::string &help, const std::string &initial);

    /// Set uppercase string value.
    extern void setUpperCaseString(Attribute &attr, const std::string &val);

    /// Create a string array attribute.
    //  Initial value is empty array.
    extern Attribute makeStringArray(const std::string &name, const std::string &help);

    /// Get string array value.
    extern std::vector<std::string> getStringArray(const Attribute &);

    /// Set string array value.
    extern void setStringArray(Attribute &, const std::vector<std::string> &);

    /// Make uppercase string array attribute.
    //  Initial value is empty array.
    extern Attribute makeUpperCaseStringArray(const std::string &name, const std::string &help);

    /// Set upper case string array value.
    extern void setUpperCaseStringArray(Attribute &, const std::vector<std::string> &);

    /// Create a table row attribute.
    //  Initial value is undefined.
    extern Attribute makeTableRow(const std::string &name, const std::string &help);

    /// Get table row value.
    extern TableRowRep getTableRow(const Attribute &);

    /// Set table row value.
    extern void setTableRow(Attribute &, const TableRowRep &);


    /// Make token list attribute.
    //  Initial value is empty list.
    extern Attribute makeTokenList(const std::string &name, const std::string &help);

    /// Return token list value.
    extern std::list<Token> getTokenList(const Attribute &attr);

    /// Set token list value.
    extern void setTokenList(Attribute &attr, const std::list<Token> &);


    /// Make token list attribute.
    //  Initial value is empty array.
    extern Attribute makeTokenListArray(const std::string &name, const std::string &help);

    /// Return token list array value.
    extern std::vector<std::list<Token> >
    getTokenListArray(const Attribute &attr);

    /// Set token list array value.
    extern void
    setTokenListArray(Attribute &attr, const std::vector<std::list<Token> > &);

};

#endif // OPAL_Attributes_HH