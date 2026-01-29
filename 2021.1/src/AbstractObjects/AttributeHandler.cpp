// ------------------------------------------------------------------------
// $RCSfile: AttributeHandler.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: AttributeHandler
//   An abstract class used to parse and print attributes.
//
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:33:34 $
// $Author: Andreas Adelmann $
//
// ------------------------------------------------------------------------

#include "AbstractObjects/AttributeHandler.h"
#include "Parser/Statement.h"
#include "Utilities/OpalException.h"

#include <boost/algorithm/string/join.hpp>

// Class AttributeHandler
// ------------------------------------------------------------------------

std::multimap<std::string,
              std::pair<AttributeHandler::OwnerType,
                        std::string>
              > AttributeHandler::attributeOwnerDictionary_s;

AttributeHandler::AttributeHandler
(const std::string &name, const std::string &help, AttributeBase *def):
    RCObject(), itsName(name), itsHelp(help), itsDefault(def),
    is_deferred(false), is_readonly(false)
{}


AttributeHandler::~AttributeHandler()
{}


AttributeHandler *AttributeHandler::clone() const {
    throw OpalException("AttributeHandler::clone()",
                        "Internal error: should not call this method.");
}


AttributeBase *AttributeHandler::getDefault() const {
    if(itsDefault.isValid()) {
        return &*itsDefault;
    } else {
        throw OpalException("AttributeHandler::getDefault()",
                            "Attribute \"" + itsName + "\" has no default value.");
    }
}


const std::string &AttributeHandler::getHelp() const {
    return itsHelp;
}


const std::string &AttributeHandler::getName() const {
    return itsName;
}


void AttributeHandler::parseComponent
(Attribute &, Statement &, bool, int) const {
    // Default behaviour.
    throw OpalException("AttributeHandler::parseComponent()",
                        "You cannot assign to a component of \"" + itsName +
                        "\" which is not a vector value.");
}


bool AttributeHandler::isDeferred() const {
    return is_deferred;
}


void AttributeHandler::setDeferred(bool flag) {
    is_deferred = flag;
}


bool AttributeHandler::isReadOnly() const {
    return is_readonly;
}


void AttributeHandler::setReadOnly(bool flag) {
    is_readonly = flag;
}

void AttributeHandler::setPredefinedValues(const std::set<std::string>& predefinedValues,
                                           const std::string& defaultValue)
{
    std::string validValues = " Valid values are " + boost::algorithm::join(predefinedValues, ", ") + ".";
    if (defaultValue != "_HAS_NO_DEFAULT_") {
        validValues += " Its default is " + defaultValue + ".";
    }
    std::string* help = const_cast<std::string*>(&itsHelp);
    *help += validValues;
}

std::multimap<AttributeHandler::OwnerType, std::string> AttributeHandler::getOwner(const std::string &att) {
    std::multimap<OwnerType, std::string> possibleOwners;

    if (attributeOwnerDictionary_s.find(att) != attributeOwnerDictionary_s.end()) {
        auto its = attributeOwnerDictionary_s.equal_range(att);

        for (auto it = its.first; it != its.second; ++ it) {
            auto owner = it->second;

            possibleOwners.insert(std::make_pair(owner.first, owner.second));
        }
    }

    return possibleOwners;
}

void AttributeHandler::addAttributeOwner(const std::string &owner,
                                         const AttributeHandler::OwnerType &type,
                                         const std::string &name) {
    attributeOwnerDictionary_s.insert(std::make_pair(name,
                                                     std::make_pair(type,
                                                                    owner)
                                                     )
                                      );
}