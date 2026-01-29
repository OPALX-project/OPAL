
// ------------------------------------------------------------------------
// $RCSfile: Object.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.4 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: Object
//   This abstract base class defines the common interface for all
//   objects defined in OPAL.
//
// ------------------------------------------------------------------------
//
// $Date: 2000/12/15 10:04:08 $
// $Author: opal $
//
// ------------------------------------------------------------------------

#include "AbstractObjects/Object.h"
#include "AbstractObjects/Attribute.h"
#include "AbstractObjects/Invalidator.h"
#include "AbstractObjects/Expressions.h"
#include "Parser/Statement.h"
#include "Utilities/Options.h"
#include "Utilities/ParseError.h"

#include <cmath>
#include <iostream>
#include <vector>

extern Inform *gmsg;


// Class Object
// ------------------------------------------------------------------------

Object::~Object() {
    // Invalidate all references to this object.
    for(std::set<Invalidator *>::iterator i = references.begin();
        i != references.end(); ++i) {
        (*i)->invalidate();
    }
}


bool Object::canReplaceBy(Object *) {
    // Default action: no replacement allowed.
    return false;
}


void Object::copyAttributes(const Object &source) {
    itsAttr = source.itsAttr;
}


void Object::execute() {
    // Default action: do nothing.
}


Attribute *Object::findAttribute(const std::string &name) {
    for(std::vector<Attribute>::iterator i = itsAttr.begin();
        i != itsAttr.end(); ++i) {
        if(i->getName() == name) return &(*i);
    }

    return 0;
}


const Attribute *Object::findAttribute(const std::string &name) const {
    for(std::vector<Attribute>::const_iterator i = itsAttr.begin();
        i != itsAttr.end(); ++i) {
        if(i->getName() == name) return &(*i);
    }

    return 0;
}


Object *Object::makeTemplate
(const std::string &name, TokenStream &, Statement &) {
    throw ParseError("Object::makeTemplate()", "Object \"" + name +
                     "\" cannot be used to define a macro.");
}


Object *Object::makeInstance(const std::string &/*name*/, Statement &, const Parser *) {
    throw ParseError("Object::makeInstance()", "Object \"" + getOpalName() +
                     "\" cannot be called as a macro.");
}


void Object::parse(Statement &stat) {
    while(stat.delimiter(',')) {
        std::string name = Expressions::parseString(stat, "Attribute name expected.");

        if(Attribute *attr = findAttribute(name)) {
            if(stat.delimiter('[')) {
                int index = int(std::round(Expressions::parseRealConst(stat)));
                Expressions::parseDelimiter(stat, ']');

                if(stat.delimiter('=')) {
                    attr->parseComponent(stat, true, index);
                } else if(stat.delimiter(":=")) {
                    attr->parseComponent(stat, false, index);
                } else {
                    throw ParseError("Object::parse()",
                                     "Delimiter \"=\" or \":=\" expected.");
                }
            } else if(stat.delimiter('=')) {
                attr->parse(stat, true);
            } else if(stat.delimiter(":=")) {
                attr->parse(stat, false);
            } else {
                attr->setDefault();
            }
        } else {
            throw ParseError("Object::parse()", "Object \"" + getOpalName() +
                             "\" has no attribute \"" + name + "\".");
        }
    }
}


void Object::parseShortcut(Statement &stat, bool eval) {
    // Only one attribute.
    if(stat.delimiter(',')) {
        stat.mark();
        std::string name;

        if(stat.word(name)) {
            if(stat.delimiter('=')) {
                if(Attribute *attr = findAttribute(name)) {
                    attr->parse(stat, eval);
                    return;
                } else {
                    throw ParseError("Object::parseShortcut()", "Object \"" + getOpalName() +
                                     "\" has no attribute \"" + name + "\".");
                }
            } else if(stat.delimiter(":=")) {
                if(Attribute *attr = findAttribute(name)) {
                    attr->parse(stat, false);
                    return;
                } else {
                    throw ParseError("Object::parseShortcut()", "Object \"" + getOpalName() +
                                     "\" has no attribute \"" + name + "\".");
                }
            }
        }

        stat.restore();
        itsAttr[0].parse(stat, false);
    }
}


void Object::print(std::ostream & msg) const {
    std::string head = getOpalName();
    Object *parent = getParent();
    if(parent != 0  &&  ! parent->getOpalName().empty()) {
        if(! getOpalName().empty()) head += ':';
        head += parent->getOpalName();
    }

    msg << head;
    int pos = head.length();

    for(std::vector<Attribute>::const_iterator i = itsAttr.begin();
        i != itsAttr.end(); ++i) {
        if(*i) i->print(pos);
    }
    msg << ';';
    msg << std::endl;
    return;
}


void Object::registerReference(Invalidator *ref) {
    references.insert(ref);
}


void Object::unregisterReference(Invalidator *ref) {
    references.erase(ref);
}

void Object::registerOwnership(const AttributeHandler::OwnerType &itsClass) const {
    if (getParent() != 0) return;

    const unsigned int end = itsAttr.size();
    const std::string name = getOpalName();
    for (unsigned int i = 0; i < end; ++ i) {
        AttributeHandler::addAttributeOwner(name, itsClass, itsAttr[i].getName());
    }
}

void Object::printHelp(std::ostream &/*os*/) const {
    *gmsg << endl << itsHelp << endl;

    if(itsAttr.size() > 0) {
        *gmsg << "Attributes:" << endl;

        size_t maxLength = 16;
        std::vector<Attribute>::const_iterator it;
        for (it = itsAttr.begin(); it != itsAttr.end(); ++ it) {
            std::string name = it->getName();
            maxLength = std::max(maxLength, name.length() + 1);
        }

        for (it = itsAttr.begin(); it != itsAttr.end(); ++ it) {
            std::string type = it->getType();
            std::string name = it->getName();
            *gmsg << '\t' << type << std::string(16 - type.length(), ' ');
            *gmsg << name << std::string(maxLength - name.length(), ' ');
            *gmsg << it->getHelp();
            if(it->isReadOnly()) *gmsg << " (read only)";
            *gmsg << endl;
        }
    }

    *gmsg << endl;
}


void Object::replace(Object *, Object *) {
    // Default action: do nothing.
}


void Object::update() {
    // Default action: do nothing.
}


bool Object::isBuiltin() const {
    return builtin;
}


bool Object::isShared() const {
    return sharedFlag;
}


void Object::setShared(bool flag) {
    sharedFlag = flag;
}


void Object::setDirty(bool dirty) {
    // The object is now different from the data base.
    modified = dirty;
}


bool Object::isDirty() const {
    return modified;
}


void Object::setFlag(bool flag) {
    flagged = flag;
}


bool Object::isFlagged() const {
    return flagged;
}

const Object *Object::getBaseObject() const {
    const Object *base = this;
    while(base->itsParent != 0) base = base->itsParent;
    return base;
}


const std::string &Object::getOpalName() const {
    return itsName;
}


Object *Object::getParent() const {
    return itsParent;
}


bool Object::isTreeMember(const Object *classObject) const {
    const Object *object = this;

    while(object != 0  &&  object != classObject) {
        object = object->itsParent;
    }

    return object != 0;
}


void Object::setOpalName(const std::string &name) {
    itsName = name;
}


void Object::setParent(Object *parent) {
    itsParent = parent;
}


void Object::clear() {
    occurrence = 0;
}


int Object::increment() {
    return ++occurrence;
}


int Object::occurrenceCount() {
    return occurrence;
}


Object::Object(int size, const char *name, const char *help):
    RCObject(), itsAttr(size), itsParent(0),
    itsName(name), itsHelp(help), occurrence(0), sharedFlag(false) {
    // NOTE: The derived classes must define the attribute handlers and
    //       any initial values for the attributes.

    // The object is an exemplar and must not be saved in the data base.
    builtin = true;
    flagged = modified = false;
}


Object::Object(const std::string &name, Object *parent):
    RCObject(), itsAttr(parent->itsAttr), itsParent(parent),
    itsName(name), itsHelp(parent->itsHelp), occurrence(0), sharedFlag(false) {
    // The object is now different from the data base.
    builtin = flagged = false;
    modified = true;
}

std::ostream &operator<<(std::ostream &os, const Object &object) {
    object.print(os);
    return os;
}