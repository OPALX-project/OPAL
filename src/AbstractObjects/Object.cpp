
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
#include "Utilities/Round.h"

#include <iostream>
//using std::cerr;
//using std::endl;
using std::vector;
using namespace std;

extern Inform *gmsg;


// Class Object
// ------------------------------------------------------------------------

Object::~Object() {
    // Invalidate all references to this object.
    for(set<Invalidator *>::iterator i = references.begin();
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


Attribute *Object::findAttribute(const string &name) {
    for(vector<Attribute>::iterator i = itsAttr.begin();
        i != itsAttr.end(); ++i) {
        if(i->getName() == name) return &(*i);
    }

    return 0;
}


const Attribute *Object::findAttribute(const string &name) const {
    for(vector<Attribute>::const_iterator i = itsAttr.begin();
        i != itsAttr.end(); ++i) {
        if(i->getName() == name) return &(*i);
    }

    return 0;
}


Object *Object::makeTemplate
(const string &name, TokenStream &, Statement &) {
    throw ParseError("Object::makeTemplate()", "Object \"" + name +
                     "\" cannot be used to define a macro.");
}


Object *Object::makeInstance(const string &/*name*/, Statement &, const Parser *) {
    throw ParseError("Object::makeInstance()", "Object \"" + getOpalName() +
                     "\" cannot be called as a macro.");
}


void Object::parse(Statement &stat) {
    while(stat.delimiter(',')) {
        string name = Expressions::parseString(stat, "Attribute name expected.");

        if(Attribute *attr = findAttribute(name)) {
            if(stat.delimiter('[')) {
                int index = int(Round(Expressions::parseRealConst(stat)));
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


void Object::parseShortcut(Statement &stat) {
    // Only one attribute.
    if(stat.delimiter(',')) {
        stat.mark();
        string name;

        if(stat.word(name)) {
            if(stat.delimiter('=')) {
                if(Attribute *attr = findAttribute(name)) {
                    attr->parse(stat, true);
                    return;
                } else {
                    throw ParseError("Object::parse()", "Object \"" + getOpalName() +
                                     "\" has no attribute \"" + name + "\".");
                }
            } else if(stat.delimiter(":=")) {
                if(Attribute *attr = findAttribute(name)) {
                    attr->parse(stat, false);
                    return;
                } else {
                    throw ParseError("Object::parse()", "Object \"" + getOpalName() +
                                     "\" has no attribute \"" + name + "\".");
                }
            }
        }

        stat.restore();
        itsAttr[0].parse(stat, false);
    }
}


void Object::print(std::ostream & msg) const {
    string head = getOpalName();
    Object *parent = getParent();
    if(parent != 0  &&  ! parent->getOpalName().empty()) {
        if(! getOpalName().empty()) head += ':';
        head += parent->getOpalName();
    }

    msg << head;
    int pos = head.length();

    for(vector<Attribute>::const_iterator i = itsAttr.begin();
        i != itsAttr.end(); ++i) {
        if(*i) i->print(pos);
    }
    msg << ';';
    msg << endl;
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
        vector<Attribute>::const_iterator it;
        for (it = itsAttr.begin(); it != itsAttr.end(); ++ it) {
            string name = it->getName();
            maxLength = std::max(maxLength, name.length() + 1);
        }

        for (it = itsAttr.begin(); it != itsAttr.end(); ++ it) {
            string type = it->getType();
            string name = it->getName();
            *gmsg << '\t' << type << string(16 - type.length(), ' ');
            *gmsg << name << string(maxLength - name.length(), ' ');
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


const string &Object::getOpalName() const {
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


void Object::setOpalName(const string &name) {
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


Object::Object(const string &name, Object *parent):
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