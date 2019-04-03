// ------------------------------------------------------------------------
// $RCSfile: Attribute.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.3.4.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class Attribute:
//   Interface a polymorphic attribute, including its parser.
//
// ------------------------------------------------------------------------
//
// $Date: 2002/12/09 15:06:07 $
// $Author: jsberg $
//
// ------------------------------------------------------------------------

#include "AbstractObjects/Attribute.h"
#include "AbstractObjects/AttributeBase.h"
#include "Utilities/Options.h"

#include <set>
#include <iostream>

extern Inform *gmsg;


// Class Attribute
// ------------------------------------------------------------------------

Attribute::Attribute():
    base(nullptr),
    handler(),
    isDefault(true)
{}


Attribute::Attribute(const Attribute &rhs):
    base(rhs.base),
    handler(rhs.handler),
    isDefault(rhs.isDefault)
{}


Attribute::Attribute(AttributeHandler *h, AttributeBase *b):
    base(b),
    handler(h),
    isDefault(true)
{}


Attribute::~Attribute()
{}


const Attribute &Attribute::operator=(const Attribute &rhs) {
    if(&rhs != this) {
        base = rhs.base;
        handler = rhs.handler;
    }

    return *this;
}


AttributeBase &Attribute::getBase() const {
    return *base;
}

bool Attribute::isBaseAllocated() const {
    return base != nullptr;
}

AttributeHandler &Attribute::getHandler() const {
    return *handler;
}


const std::string &Attribute::getHelp() const {
    return handler->getHelp();
}


std::string Attribute::getImage() const {
    return base->getImage();
}


const std::string &Attribute::getName() const {
    return handler->getName();
}


const std::string &Attribute::getType() const {
    return handler->getType();
}


bool Attribute::isDeferred() const {
    return handler->isDeferred();
}


void Attribute::setDeferred(bool flag) {
    handler->setDeferred(flag);
}


bool Attribute::isExpression() const {
    return base->isExpression();
}


bool Attribute::isReadOnly() const {
    return handler->isReadOnly();
}


void Attribute::setReadOnly(bool flag) {
    handler->setReadOnly(flag);
}


void Attribute::parse(Statement &stat, bool eval) {
    handler->parse(*this, stat, eval);
    isDefault = false;
}


void Attribute::parseComponent(Statement &stat, bool eval, int index) {
    handler->parseComponent(*this, stat, eval, index);
    isDefault = false;
}


void Attribute::set(AttributeBase *newBase) {
    base = newBase;
    isDefault = false;
}


void Attribute::setDefault() {
    base = handler->getDefault();
    isDefault = true;
}


void Attribute::print(int &pos) const {
    if(*this) {
        std::string name  = getName();
        std::string image = getImage();
        int step = name.length() + image.length() + 2;
        pos += step;

        if(pos > 74) {
            *gmsg << ',' << endl << "  ";
            pos = step + 2;
        } else {
            *gmsg << ',';
        }
        *gmsg << name << (isExpression() ? ":=" : "=") << image;
    }
}


std::ostream &operator<<(std::ostream &os, const Attribute &attr) {
    if(attr) {
        attr.getBase().print(os);
    } else {
        os << "*undefined*";
    }

    return os;
}