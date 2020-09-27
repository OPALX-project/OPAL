#ifndef OPAL_SRefAttr_HH
#define OPAL_SRefAttr_HH

// ------------------------------------------------------------------------
// $RCSfile: SRefAttr.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.2 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Template class SRefAttr<T>
//
// ------------------------------------------------------------------------
//
// $Date: 2000/04/07 12:02:54 $
// $Author: opal $
//
// ------------------------------------------------------------------------

#include "AbstractObjects/AttributeBase.h"
#include "AbstractObjects/OpalData.h"
#include "AbstractObjects/Object.h"
#include "Expressions/ADeferred.h"
#include "Expressions/SValue.h"
#include "Utilities/ParseError.h"
//#include "Utilities/Options.h"
#include "Utilities/Options.h"
#include <vector>


namespace Expressions {

    // Class SRefAttr
    // ----------------------------------------------------------------------

    /// An attribute defined as a reference to a scalar.
    //  The referred attribute may have values of type real, logical or string.
    //  When a reference is seen, the pointers to the relevant object and
    //  attribute are left zero.  When the expression value is required,
    //  the object and the attribute are searched, and the pointers cached.
    //  The reference is registered with the object.  If the object referred
    //  to is deleted, it calls the invalidate() method of all reference
    //  expressions referring to it.  This resets the pointers to zero, so
    //  that the next evaluation forces to search for a replacement object.

    template <class T>
    class SRefAttr: public AttributeBase {

    public:

        /// Constructor.
        //  Use object name [b]oName[/b] to identify the object containing
        //  the scalar, and [b]aName[/b] to identify the scalar itself.
        SRefAttr(const std::string &oName, const std::string &aName, int index);

        SRefAttr(const SRefAttr &);
        virtual ~SRefAttr();

        /// Make clone.
        virtual SRefAttr<T> *clone() const;

        /// Evaluate.
        //  Evaluate the reference and return the value.
        virtual T evaluate() const;

        /// Return real value.
        //  This function has been added for speed of access.
        virtual double getReal();

        /// Invalidate.
        //  Force re-evaluation of the reference.
        virtual void invalidate();

        /// Print the reference.
        virtual void print(std::ostream &) const;

        /// Store new value.
        //  Evaluate the reference and assign to the scalar referred to.
        virtual void set(const T &) const;

    private:

        // Not implemented.
        SRefAttr();
        void operator=(const SRefAttr &);

        // Fill in the reference.
        void fill() const;

        // The name of the type referred to.
        static const std::string typeName;

        // The referred object, attribute and index.
        const std::string obj_name;
        const std::string att_name;
        const int itsIndex;

        // The object and attribute referred to.
        mutable Object    *itsObject;
        mutable Attribute *itsAttr;
    };


    template <class T>
    inline std::ostream &operator<<(std::ostream &os, const SRefAttr<T> &a) {
        a.print(os);
        return os;
    }


    // Implementation of class SRefAttr<T>.
    // ------------------------------------------------------------------------

    template <class T>
    SRefAttr<T>::SRefAttr
    (const std::string &oName, const std::string &aName, int index):
        obj_name(oName), att_name(aName), itsIndex(index),
        itsObject(0), itsAttr(0)
    {}

    template <class T>
    SRefAttr<T>::SRefAttr(const SRefAttr &rhs):
        AttributeBase(),
        obj_name(rhs.obj_name), att_name(rhs.att_name), itsIndex(rhs.itsIndex),
        itsObject(rhs.itsObject), itsAttr(rhs.itsAttr)
    {}


    template <class T>
    SRefAttr<T>::~SRefAttr() {
        if(itsObject) itsObject->unregisterReference(this);
    }



    template <class T>
    SRefAttr<T> *SRefAttr<T>::clone() const {
        return new SRefAttr<T>(*this);
    }


    template <class T>
    T SRefAttr<T>::evaluate() const {
        fill();

        if(AttributeBase *base = &itsAttr->getBase()) {
            if(itsIndex) {
                if(ADeferred<T> *value = dynamic_cast<ADeferred<T>*>(base)) {
                    std::vector<T> array = value->evaluate();
                    if(itsIndex > int(array.size())) {
                        throw ParseError("SRefAttr::evaluate()", "Reference \"" +
                                         getImage() + "\" has index out of range.");
                    } else {
                        return array[itsIndex - 1];
                    }
                } else {
                    throw ParseError("SRefAttr::evaluate()", "Reference \"" +
                                     getImage() + "\" is not an array.");
                }
            } else {
                if(SValue<T> *value = dynamic_cast<SValue<T> *>(base)) {
                    return value->evaluate();
                } else {
                    throw ParseError("SRefAttr::evaluate()", getImage() +
                                     "\" is of the wrong type.");
                }
            }
        }

        return T(0);
    }


    template <class T>
    double SRefAttr<T>::getReal() {
        throw ParseError("SRefAttr<T>::getReal()",
                         "Attribute is not of real type.");
    }


    template <> inline
    double SRefAttr<double>::getReal() {
        return evaluate();
    }


    template <class T>
    void SRefAttr<T>::print(std::ostream &os) const {
        os << obj_name;
        if(! att_name.empty()) os << "->" << att_name;
        if(itsIndex != 0) os << '[' << itsIndex << ']';
        return;
    }


    template <class T>
    void SRefAttr<T>::invalidate() {
        itsObject = 0;
        itsAttr = 0;
    }


    template <class T>
    void SRefAttr<T>::set(const T &value) const {
        fill();

        if(AttributeBase *base = &itsAttr->getBase()) {
            if(dynamic_cast<SValue<T> *>(base)) {
                return itsAttr->set(new SValue<T>(value));
            } else {
                throw ParseError("Real::get()", "Attribute \"" +
                                 itsAttr->getName() + "\" is of the wrong type.");
            }
        }
    }


    template <class T>
    void SRefAttr<T>::fill() const {
        if(itsObject == 0) {
            itsObject = OpalData::getInstance()->find(obj_name);
            if(itsObject == 0) {
                if(att_name.empty()  &&  itsIndex <= 0) {
                    throw ParseError("SRefAttr::fill()",
                                     "\nThe <variable> \"" + obj_name + "\" is unknown.\n");
                } else {
                    throw ParseError("SRefAttr::fill()",
                                     "Object \"" + obj_name + "\" is unknown.");
                }
            }

            // Register the reference with the object, to allow invalidation
            // when the object is deleted.
            itsObject->registerReference(const_cast<SRefAttr<T>*>(this));

            if(att_name.empty()) {
                itsAttr = itsObject->findAttribute("VALUE");
                if(itsAttr == 0) {
                    throw ParseError("SRefAttr::fill()", "Object \"" + obj_name +
                                     "\" is not a variable, constant or vector.");
                }
            } else {
                itsAttr = itsObject->findAttribute(att_name);
                if(itsAttr == 0) {
                    throw ParseError("SRefAttr::fill()", "Attribute \"" + obj_name +
                                     "->" + att_name + "\" is unknown.");
                }
            }
        }
    }

}

#endif // OPAL_SRefAttr_HH