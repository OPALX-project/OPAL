#ifndef OPAL_SRefExpr_HH
#define OPAL_SRefExpr_HH

// ------------------------------------------------------------------------
// $RCSfile: SRefExpr.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.4.4.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class SRefExpr:
//
// ------------------------------------------------------------------------
//
// $Date: 2002/12/09 15:06:08 $
// $Author: jsberg $
//
// ------------------------------------------------------------------------

#include "AbstractObjects/Expressions.h"
#include "AbstractObjects/Invalidator.h"
#include "AbstractObjects/OpalData.h"
#include "Expressions/SValue.h"
#include "Utilities/ParseError.h"
#include "Utilities/Options.h"
#include <iosfwd>
#include <sstream>


namespace Expressions {

    // Class SRefExpr
    // ----------------------------------------------------------------------
    /// An expression defined as a reference to a scalar.
    //  The referred attribute may have values of type real, logical or string.
    //  Reference expressions are (re)evaluated as needed as part of the
    //  containing expression.

    template <class T>
    class SRefExpr: public Scalar<T>, public Invalidator {

    public:

        /// Constructor.
        //  Use [b]objName[/b] to identify the object containg the array, and
        //  [b]attName[/b] to identify the array itself.
        SRefExpr(const std::string &objName, const std::string &attName);

        SRefExpr(const SRefExpr<T> &rhs);
        virtual ~SRefExpr();

        /// Make clone.
        virtual Scalar<T> *clone() const;

        /// Evaluate.
        virtual T evaluate() const;

        /// Invalidate.
        //  Force re-evaluation of the reference.
        virtual void invalidate();

        /// Print expression.
        virtual void print(std::ostream &os, int precedence = 99) const;

    private:

        // Not implemented.
        void operator=(const SRefExpr &);

        // Fill in the reference.
        void fill() const;

        // Make print image.
        const std::string getImage() const;

        // The referred object and attribute names.
        const std::string obj_name;
        const std::string att_name;

        // The object and attribute referred to.
        mutable Object    *itsObject;
        mutable Attribute *itsAttr;
    };


    // Implementation
    // ------------------------------------------------------------------------

    template <class T>
    SRefExpr<T>::SRefExpr
    (const std::string &objName, const std::string &attName):
        obj_name(objName), att_name(attName),
        itsObject(0), itsAttr(0)
    {}


    template <class T>
    SRefExpr<T>::SRefExpr(const SRefExpr<T> &rhs):
        Scalar<T>(rhs),Invalidator(rhs),
        obj_name(rhs.obj_name), att_name(rhs.att_name),
        itsObject(rhs.itsObject), itsAttr(rhs.itsAttr)
    {}


    template <class T>
    SRefExpr<T>::~SRefExpr() {
        if(itsObject) itsObject->unregisterReference(this);
    }


    template <class T>
    Scalar<T> *SRefExpr<T>::clone() const {
        return new SRefExpr<T>(*this);
    }


    template <class T>
    inline T SRefExpr<T>::evaluate() const {
        fill();

        if(AttributeBase *base = &itsAttr->getBase()) {
            if(SValue<T> *value = dynamic_cast<SValue<T>*>(base)) {
                return value->evaluate();
            } else {
                throw ParseError("SRefExpr::evaluate()", "Reference \"" +
                                 getImage() + "\" is not a variable.");
            }
        }

        return T(0);
    }


    template <class T>
    const std::string SRefExpr<T>::getImage() const {
        std::ostringstream os;
        print(os);
        os << std::ends;
        return os.str();
    }


    template <class T>
    void SRefExpr<T>::invalidate() {
        itsObject = 0;
        itsAttr = 0;
    }


    template <class T>
    void SRefExpr<T>::print(std::ostream &os, int) const {
        os << obj_name;
        if(! att_name.empty()) os << "->" << att_name;
        return;
    }


    template <class T>
    void SRefExpr<T>::fill() const {
        if(itsObject == 0) {
            itsObject = OpalData::getInstance()->find(obj_name);
            if(itsObject == 0) {
                if(att_name.empty()) {
                    throw ParseError("SRefExpr::fill()",
                                     "\nThe <variable> \"" + obj_name + "\" is unknown.\n");
                } else {
                    throw ParseError("SRefExpr::fill()",
                                     "Object \"" + obj_name + "\" is unknown.");
                }
            }

            // Register the reference with the object, to allow invalidation
            // when the object is deleted.
            itsObject->registerReference(const_cast<SRefExpr<T>*>(this));

            if(att_name.empty()) {
                itsAttr = itsObject->findAttribute("VALUE");
                if(itsAttr == 0) {
                    throw ParseError("SRefExpr::fill()", "Object \"" + obj_name +
                                     "\" is not a variable, constant or vector.");
                }
            } else {
                itsAttr = itsObject->findAttribute(att_name);
                if(itsAttr == 0) {
                    throw ParseError("SRefExpr::fill()", "Attribute \"" + obj_name +
                                     "->" + att_name + "\" is unknown.");
                }
            }
        }
    }

}

#endif // OPAL_SRefExpr_HH