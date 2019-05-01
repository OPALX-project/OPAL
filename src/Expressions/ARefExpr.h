#ifndef OPAL_ARefExpr_HH
#define OPAL_ARefExpr_HH

// ------------------------------------------------------------------------
// $RCSfile: ARefExpr.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.3.4.2 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class ARefExpr:
//
// ------------------------------------------------------------------------
//
// $Date: 2004/11/18 22:52:40 $
// $Author: jsberg $
//
// ------------------------------------------------------------------------

#include "AbstractObjects/Expressions.h"
#include "AbstractObjects/Invalidator.h"
#include "Expressions/AValue.h"
#include "Expressions/SValue.h"
#include "Utilities/OpalException.h"
#include "AbstractObjects/OpalData.h"
#include <iosfwd>
#include <sstream>
#include <vector>


namespace Expressions {

    // Class ARefExpr
    // ------------------------------------------------------------------------
    /// An expression defined as a reference to an array.
    //  The components of the array may be real, logical or string.
    //  Reference expressions are (re)evaluated as needed as part of the
    //  containing expression.

    template <class T>
    class ARefExpr: public OArray<T>, public Invalidator {

    public:

        /// Constructor.
        //  Use [b]objName[/b] to identify the object containg the array, and
        //  [b]attName[/b] to identify the array itself.
        ARefExpr(const std::string &objName, const std::string &attName);

        ARefExpr(const ARefExpr<T> &rhs);
        virtual ~ARefExpr();

        /// Make clone.
        virtual OArray<T> *clone() const;

        /// Evaluate the reference and return the value referred to.
        virtual std::vector<T> evaluate() const;

        /// Print expression.
        virtual void print(std::ostream &os, int precedence = 99) const;

    private:

        // Not implemented.
        void operator=(const ARefExpr &);

        // Fill in the reference.
        void fill() const;

        // Make print image.
        const std::string getImage() const;

        // The referred object and attribute.
        const std::string obj_name;
        const std::string att_name;

        // The object and attribute referred to.
        mutable Object    *itsObject;
        mutable Attribute *itsAttr;
    };


    // Implementation
    // ------------------------------------------------------------------------

    template <class T>
    ARefExpr<T>::ARefExpr
    (const std::string &objName, const std::string &attName):
        obj_name(objName), att_name(attName),
        itsObject(0), itsAttr(0)
    {}


    template <class T>
    ARefExpr<T>::ARefExpr(const ARefExpr<T> &rhs):
        OArray<T>(rhs),Invalidator(rhs),
        obj_name(rhs.obj_name), att_name(rhs.att_name),
        itsObject(rhs.itsObject), itsAttr(rhs.itsAttr)
    {}


    template <class T>
    ARefExpr<T>::~ARefExpr() {
        if(itsObject) itsObject->unregisterReference(this);
    }


    template <class T>
    OArray<T> *ARefExpr<T>::clone() const {
        return new ARefExpr<T>(*this);
    }


    template <class T>
    inline std::vector<T> ARefExpr<T>::evaluate() const {
        fill();

        if(AttributeBase *base = &itsAttr->getBase()) {
            if(AValue<T> *value = dynamic_cast<AValue<T>*>(base)) {
                return value->evaluate();
            } else if(SValue<T> *value = dynamic_cast<SValue<T>*>(base)) {
                return std::vector<T>(1, value->evaluate());
            } else {
                throw OpalException("ARefExpr::evaluate()", "Reference \"" +
                                    getImage() + "\" is not an array.");
            }
        } else {
            return std::vector<T>();
        }
    }


    template <class T>
    const std::string ARefExpr<T>::getImage() const {
        std::ostringstream os;
        print(os);
        os << std::ends;
        return os.str();
    }


    template <class T>
    void ARefExpr<T>::print(std::ostream &os, int) const {
        os << obj_name;
        if(! att_name.empty()) os << "->" << att_name;
    }


    template <class T>
    void ARefExpr<T>::fill() const {
        if(itsObject == 0) {
            itsObject = OpalData::getInstance()->find(obj_name);
            if(itsObject == 0) {
                throw OpalException("ARefExpr::fill()",
                                    "Object \"" + obj_name + "\" is unknown.");
            }

            // Register the reference with the object, to allow invalidation
            // when the object is deleted.
            itsObject->registerReference(const_cast<ARefExpr<T>*>(this));

            if(att_name.empty()) {
                itsAttr = itsObject->findAttribute("VALUE");
                if(itsAttr == 0) {
                    throw OpalException("ARefExpr::fill()", "Object \"" + obj_name +
                                        "\" is not a variable, constant or vector.");
                }
            } else {
                itsAttr = itsObject->findAttribute(att_name);
                if(itsAttr == 0) {
                    throw OpalException("ARefExpr::fill()", "Attribute \"" + obj_name +
                                        "->" + att_name + "\" is unknown.");
                }
            }
        }
    }

}

#endif // OPAL_ARefExpr_HH
