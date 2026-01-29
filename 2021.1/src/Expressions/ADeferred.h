#ifndef OPAL_ADeferred_HH
#define OPAL_ADeferred_HH

// ------------------------------------------------------------------------
// $RCSfile: ADeferred.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.3.4.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Template class: ADeferred<T>
//
// ------------------------------------------------------------------------
//
// $Date: 2004/11/18 22:52:40 $
// $Author: jsberg $
//
// ------------------------------------------------------------------------

#include "AbstractObjects/Attribute.h"
#include "AbstractObjects/Expressions.h"
#include "Expressions/AList.h"
#include "Expressions/AValue.h"
#include "Utilities/LogicalError.h"
#include "Utilities/OpalException.h"
#include <iosfwd>
#include <vector>
#include <exception>

namespace Expressions {

    // Template class ADeferred
    // ----------------------------------------------------------------------
    /// Object attribute with a ``deferred'' array value.
    //  A deferred expression is always re-evaluated when its value is
    //  required.  This is notably needed to implement random generators.

    template <class T>
    class ADeferred: public AValue<T> {

    public:

        /// Default constructor.
        //  Construct empty array.
        ADeferred();

        /// Constructor.
        //  From array expression.
        explicit ADeferred(PtrToArray<T> expr);

        /// Constructor.
        //  From array constant.
        explicit ADeferred(const std::vector<T> &val);

        /// Constructor.
        //  From array of scalar expressions.
        explicit ADeferred(ArrayOfPtrs<T> expr);

        ADeferred(const ADeferred<T> &);
        virtual ~ADeferred();

        /// Make clone.
        virtual ADeferred<T> *clone() const;

        /// Evaluate.
        //  The expression is always re-evaluated.
        virtual std::vector<T> evaluate();

        /// Get expression flag.
        virtual bool isExpression() const;

        /// Print the attribute value.
        virtual void print(std::ostream &) const;

        /// Set a component of the value.
        void setComponent(int i, const PtrToScalar<T> expr);

    protected:

        /// The generating law for the array expression.
        //  This expression is used first to evaluate all components of
        //  the expression.  Individual values may be overwritten using
        //  [b]itsOverrides[b].
        PtrToArray<T> itsLaw;

        /// Overrides for single components.
        //  After applying [b]itsLaw[/b] to evaluate all components, any
        //  scalar expression defined in this array is used to overwrite
        //  the corresponding component.
        ArrayOfPtrs<T> itsOverrides;

    private:

        // Not implemented.
        void operator=(const ADeferred<T> &);

        // Guard against recursive definition.
        mutable bool in_evaluation;
    };


    // Implementation
    // ----------------------------------------------------------------------

    template <class T>
    ADeferred<T>::ADeferred():
        AValue<T>(),
        itsLaw(),
        itsOverrides(),
        in_evaluation(false)
    {}


    template <class T>
    ADeferred<T>::ADeferred(const ADeferred<T> &rhs):
        AValue<T>(),
        itsLaw(rhs.itsLaw->clone()),
        itsOverrides(rhs.itsOverrides.size()),
        in_evaluation(false) {
        for(typename ArrayOfPtrs<T>::size_type i = 0; i < rhs.itsOverrides.size(); ++i) {
            itsOverrides[i] = rhs.itsOverrides[i]->clone();
        }
    }


    template <class T>
    ADeferred<T>::ADeferred(const std::vector<T> &val):
        AValue<T>(val),
        itsLaw(0),
        itsOverrides(0),
        in_evaluation(false)
    {}


    template <class T>
    ADeferred<T>::ADeferred(PtrToArray<T> expr):
        AValue<T>(),
        itsLaw(expr),
        itsOverrides(),
        in_evaluation(false)
    {}


    template <class T>
    ADeferred<T>::ADeferred(ArrayOfPtrs<T> expr):
        AValue<T>(),
        itsLaw(),
        itsOverrides(expr),
        in_evaluation(false)
    {}


    template <class T>
    ADeferred<T>::~ADeferred()
    {}


    template <class T>
    ADeferred<T> *ADeferred<T>::clone() const {
        return new ADeferred<T>(*this);
    }


    template <class T>
    std::vector<T> ADeferred<T>::evaluate() {
        if(in_evaluation) {
            throw LogicalError("ADeferred::evaluate()",
                               "Recursive expression definitions found.");
        } else {
            in_evaluation = true;

            try {
                if(itsLaw.isValid()) this->value = itsLaw->evaluate();
                in_evaluation = false;
                if(this->value.size() < itsOverrides.size()) {
                    this->value.resize(itsOverrides.size());
                }

                for(typename ArrayOfPtrs<T>::size_type i = 0; i < itsOverrides.size(); ++i) {
                    if(itsOverrides[i].isValid()) {
                        this->value[i] = itsOverrides[i]->evaluate();
                    }
                }
            } catch(OpalException &ex) {
                in_evaluation = false;
                throw LogicalError(ex.where(),
                                   "Evaluating expression \"" +
                                   this->getImage() + "\": " + ex.what());
            } catch(ClassicException &ex) {
                in_evaluation = false;
                throw LogicalError(ex.where(),
                                   "Evaluating expression \"" +
                                   this->getImage() + "\": " + ex.what());
            } catch(std::exception &ex) {
                in_evaluation = false;
                throw LogicalError("ADeferred::evaluate()",
                                   "Standard C++ exception while evaluating \"" +
                                   this->getImage() + "\": " + ex.what());
            } catch(...) {
                in_evaluation = false;
                throw LogicalError("ADeferred::evaluate()",
                                   "Unknown exception while evaluating \"" +
                                   this->getImage() + "\": ");
            }
        }

        return this->value;
    }


    template <class T>
    bool ADeferred<T>::isExpression() const {
        return itsLaw.isValid() || (! itsOverrides.empty());
    }


    template <class T>
    void ADeferred<T>::print(std::ostream &os) const {
        // Print the generating law, if any.
        if(itsLaw.isValid()) itsLaw->print(os, 0);

        // Print the overrides.
        if(! itsOverrides.empty()) {
            os << "; overrides: {";
            typename ArrayOfPtrs<T>::const_iterator i = itsOverrides.begin();

            while(i != itsOverrides.end()) {
                if(i->isValid())(*i)->print(os);
                if(++i == itsOverrides.end()) break;
                os << ',';
            }

            os << '}';
        }

        return;
    }


    template <class T>
    void ADeferred<T>::setComponent
    (int index, const PtrToScalar<T> expr) {
        while(int(itsOverrides.size()) < index) {
            itsOverrides.push_back(PtrToScalar<T>());
        }

        itsOverrides[index-1] = expr;
    }

}

#endif // OPAL_ADeferred_HH
