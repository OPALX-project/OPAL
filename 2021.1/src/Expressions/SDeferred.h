#ifndef OPAL_SDeferred_HH
#define OPAL_SDeferred_HH

// ------------------------------------------------------------------------
// $RCSfile: SDeferred.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.2.4.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Template class: SDeferred<T>
//
// ------------------------------------------------------------------------
//
// $Date: 2004/11/18 22:52:40 $
// $Author: jsberg $
//
// ------------------------------------------------------------------------

#include "AbstractObjects/Attribute.h"
#include "AbstractObjects/Expressions.h"
#include "Expressions/SValue.h"
#include "Utilities/LogicalError.h"
#include "Utilities/OpalException.h"
#include <iosfwd>
#include <vector>
#include <exception>

namespace Expressions {

    // Template class SDeferred
    // ----------------------------------------------------------------------
    /// Object attribute with a ``deferred'' scalar value.
    //  An deferred expression is always re-evaluated when its value is
    //  required.  This is notably needed to implement random generators.

    template <class T>
    class SDeferred: public SValue<T> {

    public:

        /// Constructor.
        //  From scalar expression.
        explicit SDeferred(PtrToScalar<T> expr);

        SDeferred(const SDeferred<T> &);
        virtual ~SDeferred();

        /// Make clone.
        virtual SDeferred<T> *clone() const;

        /// Evaluate.
        //  The expression is always re-evaluated.
        virtual T evaluate();

        /// Return expression flag.
        //  Always true.
        virtual bool isExpression() const;

        /// Print the attribute value.
        virtual void print(std::ostream &) const;

    protected:

        /// Pointer to expression.
        PtrToScalar<T> expr_ptr;

    private:

        // Not implemented.
        SDeferred();
        void operator=(const SDeferred<T> &);

        // Guard against recursive definition.
        mutable bool in_evaluation;
    };


    // Implementation
    // ------------------------------------------------------------------------

    template <class T>
    SDeferred<T>::SDeferred(const SDeferred<T> &rhs):
        SValue<T>(), expr_ptr(rhs.expr_ptr->clone()), in_evaluation(false)
    {}


    template <class T>
    SDeferred<T>::SDeferred(PtrToScalar<T> expr):
        SValue<T>(), expr_ptr(expr), in_evaluation(false)
    {}


    template <class T>
    SDeferred<T>::~SDeferred()
    {}


    template <class T>
    SDeferred<T> *SDeferred<T>::clone() const {
        return new SDeferred<T>(*this);
    }


    template <class T>
    T SDeferred<T>::evaluate() {
        if(in_evaluation) {
            throw LogicalError("SDeferred::evaluate()",
                               "Recursive expression definitions found.");
        } else {
            in_evaluation = true;
            try {
                this->value = expr_ptr->evaluate();
                in_evaluation = false;
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
                throw LogicalError("SDeferred::evaluate()",
                                   "Standard C++ exception while evaluating \"" +
                                   this->getImage() + "\": " + ex.what());
            } catch(...) {
                in_evaluation = false;
                throw LogicalError("SDeferred::evaluate()",
                                   "Unknown exception while evaluating \"" +
                                   this->getImage() + "\": ");
            }
        }

        return this->value;
    }


    template <class T>
    bool SDeferred<T>::isExpression() const {
        return true;
    }


    template <class T>
    void SDeferred<T>::print(std::ostream &stream) const {
        expr_ptr->print(stream, 0);
        return;
    }

}

#endif // OPAL_SDeferred_HH
