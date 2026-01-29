#ifndef OPAL_SAutomatic_HH
#define OPAL_SAutomatic_HH

// ------------------------------------------------------------------------
// $RCSfile: SAutomatic.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1.4.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Template class: SAutomatic<T>
//
// ------------------------------------------------------------------------
//
// $Date: 2004/11/18 22:52:40 $
// $Author: jsberg $
//
// ------------------------------------------------------------------------

#include "Expressions/SDeferred.h"
#include "AbstractObjects/Expressions.h"
#include "AbstractObjects/OpalData.h"


namespace Expressions {

    // Template class SAutomatic
    // ----------------------------------------------------------------------
    /// Object attribute with an ``automatic'' scalar value.
    //  An automatic expression is marked as unknown and registered in a list
    //  when it is created.  When its value is required, it is evaluated, is
    //  marked as known, and the new value is cached.  Whenever a new
    //  definition is read, or an existing one is changed, all expressions
    //  are again marked as unknown.  This forces a new evaluation when an
    //  expression is used the next time.

    template <class T> class SAutomatic: public SDeferred<T> {

    public:

        /// Constructor.
        //  From scalar expression.
        explicit SAutomatic(PtrToScalar<T> expr);

        SAutomatic(const SAutomatic<T> &);
        virtual ~SAutomatic();

        /// Make clone.
        virtual SAutomatic<T> *clone() const;

        /// Evaluate.
        //  The resulting value is cached.
        virtual T evaluate();

        /// Invalidate.
        //  Mark as unknown to force re-evaluation of expression.
        virtual void invalidate();

    private:

        // Not implemented.
        SAutomatic();
        void operator=(const SAutomatic<T> &);

        // Is the expression known ?
        mutable bool is_known;
    };


    // Implementation
    // ------------------------------------------------------------------------

    template <class T>
    SAutomatic<T>::SAutomatic(const SAutomatic<T> &rhs):
        SDeferred<T>(rhs), is_known(false) {
        OpalData::getInstance()->registerExpression(this);
    }


    template <class T>
    SAutomatic<T>::SAutomatic(PtrToScalar<T> expr):
        SDeferred<T>(expr), is_known(false) {
        OpalData::getInstance()->registerExpression(this);
    }


    template <class T>
    SAutomatic<T>::~SAutomatic() {
        // Unlink expression from expression list.
        OpalData::getInstance()->unregisterExpression(this);
    }


    template <class T>
    SAutomatic<T> *SAutomatic<T>::clone() const {
        return new SAutomatic<T>(*this);
    }


    template <class T>
    T SAutomatic<T>::evaluate() {
        if(! is_known) {
            SDeferred<T>::evaluate();
            is_known = true;
        }

        return this->value;
    }


    template <class T>
    void SAutomatic<T>::invalidate() {
        is_known = false;
    }

}

#endif // OPAL_SAutomatic_HH
