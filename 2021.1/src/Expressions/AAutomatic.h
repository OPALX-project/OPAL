#ifndef OPAL_AAutomatic_HH
#define OPAL_AAutomatic_HH

// ------------------------------------------------------------------------
// $RCSfile: AAutomatic.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1.4.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Template class: AAutomatic<T>
//
// ------------------------------------------------------------------------
//
// $Date: 2004/11/18 22:52:40 $
// $Author: jsberg $
//
// ------------------------------------------------------------------------

#include "Expressions/ADeferred.h"
#include "AbstractObjects/Expressions.h"
#include "AbstractObjects/OpalData.h"


namespace Expressions {

    // Template class AAutomatic
    // ----------------------------------------------------------------------
    /// Object attribute with an ``automatic'' array value.
    //  An automatic expression is marked as unknown and registered in a list
    //  when it is created.  When its value is required, it is evaluated,
    //  marked as known, and the new value is cached.  Whenever a new
    //  definition is read, or an existing one is changed, the expressions
    //  are all marked as unknown. This forces a new evaluation when an
    //  exoression is used the next time.

    template <class T>
    class AAutomatic: public ADeferred<T> {

    public:

        /// Constructor.
        //  From array of expressions.
        explicit AAutomatic(PtrToArray<T> expr);

        /// Constructor.
        //  From array expression.
        explicit AAutomatic(ArrayOfPtrs<T> expr);

        AAutomatic(const AAutomatic &);
        virtual ~AAutomatic();

        /// Make a clone.
        virtual AAutomatic<T> *clone() const;

        /// Evaluate.
        //  The resulting value is cached.
        virtual std::vector<T> evaluate();

        /// Invalidate.
        //  Mark as unknown to force re-evaluation of expression.
        virtual void invalidate();

    private:

        // Not implemented.
        AAutomatic();
        void operator=(const AAutomatic &);

        // Is the expression known ?
        mutable bool is_known;
    };


    // Implementation
    // ----------------------------------------------------------------------

    template <class T>
    AAutomatic<T>::AAutomatic(const AAutomatic<T> &rhs):
        ADeferred<T>(rhs),
        is_known(false) {
        OpalData::getInstance()->registerExpression(this);
    }


    template <class T>
    AAutomatic<T>::AAutomatic(PtrToArray<T> expr):
        ADeferred<T>(expr),
        is_known(false) {
        OpalData::getInstance()->registerExpression(this);
    }


    template <class T>
    AAutomatic<T>::AAutomatic(ArrayOfPtrs<T> expr):
        ADeferred<T>(expr),
        is_known(false) {
        OpalData::getInstance()->registerExpression(this);
    }


    template <class T>
    AAutomatic<T>::~AAutomatic() {
        // Unlink expression from expression list.
        OpalData::getInstance()->unregisterExpression(this);
    }


    template <class T>
    AAutomatic<T> *AAutomatic<T>::clone() const {
        return new AAutomatic<T>(*this);
    }


    template <class T>
    std::vector<T> AAutomatic<T>::evaluate() {
        if(! is_known) {
            ADeferred<T>::evaluate();
            is_known = true;
        }

        return this->value;
    }


    template <class T>
    void AAutomatic<T>::invalidate() {
        is_known = false;
    }

}

#endif // OPAL_AAutomatic_HH

