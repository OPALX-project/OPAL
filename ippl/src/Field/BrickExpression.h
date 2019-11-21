// -*- C++ -*-
/***************************************************************************
 *
 * The IPPL Framework
 * 
 *
 * Visit http://people.web.psi.ch/adelmann/ for more details
 *
 ***************************************************************************/

#ifndef BRICK_EXPRESSION_H
#define BRICK_EXPRESSION_H

// include files
#include "Utility/Pooled.h"
#include "Utility/RefCounted.h"
#include "Field/AssignTags.h"

//////////////////////////////////////////////////////////////////////
/*
class BrickExpressionBase : public RefCounted, public Pooled
{
public:
  BrickExpressionBase() {}
  virtual ~BrickExpressionBase() {}
  virtual void apply() restrict =0;
};
*/
//////////////////////////////////////////////////////////////////////

// template<unsigned Dim, class LHS, class RHS, class OP>
// class BrickExpression : public BrickExpressionBase
template<unsigned Dim, class LHS, class RHS, class OP>
class BrickExpression :
  public RefCounted, public Pooled< BrickExpression<Dim,LHS,RHS,OP> >
{
public: 
  BrickExpression(const LHS& l, const RHS& r)
    : Lhs(l), Rhs(r)
      {
      }
  BrickExpression(const LHS& l, const RHS& r, const OP& o)
    : Lhs(l), Rhs(r), Op(o)
      {
      }

  virtual void apply();

private:
  LHS Lhs;
  RHS Rhs;
  OP  Op;
};

//////////////////////////////////////////////////////////////////////

#include "Field/BrickExpression.hpp"

#endif // BRICK_EXPRESSION_H

/***************************************************************************
 * $RCSfile: BrickExpression.h,v $   $Author: adelmann $
 * $Revision: 1.1.1.1 $   $Date: 2003/01/23 07:40:26 $
 * IPPL_VERSION_ID: $Id: BrickExpression.h,v 1.1.1.1 2003/01/23 07:40:26 adelmann Exp $ 
 ***************************************************************************/

