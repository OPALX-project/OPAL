// -*- C++ -*-
/***************************************************************************
 *
 * The IPPL Framework
 *
 ***************************************************************************/

#ifndef BRICK_EXPRESSION_H
#define BRICK_EXPRESSION_H

// include files
#include "Field/AssignTags.h"

// template<unsigned Dim, class LHS, class RHS, class OP>
// class BrickExpression : public BrickExpressionBase
template<unsigned Dim, class LHS, class RHS, class OP>
class BrickExpression
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
