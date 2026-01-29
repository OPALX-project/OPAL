/* 
 *  Copyright (c) 2015, Chris Rogers
 *  All rights reserved.
 *  Redistribution and use in source and binary forms, with or without 
 *  modification, are permitted provided that the following conditions are met: 
 *  1. Redistributions of source code must retain the above copyright notice,
 *     this list of conditions and the following disclaimer. 
 *  2. Redistributions in binary form must reproduce the above copyright notice, 
 *     this list of conditions and the following disclaimer in the documentation 
 *     and/or other materials provided with the distribution.
 *  3. Neither the name of STFC nor the names of its contributors may be used to 
 *     endorse or promote products derived from this software without specific 
 *     prior written permission.
 *
 *  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" 
 *  AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE 
 *  IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE 
 *  ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE 
 *  LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR 
 *  CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF 
 *  SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS 
 *  INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN 
 *  CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) 
 *  ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE 
 *  POSSIBILITY OF SUCH DAMAGE.
 */

#ifndef PolynomialCoefficient_h
#define PolynomialCoefficient_h

namespace interpolation {

/** \class PolynomialCoefficient
 * \brief PolynomialCoefficient represents a coefficient in a multi-dimensional
 *  polynomial
 *
 *  PolynomialCoefficient has three member data
 *  \param inVarByVec_m is the x power to which the coefficient pertains,
 *  indexed by vector e.g. \f$x_1^4 x_2^3 =\f$ {1,1,1,1,2,2,2}
 *  \param outVar_m is the output y variable to which the coefficient pertains
 *  \param coefficient is the value of the coefficient
 */
class PolynomialCoefficient
{
  public:
    /** Construct the coefficient
     *  \param inVariablesByVector x power indexed like e.g.
     *  \param outVariable y index that the coefficient pertains to
     *  \param coefficient value of the coefficient
     */
    PolynomialCoefficient(std::vector<int> inVariablesByVector,
                          int outVariable,
                          double coefficient) 
        : _inVarByVec(inVariablesByVector), _outVar(outVariable),
          _coefficient(coefficient) {
    }


    /** Copy constructor */
    PolynomialCoefficient(const PolynomialCoefficient& pc)
        : _inVarByVec(pc._inVarByVec), _outVar(pc._outVar),
          _coefficient(pc._coefficient) {
    }

    /** Return the vector of input variables, indexed by vector */
    std::vector<int> InVariables() const {return _inVarByVec;}

    /** Return the output variable */
    int              OutVariable() const {return _outVar;}

    /** Return the coefficient value */
    double           Coefficient() const {return _coefficient;}

    /** Set the vector of input variables, indexed by vector */
    std::vector<int> InVariables(std::vector<int> inVar ) {_inVarByVec  = inVar;  return _inVarByVec;}

    /** Set the output variable */
    int              OutVariable(int              outVar) {_outVar      = outVar; return _outVar;}
    /** Set the coefficient */
    double           Coefficient(double           coeff ) {_coefficient = coeff;  return _coefficient;}

    /** Transform coefficient from subspace space_in to subspace space_out, both
     *  subspaces of some larger space
     *
     *  \param spaceIn  Describes the input subspace
     *  \param spaceOut  Describes the output subspace
     *  Throws an exception if the coefficient is not in the input or output
     *  subspace
     *  So for coeff({1,2},0,1.1), coeff.space_transform({0,2,3,5}, {4,7,1,2,3,0})
     *  would return coeff({3,4},5,1.1)
     */
    void             SpaceTransform(std::vector<int> spaceIn,
                                    std::vector<int> spaceOut);
    /** Return true if var is in inVariables */
    bool             HasInVariable(int var) const {for(unsigned int i=0; i<_inVarByVec.size(); i++) if(_inVarByVec[i] == var) return true; return false;}

  private:
    std::vector<int> _inVarByVec;
    int              _outVar;
    double           _coefficient;
};

std::ostream& operator << (std::ostream&, const PolynomialCoefficient&);

}

#endif
