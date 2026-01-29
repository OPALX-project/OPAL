/*
 *  Copyright (c) 2018, Martin Duy Tat
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

#ifndef TANH_INTEGRAL
#define TANH_INTEGRAL

/** ---------------------------------------------------------------------
  *
  * Integrate performs a contour integral to find the nth derivative of the \n
  * Tanh model fringe function. It uses Cauchy's integral formula to do so. \n
  *
  * ---------------------------------------------------------------------
  *
  * Class category: AbsBeamline \n
  * Author: Martin Duy Tat\n
  *
  * ---------------------------------------------------------------------
  */

namespace tanhderiv {

/** Integrand
 *  \param x -> Independent parameter to be integrated
 *  \param p -> Pointer to struct of parameters
 */
double my_f (double x, void *p);
/** Perform Cauchy's integral to find derivative
 *  \param a -> Point of differentiation
 *  \param s0 -> Centre fringe length
 *  \param lambdaleft -> Left fringe field length
 *  \param lambdaright -> Left fringe field length
 *  \param n -> nth derivative
 */
double integrate(const double &a,
                 const double &s0,
                 const double &lambdaleft,
                 const double &lambdaright,
                 const int &n);

}

#endif
