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

#include <gsl/gsl_integration.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_sf_pow_int.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_sf_gamma.h>
#include "tanhDeriv.h"

namespace tanhderiv {

struct my_f_params {
    double a;
    double s0;
    double lambdaleft;
    double lambdaright;
    double r;
    int n;
};

double my_f (double x, void *p) {
    struct my_f_params *params = (struct my_f_params *)p;
    gsl_complex z = gsl_complex_add(gsl_complex_rect(params->a, 0),
                                    gsl_complex_polar(params->r, x));
    gsl_complex z1 = gsl_complex_div(gsl_complex_add(z,
                                     gsl_complex_rect(params->s0, 0)),
                                     gsl_complex_rect(params->lambdaleft, 0));
    gsl_complex z2 = gsl_complex_div(gsl_complex_sub(z,
                                     gsl_complex_rect(params->s0, 0)),
                                     gsl_complex_rect(params->lambdaright, 0));
    gsl_complex func = gsl_complex_div(gsl_complex_sub(gsl_complex_tanh(z1),
                                       gsl_complex_tanh(z2)),
                                       gsl_complex_rect(2, 0));
    func = gsl_complex_mul(func, gsl_complex_polar(1, -params->n * x));
    return gsl_sf_fact(params->n) * GSL_REAL(func)
           / (2 * M_PI * gsl_sf_pow_int(params->r, params->n));
}

double integrate(const double &a,
                 const double &s0,
                 const double &lambdaleft,
                 const double &lambdaright,
                 const int &n) {
    gsl_function F;
    double radius = gsl_hypot(a - 2, lambdaright * M_PI / 2) - 0.01;
    double radius2 = gsl_hypot(a + 2, lambdaleft * M_PI / 2) - 0.01;
    if (radius > radius2) radius = radius2;
    struct my_f_params params = {a, s0, lambdaleft, lambdaright, radius, n};
    F.function = &my_f;
    F.params = &params;
    gsl_integration_workspace *w = gsl_integration_workspace_alloc(100);
    double error = gsl_sf_pow_int(10, -12);
    double *result = new double;
    double *abserr = new double;
    gsl_set_error_handler_off();
    int status = gsl_integration_qag(&F, 0, 2 * M_PI, 0, error,
                                     100, 6, w, result, abserr);
    gsl_integration_workspace_free(w);
    double finalResult = *result;
    delete result;
    delete abserr;
    if (status) {
        return 0;
    }
    else return finalResult;
}

}
