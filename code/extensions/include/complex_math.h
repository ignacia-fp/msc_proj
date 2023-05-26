#ifndef complex_math_h
#define complex_math_h

#include "config.h"
#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_complex.h>
#include <complex.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>
#include <string.h>


gsl_complex gsl_complex_div_real(gsl_complex a, double b);

gsl_complex gsl_complex_conjugate(gsl_complex a);

gsl_complex gsl_complex_mul(gsl_complex a, gsl_complex b);

gsl_complex gsl_complex_add(gsl_complex a, gsl_complex b);

gsl_complex gsl_complex_sub(gsl_complex a, gsl_complex b);

gsl_complex gsl_complex_negative(gsl_complex a);

gsl_complex gsl_complex_sub_real(gsl_complex a, double x);

gsl_complex gsl_complex_rect (double x, double y);

gsl_complex gsl_complex_inverse (gsl_complex x);

double gsl_complex_abs2(gsl_complex a);

double gsl_complex_abs(gsl_complex a);

int gsl_complex_eq(gsl_complex a, gsl_complex b);

size_t  gsl_vector_complex_max_index(const gsl_vector_complex *v);
#endif
