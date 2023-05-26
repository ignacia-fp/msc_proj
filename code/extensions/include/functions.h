#ifndef functions_h
#define functions_h
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include <omp.h>
#include <string.h>
#include "complex_math.h"
//#include <gsl/gsl_math.h>
// LU decomoposition of a general matrix
void dgetrf_(int* M, int *N, double* A, int* lda, int* IPIV, int* INFO);
// generate inverse of a matrix given its LU decomposition
void dgetri_(int* N, double* A, int* lda, int* IPIV, double* WORK, int* lwork, int* INFO);

int check_if_equal( int *v1, int *v2 );

void inverse(int N, double *A);

double norm(double *x);

void vec_prod(int sz, double *a, double *b, double *c, double *res);

void vec_prod_maxwell(int sz,double *a, double *b,double *phi, double *res); 

gsl_vector_complex *vec_prod2(int sz,double *a, double *b,double *c);

void matrix_prod(int s_r_A, int s_c_A, int s_r_B, int s_c_B, double *A, double *B);

double dot_prod(int sz, int idx1, int idx2, double *a, double *b);

double dot_prod2(int sz, double *a, double *b);

void  mat_vec_prod(int s_r_A, int s_c_A, int idx,double *A, double *b, double *res);

void transpose( int row, int col, double *src, double *dst);

void sust_add_vec(double *a, double *b, int opt, int index1, int index2, double *res);

void sust_add_vec_mat(int col_A, int opt, double *A, double *b, int row_A, double *res);

double sum_v(double *v, int sz);

int is_elmt(int e, int *array, int sz);

#endif /* functions_h */

