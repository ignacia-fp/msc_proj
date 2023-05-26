#ifndef hmat_arithmetics_h
#define hmat_arithmetics_h
#include "complex_math.h"
#include "block_cluster_tree_elements.h"
#include <sys/time.h>
#include <time.h>
#include <mpi.h>
gsl_complex mult_mat_vec(gsl_vector_complex *a, gsl_vector_complex *b, int *in);

gsl_complex mult_mat_vec2(gsl_vector_complex *a, gsl_vector_complex *b);

void full_mult(const pfullmatrix f, pcluster rows, pcluster cols, gsl_vector_complex *in, gsl_vector_complex *out);

void low_rank_mult(const prkmatrix r, pcluster rows, pcluster cols,gsl_vector_complex *in, gsl_vector_complex *out);

void full_mult2(const pfullmatrix f, gsl_vector_complex *in, gsl_vector_complex *out, int *inv, int *outv);

void low_rank_mult2(const prkmatrix r, gsl_vector_complex *in, gsl_vector_complex *out, int *inv, int *outv);

void set_vec(psupermatrix s, const gsl_vector_complex * x);

void mat_vec(const psupermataux s, const gsl_vector_complex * x, gsl_vector_complex * res, double *t);

void mat_vec2(const psupermataux s, const gsl_vector_complex * x, gsl_vector_complex * res, double *t); 

void mat_mat(psupermataux A, gsl_matrix_complex *B, double *t);
#endif 




