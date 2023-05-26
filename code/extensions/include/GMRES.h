#ifndef GMRES_complex_h
#define GMRES_complex_h
#include <gsl/gsl_spmatrix.h>                                                                                 
#include <gsl/gsl_spblas.h> 
#include "config.h"
#include "hmat_arithmetics.h"
typedef struct
{
    size_t n;        /* size of linear system */
    size_t m;        /* dimension of Krylov subspace K_m */
    gsl_vector_complex *r;   /* residual vector r = b - A*x */
    gsl_matrix_complex *H;   /* Hessenberg matrix n-by-(m+1) */
    gsl_vector_complex *tau; /* householder scalars */
    gsl_vector_complex *y;   /* least squares rhs and solution vector */
    
    double complex *c;       /* Givens rotations */
    double complex *s;
    
    double normr;    /* residual norm ||r|| */
} gmres_state_t_c;

typedef struct
{
    const char *name;
    void * (*alloc) (const size_t n, const size_t m);
    int (*iterate) (const gsl_matrix_complex *A, const gsl_vector_complex *b,const double tol, gsl_vector_complex *x, void *);
    double (*normr)(const void *);
    void (*free) (void *);
} gsl_linalg_itersolve_type_c;

typedef struct
{
    const gsl_linalg_itersolve_type_c * type;
    double normr; /* current residual norm || b - A x || */
    gmres_state_t_c * state;
} gsl_linalg_itersolve_c;


gsl_complex gsl_linalg_complex_householder_transform2 (gsl_vector_complex * v);

int gsl_linalg_complex_householder_hv2 (gsl_complex tau, const gsl_vector_complex * v, gsl_vector_complex *  w);

void gsl_linalg_complex_givens (const gsl_complex a, const gsl_complex b, complex double *c, complex double *s);

void gsl_linalg_complex_givens_gv (gsl_vector_complex * v, const size_t i, const size_t j, double complex a, double complex b);

void gsl_linalg_itersolve_alloc_c(const gsl_linalg_itersolve_type_c *T, gsl_linalg_itersolve_c *w,
                                  const size_t n, const size_t m);

void gsl_linalg_itersolve_free_c(gsl_linalg_itersolve_c *w);

gmres_state_t_c *gmres_alloc_c(const size_t n, const size_t m);

void gmres_free_c(gmres_state_t_c *state);

int gmres_iterate_c(const gsl_matrix_complex *A, const gsl_matrix_complex *P, int prec, const gsl_vector_complex *b, const double tol, gsl_vector_complex *x, gmres_state_t_c *state, double *res);

int gmres_iterate_prec(const gsl_matrix_complex *A, psupermataux sa, psupermatrix P,  gsl_matrix_complex *MM, const gsl_vector_complex *b, const double tol, gsl_vector_complex *x, gmres_state_t_c *state, double *res, double *tm);

double gmres_normr_c(const gmres_state_t_c *state);



#endif
