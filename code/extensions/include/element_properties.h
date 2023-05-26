#ifndef element_properties_h
#define element_properties_h

#include "functions.h"
#define PI 3.14159265358979323846

int compare_vert(double *v1, double *v2);

double len_edge(double *v1, double *v2);

double deriv_t( double x1, double x2, double x3, int opt );

double deriv_q( double x1, double x2, double x3,double x4,double eta, double chi, int opt );

double area( double *M);

void evaluate_kernel(double *x, double *y, double k, int npoints_x, int npoints_y, double *di, double *dj, double *res);

void evaluate_kernel_maxwell( double *x, double *y, double kappa, int npoints_x, int npoints_y,double *di, double *dj, double *res );

void evaluate_kernel_maxwell_T( double *x, double *y, double kappa, int npoints_x, int npoints_y,double *di, double *dj, double *di2, double *dj2, double *cres1, double *cres2 );

gsl_complex HP_kern( double kappa, double *weights_x, double *weights_y, double *phii, double *phij, double *points_x, double *points_y, double *IE_x, double *IE_y, gsl_vector *n_x, gsl_vector *n_y, gsl_vector *curl_x, gsl_vector *curl_y, int sz_x, int sz_y );

gsl_complex HP_kern2( double kappa, double *weights_x, double *weights_y, double *phii, double *phij, double *points_x, double *points_y, double *IE_x, double *IE_y, gsl_vector *n_x, gsl_vector *n_y, gsl_matrix *curl_x, gsl_matrix *curl_y, int sz_x, int sz_y );

gsl_complex HP_kern3( double kappa, double *weights_x, double *weights_y, double *phii, double *phij, double *points_x, double *points_y, double *IE_x, double *IE_y, gsl_vector *n_x, gsl_vector *n_y, gsl_vector *curl_x, gsl_vector *curl_y, int sz_x, int sz_y );

gsl_complex HP_kern4( double kappa, double *weights_x, double *weights_y, double *phii, double *phij, double *points_x, double *points_y, double *IE_x, double *IE_y, gsl_vector *n_x, gsl_vector *n_y, gsl_matrix *curl_x, gsl_vector *curl_y, int sz_x, int sz_y );

void cross_product(const gsl_vector *u, const gsl_vector *v, gsl_vector *product);

void cross_product2( double *u, gsl_vector *v, int sz );

gsl_vector *get_normal(gsl_vector *p1, gsl_vector *p2, gsl_vector *p3);

gsl_vector *normal(double *points, double *v);

gsl_matrix *invert_a_matrix(gsl_matrix *matrix);

void divergence(double a,double len, double *v,char *tp, int opt, double *div, int sz);

gsl_vector *curl_bf(double *v, int dof, gsl_vector *n);

gsl_vector *curl_bf2(double *v, int dof, gsl_vector *n, double eta, double chi);

gsl_matrix *curl_bfq(double *v, double *eta, double *chi, int sz, int dof, gsl_vector *n);

gsl_matrix *curl_bfq2(double *v, double *eta, double *chi, int sz, int dof, gsl_vector *n);

void get_ie_t( int sz_eta, double *M, double *IE );

void get_ie_q( int sz_eta,int sz_chi,double *M,double *eta, double *chi, double *IE );

void points_t( double *M, int sz_eta, double *eta, double *chi, double *points_aux);

void points_q( double *M, int sz_eta, int sz_chi, double *eta, double *chi, double *points_aux);

double nodal_function_t( double x1, double x2, double x3, double eta, double chi );

double nodal_function_q( double x1, double x2, double x3, double x4, double eta, double chi );

void quad_rules_PRIMAL( double *eta, double *chi, int *sz, int order );

void quad_rules_DUAL( double *eta, double *chi, int *sz, int order );

void weights_DUAL(double *weights, int order, int sz);

void weights_PRIMAL( double *weights, int order, int sz );

void weights_int( double *weights, int order );

void ev_basis_function( int sz, int dof, double *phi, double *eta, double *chi, char *tp );

void ev_basis_function_maxwell(double *points, double *phi, double *p0, double A, double len, gsl_vector *n, int sz, int opt, char *tp);

void evaluate_sing_kern(double *y,  double *M, double *V, int order, double *res, int dof, double kappa, char *tp);

void evaluate_sing_kern_maxwell_weakly(double *y,  double *M, double *p0, double A, double len, gsl_vector *n, int opt, int order, double *res, double kappa);

void evaluate_sing_kern_maxwell_hypersingular(double *y,  double *M, double *V, double a, double len, int order, int opt, double *res, double kappa, char *tp);

void singular(double *P0, double *V, double *weights, double *res, double *IE, int *sz, int order, int dof, double kappa, char *tp, double phi);

void singular_maxwell_weakly(double *P0, double *V, double *p0, double *weights, double A_t, double len, gsl_vector *n, double *res, int *sz,  int order, int opt, double kappa,  double *phi);

void singular_maxwell_hypersingular(double *P0, double *V, double *weights, double A_t, double len, double *res, double *IE, int *sz,  int order, int opt, double kappa, char *tp, double phi);

void singular_maxwell_weakly2(double kappa, int *order, double *v, double *p01, double *p02, double *res, gsl_vector *n, double A_t, double len1, double len2, int opt1, int opt2);

void singular_maxwell_weakly_ce(double kappa, int *order, double *v1, double *v2,  double *p01, double *p02, double *res, gsl_vector *n, double len1, double len2, double A_t1, double A_t2, int *edges1, int *edges2, int opt1, int opt2);

void singular_maxwell_hypersingular2(double kappa, int *order, double *v, double *p01, double *p02, double *res, gsl_vector *n, double A_t, double len1, double len2, int opt1, int opt2);

void singular_maxwell_hypersingular_ce(double kappa, int *order, double *v1, double *v2,  double *p01, double *p02, double *res, gsl_vector *n, double len1, double len2, double A_t1, double A_t2, int *edges1, int *edges2, int opt1, int opt2);

void singular_maxwell_weakly_cv(double kappa, int *order, double *v1, double *v2,  double *p01, double *p02, double *res, gsl_vector *n, double len1, double len2, double A_t1, double A_t2, int *dofs1, int *dofs2, int opt1, int opt2);

void singular_maxwell_hypersingular_cv(double kappa, int *order, double *v1, double *v2,  double *p01, double *p02, double *res, gsl_vector *n, double len1, double len2, double A_t1, double A_t2, int *dofs1, int *dofs2, int opt1, int opt2);

void integration_elements_PRIMAL( int sz_eta, int sz_chi, double *eta, double *chi, double *corners, double *IE, double *P0 );

void integration_elements_DUAL( int sz_eta, int sz_chi, double *eta, double *chi, double *corners, double *IE, double *P0 );

void integration_elements_DUAL_PRIMAL( int sz_eta, int sz_chi, double *corners_t, double *corners_q, double *eta, double *chi,double *eta_p, double *chi_p, double *IE, double *P0 );

void integration_elements_DUAL(int sz_eta, int sz_chi, double *eta, double *chi, double *corners,double *IE,double *P0);

void change_of_basis(double *corners_t, double *eta_p, double *chi_p, int id, double *P0 );
#endif


