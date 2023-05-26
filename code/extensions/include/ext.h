#ifndef ext_h
#include "block_cluster_tree.h"
#include "GMRES.h"
#include "assembly.h"
#include <gsl/gsl_spmatrix.h>                                                                               

pMeshes new_space( double *PyVp, double *PyVb, int *PyDp, int *PyDb, int *PyDq, double *PyTp, int *PyEdgep, int *PyEdgeb, int size_of_Vp, int size_of_Vb, int size_of_Dp, int size_of_Db,  int size_of_Dq, int size_of_Tp, int size_of_Edgep, int size_of_Edgeb, char *trialbf, char *testbf );

pMesh new_mesh( double *PyVp, double *PyVb, int *PyDp, int *PyDb, int *PyDq, double *PyTp, double *PyTb, int *PyEdgep, int *PyEdgeb, int *PyEdgep2, int *PyEdgeb2, int size_of_Vp, int size_of_Vb, int size_of_Dp, int size_of_Db,  int size_of_Dq, int size_of_Tp, int size_of_Tb, int size_of_Edgep, int size_of_Edgeb);

psupermatrix hierarchical_mat_slp(pMeshes Mesh, double eta, int ls, double eps, int *order, double kappa, double *time);

psupermatrix hierarchical_mat_hlp(pMeshes Mesh, double eta, int ls, double eps, int *order, double kappa, double *time);

void del_hierarchical_mat(psupermatrix s);

double get_gsl_real( gsl_matrix_complex *M, int row, int col );

double get_gsl_imag( gsl_matrix_complex *M, int row, int col );

void del_gsl_mat(gsl_matrix_complex *M);

void del_potential(double *MR, double *MI);

void del_id(double *M);

void supermat2supermataux(psupermatrix s, psupermataux sa);

//double *gmres_prec_hmat( double *MR, double *MI, double *MMR, double *MMI, psupermatrix P,  double *br, double *bi, double *xr, double *xi,  double tol, int max_iter, int size_of_b, double *time );

//double *gmres_prec( double *MR, double *MI, double *MMR, double *MMI, double *PR, double *Pi, double *br, double *bi, double *xr, double *xi,  double tol, int max_iter, int size_of_b, double *time );

//double *gmres( double *MR, double *MI, double *br, double *bi, double *xr, double *xi,  double tol, int max_iter, int size_of_b, double *time );

void helmholtz_assemble_matrix_sl (double *MR, double *MI, pMeshes Mesh, double kappa, int *order, int opt, double *time);

void helmholtz_assemble_matrix_hl (double *MR, double *MI, pMeshes Mesh, double kappa, int *order, double *time, int opt);

void helmholtz_assemble_matrix_id (double *MR, pMeshes Mesh,  int order, int opt, double *time);

void helmholtz_project_mat(double *MRb, double *MIb, double *MR, double *MI,pMeshes Mesh, double *time);

void maxwell_weakly(double *MR, double *MI, pMeshes Mesh, double kappa, int *order, double *time, int opt);

void maxwell_hypersingular(double *MR, double *MI, pSpace space_trial, pSpace space_test, double kappa, int *order, double *time, int opt);
//void maxwell_hypersingular(double *MR, double *MI, pMeshes Mesh, double kappa, int *order, double *time, int opt);

void maxwell_id(double *MR, pMeshes Mesh, double kappa, int *order, double *time, int opt);

void maxwell_lb(double *MR, pMeshes Mesh, double kappa, int *order, double *time, int opt);

void mv(psupermatrix P, double *br, double *bi, double *t, int size_of_b);

gsl_matrix_complex *mtmt(psupermatrix P, double *MR, double *MI, double *t, int size_of_b);

void mv2(psupermatrix P, double *br, double *bi, double *t, int size_of_b);

void supermat2supermataux2(psupermatrix s, psupermataux sa, double complex *in, double complex *out); 
#endif
