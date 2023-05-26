//
//  Block_cluster_tree.h
//  Block cluster tree
//
//  Created by María Ignacia Fierro on 11-12-17.
//  Copyright © 2017 María Ignacia Fierro. All rights reserved.
//

#ifndef block_cluster_tree_h
#define block_cluster_tree_h
#include "ACA.h"
#include "hmat_arithmetics.h"
#define HLIB_BLOCK_MINADM 3
#define HLIB_BLOCK_WEAKADM 1
#define HLIB_BLOCK_STRONGADM 7
#define MIN(a,b) (((a)<(b))?(a):(b))
#define MAX(a,b) (((a)>(b))?(a):(b))

double distance_cluster(pcluster row, pcluster col);

double diameter_cluster(pcluster c);

pcluster new_cluster(int start, int size, int d, int sons);

pclustertree new_clustertree(int n, int nidx);

pclusterfactory new_clusterfactory( int ndofp, int nidxp, int ndofb, int nidxb, int d);

pblockcluster new_blockcluster(pcluster row, pcluster col, int block_rows, int block_cols, unsigned type);

pfullmatrix new_fullmatrix_P0D(pcluster rows, pcluster cols, int *index, double kappa, int *order, double eps);

pfullmatrix new_fullmatrix_P1D(gsl_matrix_complex *M, pcluster rows, pcluster cols, int *index, double kappa, int *order, double eps);

prkmatrix new_rkmatrix_P0D(pcluster rows, pcluster cols, int *index, double kappa, int *order, double eps);

prkmatrix new_rkmatrix_P1D(gsl_matrix_complex *M, pcluster rows, pcluster cols, int *index, double kappa, int *order, double eps);

prkmatrix new_rkmatrix_P1D_hlp(gsl_matrix_complex *M, pcluster rows, pcluster cols, int *index, double kappa, int *order, double eps);

pfullmatrix new_fullmatrix_P1D_hlp(gsl_matrix_complex *M, pcluster rows, pcluster cols, int *index, double kappa, int *order, double eps);

psupermatrix new_supermatrix(int block_rows, int block_cols, pcluster rows, pcluster cols,  psupermatrix *s, prkmatrix rk, pfullmatrix full);

void del_pcluster(pcluster c);

void del_clustertree(pclustertree c);

void del_clusterfactory(pclusterfactory c);

void del_rkmatrix(prkmatrix r);

void del_fullmatrix(pfullmatrix f);

void del_supermatrix(psupermatrix s);

pcluster split_boundingbox_P0D(pclusterfactory factory, pMeshes Mesh, int d, int *indexp, int leafsize, int startp, int sizep);

pcluster split_boundingbox_P1D(pclusterfactory factory, pMeshes Mesh, int d, int *indexp, int *indexb, int *notriangle, int leafsize, int startp, int sizep, int startb, int sizeb, int lnt);

pclustertree create_subclustertree(pclusterfactory factory, pMeshes Mesh, int d, const int *indexp, int ndp, int leafsize);

pclustertree create_subclustertree_PD1(pclusterfactory factory, pMeshes Mesh, int d, const int *indexp,int ndp, const int *indexb,int ndb,int leafsize);
//void find_boundingbox(pcluster tau, pclusterfactory factory, int *index, int size);

pblockcluster build_blockcluster(pcluster row, pcluster col, double eta);

psupermatrix build_supermatrix_from_blockcluster_P0D(pblockcluster bc, int *index, double kappa, int *order, double eps);

psupermatrix build_supermatrix_from_blockcluster_P1D(gsl_matrix_complex *M, pblockcluster bc, int *index, double kappa, int *order, double eps, int opt);

psupermatrix create_hmat(pMeshes Mesh, double eta, int ls, double eps, int *order, double kappa);

psupermatrix create_hmat_hlp(pMeshes Mesh, double eta, int ls, double eps, int *order, double kappa);
#endif /* Block_cluster_tree_h */
