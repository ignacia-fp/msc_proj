//
//  Block_cluster_tree.h
//  Block cluster tree
//
//  Created by María Ignacia Fierro on 11-12-17.
//  Copyright © 2017 María Ignacia Fierro. All rights reserved.
//

#ifndef block_cluster_tree_elements_h
#define block_cluster_tree_elements_h
#include "basis_functions.h"
#include "complex_math.h"

typedef struct _clusterfactory clusterfactory;                                                                
typedef clusterfactory *pclusterfactory;                                                                      
struct _clusterfactory {                                                                                      
    double **Vp;
    double **Vb;
    int ndofp;
    int nidxp;
    int ndofb;                                                                                                
    int nidxb;
    int d;                                                                                                    
    double *bmin;                                                                                             
    double *bmax;
};                                                                                                            
                                                                                                              
typedef struct _cluster cluster;                                                                              
typedef cluster *pcluster;                                                                                    
struct _cluster {                                                                                             
    int start; /* First element of the associated index set */                                                
    int size; /* Size of the associated index set */                                                          
    int sons; /* Number of sons */
    int d;
    double *bmin;                                                                                             
    double *bmax;
    pMeshes M;
    pcluster *son; /* Array of sons */                                                                        
};                                                                                                            
                                                                                                              
typedef struct _clustertree clustertree;                                                                      
typedef clustertree *pclustertree;                                                                            
struct _clustertree {                                                                                         
    int ndof;                                                                                                 
    int nidx;                                                                                                                                                                                             
    pcluster root;                                                                                            
};

typedef struct _fullmatrix fullmatrix;
typedef fullmatrix* pfullmatrix;
struct _fullmatrix{
    int rows;
    int cols;
    gsl_matrix_complex *e;
};

typedef struct _blockcluster blockcluster;
typedef blockcluster *pblockcluster;
struct _blockcluster {
    pcluster row;
    pcluster col;
    int type;
    pblockcluster *son;
    int block_rows;
    int block_cols;
};

typedef struct _vectorcluster vectorcluster;
typedef vectorcluster* pvectorcluster;
struct _vectorcluster{
    gsl_vector_complex *v;
    pvectorcluster *vc;
    int start;
    int size;
    int sons;
};

typedef struct _hvector hvector;
typedef hvector* phvector;
struct _hvector{
    int *hf;
    pvectorcluster in;
    pvectorcluster out;
};

typedef struct _rkmatrix rkmatrix;
typedef rkmatrix* prkmatrix;
struct _rkmatrix{
    int k; /* maximal rank */
    int kt; /* current rank */
    int rows;
    int cols;
    gsl_matrix_complex *a;
    gsl_matrix_complex *b;
};

typedef struct _supermatrix supermatrix;
typedef supermatrix* psupermatrix;
struct _supermatrix{
    pcluster rows;
    pcluster cols;
    int block_rows;
    int block_cols;
    prkmatrix rk;
    pfullmatrix full;
    psupermatrix *s;
    gsl_vector_complex *in;
    gsl_vector_complex *out;
    int depth;
    int nv;
    int *index;
    int *hfi;
    int *hfo;
    double complex **vecin;                                                                                   
    double complex **vecout;
};

typedef struct _supermataux supermataux;                                                                      
typedef supermataux* psupermataux;                                                                            
struct _supermataux{                                                                                                                                                                               
    psupermatrix *s;
    int *index;
    int size; 
};                                                                                                            
  
#endif /* Block_cluster_tree_elements_h */
