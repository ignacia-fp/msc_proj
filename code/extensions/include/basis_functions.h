#ifndef basis_functions_h
#define basis_functions_h

#include "element_properties.h"
                                                                                                              
typedef struct _Triangle Triangle;
typedef Triangle *pTriangle;
struct _Triangle {
    double v[9];
    double   a;
    int dofs[3];
    int edges[3];
    gsl_vector *n; 
    gsl_vector *c;
};

typedef struct _Trianglep Trianglep;                                                                            
typedef Trianglep *pTrianglep;                                                                                  
struct _Trianglep {
    pTriangle *T;
    double v[9];                                                                                              
    double   a; 
    int *num;
    int dofs[3];  
    gsl_vector *n;
    gsl_vector *c;
};

typedef struct _Quadrilateral Quadrilateral;
typedef Quadrilateral *pQuadrilateral;
struct _Quadrilateral {
    double v[12];
    double     a;
    int  dofs[4];
    int dofsp[3];
    gsl_vector *n;
    gsl_matrix *c[3];
};

typedef struct _sP0P sP0P;
typedef sP0P *psP0P;
struct _sP0P {
    pTriangle T;
    int *ntb;
};

typedef struct _sP1P sP1P;
typedef sP1P *psP1P;
struct _sP1P {
    pTriangle *T;
    int len;
    int *num;
    int dv;
};

typedef struct _sPD0 sPD0;
typedef sPD0 *psPD0;
struct _sPD0 {
    pQuadrilateral *Q;
    int *num;
    int len;
    int dv;
};

typedef struct _sPD1 sPD1;
typedef sPD1 *psPD1;
struct _sPD1 {
    pTrianglep *T;
    int *num;
    int len;
};

typedef struct _sPD1b2 sPD1b2;                                                                                  
typedef sPD1b2 *psPD1b2;                                                                                        
struct _sPD1b2 {                                                                                               
    pQuadrilateral *Q;                                                                                        
    pTrianglep *T;  
    int len;                                                                                                  
};

typedef struct _sPD1b sPD1b;                                                                                  
typedef sPD1b *psPD1b;                                                                                        
struct _sPD1b {                                                                                               
    int idx[7];                                                                                               
    int idx2[7]; 
    psPD1b2 *sbf;                                                                                             
};

typedef struct _sRWG sRWG;                                                                                    
typedef sRWG *psRWG;                                                                                          
struct _sRWG {                                                                                                
    pTriangle T1;
    pTriangle T2;
    int num[2];
    double len;
    double p1[3];
    double p2[3]; 
    int edge;                                                                                                 
};

typedef struct _sBC sBC;                                                                                    
typedef sBC *psBC;                                                                                          
struct _sBC {  
    int len1;
    int len2;
    int len3;
    int len4;
    int len5;
    int edge;
    psRWG *sbf1;  
    psRWG *sbf2;
    psRWG *sbf3;
    psRWG *sbf4;
    psRWG *sbf5; 
    gsl_vector_complex *coef1; 
    gsl_vector_complex *coef2;
    gsl_vector_complex *coef3;                                                                                            
    gsl_vector_complex *coef4; 
    gsl_vector_complex *coef5;                                                                                                                                                                                         
}; 

/*typedef struct _sBC sBC;                                                                                    
typedef sRWG *psBC;                                                                                          
struct _sRWG {                                                                                                                        
    psRWG *sbf1;
    psRWG *sbf2;
    double *coef1;
    double *coef2;
    int len1;
    int len2; 
    int edge;                                                                                                 
};*/

typedef struct _P0Pmesh P0Pmesh;
typedef P0Pmesh *pP0Pmesh;
struct _P0Pmesh {
    psP0P *sbf;
    int nt;
    char *type;
    int type2;
};

typedef struct _P1Pmesh P1Pmesh;
typedef P1Pmesh *pP1Pmesh;
struct _P1Pmesh {
    psP1P *sbf;
    int ndof;
    char *type;
    int type2;
    int *globaldofs;
};

typedef struct _PD0mesh PD0mesh;
typedef PD0mesh *pPD0mesh;
struct _PD0mesh {
    psPD0 *sbf;
    int ndof;
    int nt;
    char *type;
    int type2;
};

typedef struct _PD1mesh PD1mesh;
typedef PD1mesh *pPD1mesh;
struct _PD1mesh {
    psPD1 *sbf;
    psPD1b *sbfb;
    pTrianglep *T;
    int *triangles;
    int **lensbf;
    int *dofsb;
    int ndofb;
    int ndof;
    int nt; 
    char *type;
    int type2;
};

typedef struct _RWGmesh RWGmesh;                                                                              
typedef RWGmesh *pRWGmesh;                                                                                    
struct _RWGmesh {                                                                                             
    psRWG *sbf;                                                                                               
    int nedge;                                                                                                
    char *type;                                                                                               
    int type2;
    int *sing1;
    int *sing2;
    int *sing3;
    int *sing4;
    int *idx_sing;
    int sz_sing;
    int *idx1;
    int sz_idx;  
}; 

typedef struct _BCmesh BCmesh;                                                                              
typedef BCmesh *pBCmesh;                                                                                    
struct _BCmesh {                                                                                             
    psBC *sbf;                                                                                               
    int nedge;                                                                                                
    char *type;                                                                                                                                                                                            
};

typedef struct _Meshes Meshes;
typedef Meshes *pMeshes;
struct _Meshes {
    pP0Pmesh P0P;
    pP1Pmesh P1P;
    pPD0mesh PD0;
    pPD1mesh PD1;
    pRWGmesh RWG;
    int ntp;
    int ndp;
    int ndb;
    int nvp;
    int nvb;
    int nep;
    int neb;
    int ndofb;
    double **Vp;
    double **Vb;
    char *trialbf;
    char *testbf;
};

typedef struct _Mesh Mesh;                                                                                
typedef Mesh *pMesh;                                                                                      
struct _Mesh {                                                                                                                                                                                        
    int ntp;                                                                                                  
    int ntb;
    int ndp;                                                                                                  
    int ndb;                                                                                                  
    int nvp;                                                                                                  
    int nvb;                                                                                                  
    int nep;                                                                                                  
    int neb;                                                                                                  
    int ndofb;                                                                                                
    double **Vp;                                                                                              
    double **Vb;  
    double **Tp;                                                                                              
    double **Tb;
    int **Dp;                                                                                              
    int **Db;
    int **Dq;
    int **Ep;                                                                                              
    int **Eb;  
    int **Ep2;                                                                                                 
    int **Eb2;
};

typedef struct _Space Space;                                                                                  
typedef Space *pSpace;                                                                                        
struct _Space {                                                                                               
    pP0Pmesh P0P;                                                                                             
    pP1Pmesh P1P;                                                                                             
    pPD0mesh PD0;                                                                                             
    pPD1mesh PD1;                                                                                             
    pRWGmesh RWG;                                                                                             
    pBCmesh BC;
    pMesh M;  
    int space;
};
void delete_Triangle(pTriangle T);

void delete_Trianglep(pTrianglep T);

void delete_Quadrilateral(pQuadrilateral Q);

void delete_sP0P(psP0P s);

void delete_sP1P(psP1P s);

void delete_sPD0(psPD0 s);

void delete_sPD1(psPD1 s);

void delete_P0Pmesh(pP0Pmesh M);

void delete_P1Pmesh(pP1Pmesh M);

void delete_PD0mesh(pPD0mesh M);

void delete_PD1mesh(pPD1mesh M);

void del_mesh(pMeshes M);

pTriangle new_Triangle(double *vertices, int *dofs);

pTrianglep new_Trianglep(double *vertices, int *dofs, double **verticesb, int **dofsb, int nt);

pQuadrilateral new_Quadrilateral(double *vertices, int *dofs, int *dofsp);

psP0P new_sP0P(double *vertices, int *dofs);

psP1P new_sP1P(double **vertices, int **dofs, int len, int *num);

psPD0 new_sPD0(double **vertices, int **dofs, int **dofsp, int len, int *num);

psPD1 new_sPD1(double **vertices, double **verticesb, int **dofsp, int ***dofsb, int len, int *num);

psRWG new_sRWG(double *v1, double *v2, int *dofs1, int *dofs2, int *E, int i);

void complete_sRWG( pTriangle T1, pTriangle T2, int sz, int *E );

pP0Pmesh new_P0Pmesh(double **triangles, int **dofs, int nt);

pP1Pmesh new_P1Pmesh(double **vertices, int **dofs, int ndof, int nt);

pPD0mesh new_PD0mesh(double **verticesp, double **verticesb, int **dofsp, int **dofsb, int **dofsq, int ndof,  int ntb, int ntp);

pPD1mesh new_PD1mesh(double **vertices, double **verticesb, int **dofs, int **dofsb, int ndof, int nt, int ndofb, int ntb);

pRWGmesh new_RWGmesh(double **triangles, int **edges, int **dofs,int nt, int ne);

pMeshes new_Mesh(double **Vp, double **Vb, int **Dp, int **Db, int **Dq, double **Tp, int **Ep, int **Eb, int size_of_Dp, int size_of_Db, int size_of_Vp, int size_of_Vb, int size_of_Edgep, int size_of_Edgeb, int *opt);

pSpace new_Space(double **Vp, double **Vb, int **Dp, int **Db, int **Dq, double **Tp, double **Tb, int **Ep, int **Eb, int **Ep2, int **Eb2, int size_of_Dp, int size_of_Db, int size_of_Vp, int size_of_Vb, int size_of_Edgep, int size_of_Edgeb, int opt);

int inlist(int e, int *list, int len);

int ntriangles(pPD1mesh PD1, int *notriangles, int len, int lnt);

pMeshes sub_Mesh_P0D(pMeshes Mesh, char *type, int *indexp, int startp, int sizep);

pMeshes sub_Mesh_P1D(pMeshes Mesh, char *type, int *indexp, int *notriangle, int startp, int sizep, int lnt);

pMeshes sub_Mesh_P1D2(pMeshes Mesh, char *type, int *indexp, int *indexb, int *notriangle, int startp, int sizep, int startb, int sizeb, int lnt);
#endif
