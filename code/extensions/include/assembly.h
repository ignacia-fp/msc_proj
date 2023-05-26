#ifndef fasthelm_h
#define fasthelm_h
#include "block_cluster_tree_elements.h"
#include <time.h> 


void assemble_SLD(
                  // Input:
                  // V: vertices of each element
                  pMeshes Mesh,
                  // wavenumber
                  double kappa,
                  // Order of quadrature
                  int *order,
                  // Output: complex and real part
                  double *M1,
                  double *M2,
                  char *trialbf,
                  char *testbf
                  );

void assemble_SLP(
                  // Input:
                  // Space data
                  pMeshes Mesh,
                  // wavenumber
                  double kappa,
                  // Order of quadrature
                  int *order,
                  // Output: complex and real part
                  double *M1,
                  double *M2,
                  char *trialbf,
                  char *testbf
                  );

void assemble_HP(
                 // Input:
                 // Space data
                 pMeshes Mesh,
                 // wavenumber
                 double kappa,
                 // Order of quadrature
                 int *order,
                 // Output: complex and real part
                 double *M1,
                 double *M2,
                 char *trialbf,
                 char *testbf
                 );

void assemble_HP2(                                                                                            
                 // Input:                                                                                    
                 // Space data                                                                                
                 pMeshes Mesh,                                                                                
                 // wavenumber                                                                                
                 double kappa,                                                                                
                 // Order of quadrature                                                                       
                 int *order,                                                                                  
                 // Output: complex and real part                                                             
                 double *M1,                                                                                  
                 double *M2,                                                                                  
                 char *trialbf,                                                                               
                 char *testbf                                                                                 
                 );

void assemble_maxwell_weakly(// Input:                                                                                    
                 // Space data                                                                                
                 pMeshes Mesh,                                                                                
                 // wavenumber                                                                                
                 double kappa,                                                                                
                 // Order of quadrature                                                                       
                 int *order,                                                                                  
                 // Output: complex and real part                                                             
                 double *M1,                                                                                  
                 double *M2,                                                                                  
                 char *trialbf,                                                                               
                 char *testbf );

void assemble_maxwell_T(// Input:                                                                                    
                 // Space data                                                                                
                 pSpace trial,                                                                                
                 pSpace test,                                                                                 
                 // wavenumber                                                                                
                 double kappa,                                                                                
                 // Order of quadrature                                                                       
                 int *order,                                                                                  
                 // Output: complex and real part                                                             
                 double *M1,                                                                                  
                 double *M2,                                                                                  
                 double *M3,                                                                                  
                 double *M4 );

void assemble_maxwell_hypersingular(// Input:                                                                                    
                 // Space data                                                                                
                 pSpace trial,                                                                                
                 pSpace test,                                                                                 
                 // wavenumber                                                                                
                 double kappa,                                                                                
                 // Order of quadrature                                                                       
                 int *order,                                                                                  
                 // Output: complex and real part                                                             
                 double *M1,                                                                                  
                 double *M2 );

void assemble_maxwell_identity(// Input:                                                                                    
                 // Space data                                                                                
                 pMeshes Mesh,                                                                                
                 // wavenumber                                                                                
                 double kappa,                                                                                
                 // Order of quadrature                                                                       
                 int *order,                                                                                  
                 // Output: complex and real part                                                             
                 double *M1,                                                                                  
                 char *trialbf,                                                                               
                 char *testbf );

void assemble_maxwell_laplace_beltrami(// Input:                                                                                    
                 // Space data                                                                                
                 pMeshes Mesh,                                                                                
                 // wavenumber                                                                                
                 double kappa,                                                                                
                 // Order of quadrature                                                                       
                 int *order,                                                                                  
                 // Output: complex and real part                                                             
                 double *M1,                                                                                  
                 char *trialbf,                                                                               
                 char *testbf );

gsl_complex assemble_P1_HP_entries(int *order, double kappa, int i, int j, double *sing, pMeshes Mesh);

gsl_complex s_part(double kappa, int *order, pMeshes Mesh, int i, int k, int j, int l, int opt);

void assembly_P1d1( pMeshes Mesh, int i, int dof, int k, int *sz, double *P0, double *IE, double *eta, double *chi);

void assembly_P1d2(pMeshes Mesh, int i, int k, int *sz, double *P0, double *IE, double *eta, double *chi);

void assembly_P1d3(pMeshes Mesh, int i, int dof, int k, int *sz, double *P0, double *IE, double *eta, double *chi);

gsl_complex assemble_P1_HP_entriesd(int *order, double kappa, int i, int j, int idx, int idy, double *sing, pMeshes Mesh, int dofx, int dofy);

void project_mat(                                                                                             
                 // Input:                                                                                    
                 // Space data                                                                                
                 pMeshes Mesh,                                                                                                                                                                
                 // Output: complex and real part                                                             
                 double *Matr,                                                                                
                 double *Mati,                                                                                
                 double *M1,                                                                                  
                 double *M2,                                                                                  
                 char *trialbf,                                                                               
                 char *testbf                                                                                 
                 );
void assemble_IDP(
                  // Input:
                  // V: vertices of each element
                  pMeshes Mesh,
                  // wavenumber
                  // Order of quadrature
                  int order,
                  // Output: complex and real part
                  double* M,
                  char *trialbf,
                  char *testbf
                  );

void assemble_IDD(
                  // Input:
                  pMeshes Mesh,
                  // wavenumber
                  // Order of quadrature
                  int order,
                  // Output: complex and real part
                  double* M,
                  char *trialbf,
                  char *testbf
                  );

void assemble_IDDP(pMeshes Mesh,                                                                              
                  // wavenumber                                                                               
                  // Order of quadrature                                                                      
                  int order,                                                                                  
                  // Output                                                                                   
                  double* M,                                                                                  
                  char *trialbf,                                                                              
                  char *testbf);

void assemble_IDM(
                  pMeshes Mesh,
                  // wavenumber
                  // Order of quadrature
                  int order,
                  // Output
                  double* M,
                  char *trialbf,
                  char *testbf
                  );

void assemble_IDDB(
                  pMeshes Mesh,
                  // Order of quadrature
                  int order,
                  // Output
                  double* M,
                  char *trialbf,
                  char *testbf
                  );

void average_matrix_full(
                   pMeshes Mesh,
                   // Output
                   double* M,
                   int opt
                   );


void projection_matrix(                                                                                       
                   pMeshes Mesh,                                                                              
                   // Order of quadrature                                                                     
                   // Output                                                                                                                                                                
                   double* M,                                                                                 
                   int opt                                                                                    
                   );


#endif /* functions_h */




