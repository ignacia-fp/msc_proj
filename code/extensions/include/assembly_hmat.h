#ifndef fasthelm_hmat_h
#define fasthelm_hmat_h
#include "assembly.h"


gsl_matrix_complex *assemble_SLD_f(
                                  // Input:
                                  // V: vertices of each element
                                  pcluster row,
                                  // wavenumber
                                  pcluster col,
                                  int *index,
                                  double kappa,
                                  // Order of quadrature
                                  int *order,
                                  // Output: complex and real part
                                  char *trialbf,
                                  char *testbf
                                  );

gsl_complex assemble_P1_HP_entries2(int *order, double kappa, int i, int j, int idx, int idy, double *sing, pMeshes Mesh1, pMeshes Mesh2, int dofx, int dofy);
                                                                                                              
gsl_complex assemble_P1_HP_entriesd2(int *order, double kappa, int i, int j, double *sing, pMeshes Mesh1, pMeshes Mesh2);

gsl_matrix_complex *assemble_HLP_lr2(gsl_matrix_complex* Mat,                                                 
                    gsl_matrix_complex* m,                                                                    
                   pMeshes Mesh1,                                                                             
                   pMeshes Mesh2,                                                                             
                    char *mode,                                                                               
                   int *order,                                                                                
                   double kappa,                                                                              
                   // Order of quadrature                                                                     
                   // Output                                                                                           
                   int urow,                                                                                  
                   int ucol,                                                                                  
                   int it                                                                                     
                  );

gsl_matrix_complex *assemble_HLP_f(gsl_matrix_complex* Mat,                                                   
                   pMeshes Mesh1,                                                                             
                   pMeshes Mesh2,                                                                             
                   int *order,                                                                                 
                   double kappa,                                                                              
                   // Order of quadrature                                                                     
                   // Output                                                                                  
                   int opt                                                                                    
                  );

void assemble_SLD_lr(
                                   // Input:
                                   // V: vertices of each element
                                   gsl_matrix_complex *M,
                                   gsl_matrix_complex *m,
                                   pcluster row,
                                   // wavenumber
                                   pcluster col,
                                   int *index,
                                   double kappa,
                                   // Order of quadrature
                                   int *order,
                                   char *trialbf,
                                   char *testbf,
                                   char *mode,
                                   int urow,
                                   int ucol,
                                   int it
                                   );

gsl_matrix_complex *average_matrix(                                                                           
                   pMeshes Mesh,                                                                              
                   // Order of quadrature                                                                     
                   // Output                                                                                  
                   int opt                                                                                    
                   );

gsl_complex assemble_P1_entries(int *order,double kappa, int i, int j, pMeshes Mesh1, pMeshes Mesh2);

gsl_complex assemble_P1_entries2(int *order,double kappa, int i, int j, pMeshes Mesh1, pMeshes Mesh2); 

gsl_matrix_complex *assemble_SLD_P1_f( 
                    gsl_matrix_complex *M,
                   pMeshes Mesh1,                                                                             
                   pMeshes Mesh2,                                                                             
                   int *order,    
                   double kappa,
                   // Order of quadrature                                                                     
                   // Output                                                                                  
                   int opt                                                                                    
                   );

void assemble_HLP_lr(                                                                                         
                    // Input:                                                                                 
                    // V: vertices of each element                                                            
                    gsl_matrix_complex *M,                                                                    
                    gsl_matrix_complex *m,                                                                    
                    pcluster row,                                                                             
                    // wavenumber                                                                             
                    pcluster col,                                                                             
                    int *index,                                                                               
                    double kappa,                                                                             
                    // Order of quadrature                                                                    
                    int *order,                                                                                
                    char *trialbf,                                                                            
                    char *testbf,                                                                             
                    char *mode,                                                                               
                    int urow,                                                                                 
                    int ucol,                                                                                 
                    int it                                                                                    
                    );
#endif 




