#ifndef ACA_h
#define ACA_h

#include "block_cluster_tree_elements.h"
#include "assembly_hmat.h"
void ACA_it( gsl_matrix_complex *a, gsl_matrix_complex *b, int it, int opt, int idx );

prkmatrix ACA(
              gsl_matrix_complex *M,
              pcluster row,
              // wavenumber
              pcluster col,
              int *index,
              double kappa,
              // Order of quadrature
              int *order,
              char *trialbf,
              char *testbf,
              int rows, int cols, double eps);

prkmatrix ACA_hlp(                                                                                                
              gsl_matrix_complex *M,                                                                          
              pcluster row,                                                                                   
              // wavenumber                                                                                   
              pcluster col,                                                                                   
              int *index,                                                                                     
              double kappa,                                                                                   
              // Order of quadrature                                                                          
              int *order,                                                                                      
              char *trialbf,                                                                                  
              char *testbf,                                                                                   
              int rows, int cols, double eps);   
#endif
