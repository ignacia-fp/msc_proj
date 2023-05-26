
#include "ACA.h"

#define MAX(x, y) (((x) > (y)) ? (x) : (y))
#define MIN(x, y) (((x) < (y)) ? (x) : (y))


void ACA_it( gsl_matrix_complex *a, gsl_matrix_complex *b, int it, int opt, int idx )
{
    int i;
    gsl_complex mone;
    GSL_SET_COMPLEX(&mone,-1,0);
    if(opt == 0)
    {
        gsl_vector_complex_view resa = gsl_matrix_complex_row(a,it);
        for(i=0;i < it; i++)
        {
            gsl_vector_complex_view auxa = gsl_matrix_complex_row(a,i);
            gsl_vector_complex_view auxb = gsl_matrix_complex_row(b,i);
            gsl_blas_zaxpy(gsl_complex_mul(mone,gsl_vector_complex_get(&auxb.vector,idx)),&auxa.vector,&resa.vector);
        }
    }
    else
    {
        gsl_vector_complex_view resb = gsl_matrix_complex_row(b,it);
        for(i=0;i < it; i++)
        {
            gsl_vector_complex_view auxb = gsl_matrix_complex_row(b,i);
            gsl_vector_complex_view auxa = gsl_matrix_complex_row(a,i);
        gsl_blas_zaxpy(gsl_complex_mul(mone,gsl_vector_complex_get(&auxa.vector,idx)),&auxb.vector,&resb.vector);
        }
    }
    
}


prkmatrix ACA(gsl_matrix_complex* M, pcluster row,
              // wavenumber
              pcluster col,
              int *index,
              double kappa,
              // Order of quadrature
              int *order,
              char *trialbf,
              char *testbf,
              int rows, int cols, double eps)
{
    prkmatrix r = (prkmatrix) malloc(sizeof(rkmatrix));
    gsl_complex mone;
    gsl_complex zero;
    gsl_complex max_b;
    GSL_SET_COMPLEX(&zero,0,0);
    GSL_SET_COMPLEX(&mone,-1,0);
    r->rows = rows;
    r->cols = cols;
    r->a = gsl_matrix_complex_alloc(1, rows);
    r->b = gsl_matrix_complex_alloc(1, cols);
    int *Index_a, *Index_b, it;
    gsl_complex delta;
    double err;
    int i;
    Index_a = (int*) malloc(sizeof(int));
    Index_b = (int*) malloc(sizeof(int));
    it = 0;
    Index_b[0] = 0;
    assemble_SLD_lr(M, r->b, row, col, index, kappa, order, trialbf, testbf, "row", Index_b[0], 0, 0);
    gsl_vector_complex_view auxb = gsl_matrix_complex_row(r->b, 0);
    Index_a[0] = gsl_vector_complex_max_index(&auxb.vector);
    assemble_SLD_lr(M, r->a, row, col, index, kappa, order, trialbf, testbf, "col", 0, Index_a[0], 0);
    gsl_vector_complex_view auxa = gsl_matrix_complex_row(r->a, 0);
    delta = gsl_vector_complex_get(&auxb.vector,Index_a[0]); 
    gsl_vector_complex_scale(&auxa.vector, gsl_complex_inverse(delta));
   
    err = gsl_blas_dznrm2(&auxa.vector)*gsl_blas_dznrm2(&auxb.vector);
    it++;
    
    while(err > eps && it<(MIN(rows,cols)))
    {
        gsl_vector_complex_view auxit = gsl_matrix_complex_row(r->a, it-1);
        gsl_vector_complex *auxitb = gsl_vector_complex_alloc((&auxit.vector)->size);
        gsl_vector_complex_memcpy(auxitb,&auxit.vector);
        
        for(i=0;i<it;i++)
        {
            gsl_vector_complex_set(auxitb,Index_b[i],zero);
        }
        
        Index_b = (int*)realloc(Index_b,sizeof(int)*(it+1));
        size_t sb = r->b->tda;
        r->b->block->data = (gsl_complex *)realloc(r->b->block->data, sizeof(gsl_complex)*(r->b->block->size + sb));
        r->b->block->size += sb;
        r->b->size1++;
        r->b->data = r->b->block->data;
        Index_b[it] = gsl_vector_complex_max_index(auxitb);
        gsl_vector_complex_view auxb = gsl_matrix_complex_row(r->b, it);
        assemble_SLD_lr(M, r->b, row, col, index, kappa, order, trialbf, testbf, "row", Index_b[it], 0, it);
        ACA_it( r->a, r->b, it, 1, Index_b[it] );
        max_b = gsl_vector_complex_get(&auxb.vector,gsl_vector_complex_max_index(&auxb.vector));
        
        while(GSL_REAL(max_b) == 0.0 && GSL_IMAG(max_b) == 0.0)
        {
            gsl_vector_complex_set(auxitb,Index_b[it],zero);
            Index_b[it] = gsl_vector_complex_max_index(auxitb);
            assemble_SLD_lr(M, r->b, row, col, index, kappa, order, trialbf, testbf, "row", Index_b[it], 0, it);
            ACA_it( r->a, r->b, it, 1, Index_b[it] );
            max_b = gsl_vector_complex_get(&auxb.vector,gsl_vector_complex_max_index(&auxb.vector));
        }
        
        Index_a = (int*)realloc(Index_a,sizeof(int)*(it+1));
        size_t sa = r->a->tda;
        r->a->block->data = (gsl_complex *)realloc(r->a->block->data, sizeof(gsl_complex)*(r->a->block->size + sa));
        r->a->block->size += sa;
        r->a->size1++;
        r->a->data = r->a->block->data;
        gsl_vector_complex *auxita = gsl_vector_complex_alloc((&auxb.vector)->size);
        gsl_vector_complex_memcpy(auxita,&auxb.vector);
        
        for(i=0;i<it;i++)
        {
            gsl_vector_complex_set(auxita,Index_a[i],zero);
        }
        
        Index_a[it] = gsl_vector_complex_max_index(auxita);
        assemble_SLD_lr(M, r->a, row, col, index, kappa, order, trialbf, testbf, "col", 0, Index_a[it], it);
        delta = gsl_vector_complex_get(&auxb.vector,Index_a[it]);
        gsl_vector_complex_view auxa = gsl_matrix_complex_row(r->a, it);
        ACA_it( r->a, r->b, it, 0, Index_a[it] );
        gsl_vector_complex_scale(&auxa.vector, gsl_complex_inverse(delta));
        err = gsl_blas_dznrm2(&auxa.vector)*gsl_blas_dznrm2(&auxb.vector);
        it++;
    }
    free(Index_a);
    free(Index_b);
    r->kt = it;


    return r;
}


prkmatrix ACA_hlp(gsl_matrix_complex* M, pcluster row,                                                            
              // wavenumber                                                                                   
              pcluster col,                                                                                   
              int *index,                                                                                     
              double kappa,                                                                                   
              // Order of quadrature                                                                          
              int *order,                                                                                      
              char *trialbf,                                                                                  
              char *testbf,                                                                                   
              int rows, int cols, double eps)                                                                 
{                                                                                                             
    prkmatrix r = (prkmatrix) malloc(sizeof(rkmatrix));                                                       
    gsl_complex mone;                                                                                         
    gsl_complex zero;                                                                                         
    gsl_complex max_b;                                                                                        
    GSL_SET_COMPLEX(&zero,0,0);                                                                               
    GSL_SET_COMPLEX(&mone,-1,0);                                                                              
    r->rows = rows;                                                                                           
    r->cols = cols;                                                                                           
    r->a = gsl_matrix_complex_alloc(1, rows);                                                                 
    r->b = gsl_matrix_complex_alloc(1, cols);                                                                 
    int *Index_a, *Index_b, it;                                                                               
    gsl_complex delta;                                                                                        
    double err;                                                                                               
    int i;                                                                                                    
    Index_a = (int*) malloc(sizeof(int));                                                                     
    Index_b = (int*) malloc(sizeof(int));                                                                     
    it = 0;                                                                                                   
    Index_b[0] = 0;                                                                                           
    //assemble_HLP_lr2(M, r->b, row->M, col->M, "row", order, kappa, Index_b[0], 0,0);
    assemble_HLP_lr(M, r->b, row, col, index, kappa, order, trialbf, testbf, "row", Index_b[0], 0, 0);        
    gsl_vector_complex_view auxb = gsl_matrix_complex_row(r->b, 0);                                           
    Index_a[0] = gsl_vector_complex_max_index(&auxb.vector);                                                  
    //assemble_HLP_lr2(M, r->a, row->M, col->M, "col", order, kappa, 0, Index_a[0], 0);
    assemble_HLP_lr(M, r->a, row, col, index, kappa, order, trialbf, testbf, "col", 0, Index_a[0], 0);        
    gsl_vector_complex_view auxa = gsl_matrix_complex_row(r->a, 0);                                           
    delta = gsl_vector_complex_get(&auxb.vector,Index_a[0]);                                                  
    gsl_vector_complex_scale(&auxa.vector, gsl_complex_inverse(delta));                                       
                                                                                                              
    err = gsl_blas_dznrm2(&auxa.vector)*gsl_blas_dznrm2(&auxb.vector);                                        
    it++;                                                                                                     
                                                                                                              
    while(err > eps && it<(MIN(rows,cols)))                                                                   
    {                                                                                                         
        gsl_vector_complex_view auxit = gsl_matrix_complex_row(r->a, it-1);                                   
        gsl_vector_complex *auxitb = gsl_vector_complex_alloc((&auxit.vector)->size);                         
        gsl_vector_complex_memcpy(auxitb,&auxit.vector);                                                      
                                                                                                              
        for(i=0;i<it;i++)                                                                                     
        {                                                                                                     
            gsl_vector_complex_set(auxitb,Index_b[i],zero);                                                   
        }                                                                                                     
                                                                                                              
        Index_b = (int*)realloc(Index_b,sizeof(int)*(it+1));                                                  
        size_t sb = r->b->tda;                                                                                
        r->b->block->data = (gsl_complex *)realloc(r->b->block->data, sizeof(gsl_complex)*(r->b->block->size + sb));
        r->b->block->size += sb;                                                                              
        r->b->size1++; 
        r->b->data = r->b->block->data;                                                                       
        Index_b[it] = gsl_vector_complex_max_index(auxitb);                                                   
        gsl_vector_complex_view auxb = gsl_matrix_complex_row(r->b, it);                                      
        //assemble_HLP_lr2(M, r->b, row->M, col->M, "row", order, kappa, Index_b[it], 0, it);  
        assemble_HLP_lr(M, r->b, row, col, index, kappa, order, trialbf, testbf, "row", Index_b[it], 0, it);  
        ACA_it( r->a, r->b, it, 1, Index_b[it] );                                                             
        max_b = gsl_vector_complex_get(&auxb.vector,gsl_vector_complex_max_index(&auxb.vector));              
                                                                                                              
        while(GSL_REAL(max_b) == 0.0 && GSL_IMAG(max_b) == 0.0)                                               
        {                                                                                                     
            gsl_vector_complex_set(auxitb,Index_b[it],zero);                                                  
            Index_b[it] = gsl_vector_complex_max_index(auxitb);                                               
            //assemble_HLP_lr2(M, r->b, row->M, col->M, "row", order, kappa, Index_b[it], 0, it);  
            assemble_HLP_lr(M, r->b, row, col, index, kappa, order, trialbf, testbf, "row", Index_b[it], 0, it);
            ACA_it( r->a, r->b, it, 1, Index_b[it] );                                                         
            max_b = gsl_vector_complex_get(&auxb.vector,gsl_vector_complex_max_index(&auxb.vector));          
        }                                                                                                     
                                                                                                              
        Index_a = (int*)realloc(Index_a,sizeof(int)*(it+1));                                                  
        size_t sa = r->a->tda;                                                                                
        r->a->block->data = (gsl_complex *)realloc(r->a->block->data, sizeof(gsl_complex)*(r->a->block->size + sa));
        r->a->block->size += sa;                                                                              
        r->a->size1++;                                                                                        
        r->a->data = r->a->block->data;                                                                       
        gsl_vector_complex *auxita = gsl_vector_complex_alloc((&auxb.vector)->size);                          
        gsl_vector_complex_memcpy(auxita,&auxb.vector);                                                       
                                                                                                              
        for(i=0;i<it;i++)                                                                                     
        {                                                                                                     
            gsl_vector_complex_set(auxita,Index_a[i],zero);                                                   
        }                                                                                                     
                                                                                                              
        Index_a[it] = gsl_vector_complex_max_index(auxita);                                                   
        //assemble_HLP_lr2(M, r->a, row->M, col->M, "col", order, kappa, 0, Index_a[it], it);  
        assemble_HLP_lr(M, r->a, row, col, index, kappa, order, trialbf, testbf, "col", 0, Index_a[it], it);  
        delta = gsl_vector_complex_get(&auxb.vector,Index_a[it]);                                             
        gsl_vector_complex_view auxa = gsl_matrix_complex_row(r->a, it);                                      
        ACA_it( r->a, r->b, it, 0, Index_a[it] );                                                             
        gsl_vector_complex_scale(&auxa.vector, gsl_complex_inverse(delta));                                   
        err = gsl_blas_dznrm2(&auxa.vector)*gsl_blas_dznrm2(&auxb.vector);                                    
        it++;                                                                                                 
    }                                                                                                         
    free(Index_a);                                                                                            
    free(Index_b);                                                                                            
    r->kt = it;                                                                                               
                                                                                                              
                                                                                                              
    return r;                                                                                                 
}
