#include "hmat_arithmetics.h"

#define MAX(x, y) (((x) > (y)) ? (x) : (y))
#define MIN(x, y) (((x) < (y)) ? (x) : (y))

gsl_complex mult_mat_vec(gsl_vector_complex *a, gsl_vector_complex *b, int *in)
{
    int i;
    gsl_complex res = gsl_complex_rect(0,0); 
    for(i=0; i<a->size; i++)
    {
        res = gsl_complex_add(res, gsl_complex_mul(gsl_vector_complex_get(a,i), gsl_vector_complex_get(b,in[i]))); 
    }
    return res;
}

gsl_complex mult_mat_vec2(gsl_vector_complex *a, gsl_vector_complex *b)                               
{                                                                                                             
    int i;                                                                                                    
    gsl_complex res = gsl_complex_rect(0,0);                                                                  
    for(i=0; i<a->size; i++)                                                                                  
    {                                                                                                         
        res = gsl_complex_add(res, gsl_complex_mul(gsl_vector_complex_get(a,i), gsl_vector_complex_get(b,i)));
    }                                                                                                         
    return res;                                                                                               
} 

void full_mult(const pfullmatrix f, pcluster rows, pcluster cols, gsl_vector_complex *in, gsl_vector_complex *out)
{ 
    int i, j;
    gsl_complex one = gsl_complex_rect(1,0);
    gsl_blas_zgemv (CblasNoTrans, one, f->e, in, gsl_complex_rect(0,0), out);

}

void low_rank_mult(const prkmatrix r, pcluster rows, pcluster cols, gsl_vector_complex *in, gsl_vector_complex *out)
{
    gsl_vector_complex * auxres = gsl_vector_complex_alloc (r->b->size1);
    gsl_complex one = gsl_complex_rect(1,0);
    gsl_complex zero = gsl_complex_rect(0,0);
    gsl_blas_zgemv(CblasNoTrans,one,r->b, in, zero,auxres);
    gsl_blas_zgemv(CblasTrans,one,r->a, auxres, zero,out);
    gsl_vector_complex_free(auxres);
}

void full_mult2(const pfullmatrix f, gsl_vector_complex *in, gsl_vector_complex *out, int *inv, int *outv)
{                                                                                                             
    int i;                  
    for(i=0; i<f->e->size1; i++)                                                                                  
    { 
        gsl_vector_complex_view a = gsl_matrix_complex_row(f->e,i);
        gsl_vector_complex_set(out, outv[i], mult_mat_vec(&a.vector, in, inv));
    }                                                                                                                                        
}                                                                                                             
                                                                                                              
void low_rank_mult2(const prkmatrix r, gsl_vector_complex *in, gsl_vector_complex *out, int *inv, int *outv)
{                                                                                                             
    gsl_vector_complex * auxres = gsl_vector_complex_alloc (r->b->size1);                                      
    int i;                                                                                                    
    for(i=0; i<r->b->size1; i++)                                                                              
    {                                                                                                         
        gsl_vector_complex_view b = gsl_matrix_complex_row(r->b,i);                                           
        gsl_vector_complex_set(auxres, i, mult_mat_vec(&b.vector, in, inv));                                       
    }                                                   
    for(i=0; i<r->a->size2; i++)                                                                              
    {                                                                                                         
        gsl_vector_complex_view a = gsl_matrix_complex_column(r->a,i);                                           
        gsl_vector_complex_set(out, outv[i], mult_mat_vec2(&a.vector, auxres));                                          
    }
    gsl_vector_complex_free(auxres);                                                                          
}  

void set_vec(psupermatrix s, const gsl_vector_complex * x)
{
    int j;
    for(j = 0; j< s->in->size; j++)                                                                  
    {                                                                                                     
         gsl_vector_complex_set(s->in,j,gsl_vector_complex_get(x,s->hfi[j]));                  
    } 
}


void mat_vec(const psupermataux s, const gsl_vector_complex * x,  gsl_vector_complex * res, double *t)                    
{                                                                                                                                                                                      
    int i,j;
    gsl_vector_complex *ones = gsl_vector_complex_alloc(s->size);
    struct timeval st, e;

    gettimeofday(&st, NULL);
    for(i=0;i<s->size;i++)                                                                                    
    {                                                                                                         
        set_vec(s->s[i],x);                   
        gsl_vector_complex_set(ones,i,gsl_complex_rect(1,0));
    }
    gsl_matrix_complex *A = gsl_matrix_complex_calloc(res->size, s->size);
    gettimeofday(&e, NULL);                                                                               
    t[0]+= (e.tv_sec - st.tv_sec) * 1000.0;                                               
    t[0]+= (e.tv_usec - st.tv_usec) / 1000.0;
    
    for(i=0;i<s->size;i++)                                                                                    
    {                
        gsl_vector_complex_view R = gsl_matrix_complex_column(A,i);                                                                                           
        gettimeofday(&st, NULL);
        if(s->s[i]->full == NULL)                                                                             
        {                                                                                                     
            low_rank_mult(s->s[i]->rk, s->s[i]->rows, s->s[i]->cols, s->s[i]->in, s->s[i]->out);              
        }                                                                                                     
        else                                                                                                  
        {                                                                                                     
            full_mult(s->s[i]->full, s->s[i]->rows, s->s[i]->cols, s->s[i]->in, s->s[i]->out);                
        }    
        gettimeofday(&e, NULL);
        t[1]+= (e.tv_sec - st.tv_sec) * 1000.0; 
        t[1]+= (e.tv_usec - st.tv_usec) / 1000.0; 
        gettimeofday(&st, NULL); 
        for(j=0;j<s->s[i]->out->size; j++)                                                                    
        {                                                                                                     
            gsl_vector_complex_set(&R.vector,s->s[i]->hfo[j], gsl_vector_complex_get(s->s[i]->out,j));                                                
        }                        
        gettimeofday(&e, NULL);
        t[0]+= (e.tv_sec - st.tv_sec) * 1000.0;                                                          
        t[0]+= (e.tv_usec - st.tv_usec) / 1000.0; 
    } 

    gettimeofday(&st, NULL);
    gsl_blas_zgemv (CblasNoTrans, gsl_complex_rect(1,0), A, ones, gsl_complex_rect(0,0), res); 
    gettimeofday(&e, NULL);                                                                               
    t[0]+= (e.tv_sec - st.tv_sec) * 1000.0;                                                          
    t[0]+= (e.tv_usec - st.tv_usec) / 1000.0;
    gsl_vector_complex_free(ones);
    gsl_matrix_complex_free(A);
}


void mat_vec2(const psupermataux s, const gsl_vector_complex * x,  gsl_vector_complex * res, double *t)        
{                                                                                                             
    int i,j;                                                                                                  
    gsl_vector_complex *ones = gsl_vector_complex_alloc(s->size);                                             
    struct timeval st, e;                                                                                     
                                                                                                              
    gettimeofday(&st, NULL);                                                                                  
    for(i=0;i<s->size;i++)                                                                                    
    {                                                                                                         
        set_vec(s->s[i],x);                                                                                   
        gsl_vector_complex_set(ones,i,gsl_complex_rect(1,0));                                                 
    }                                                                                                         
    //gsl_matrix_complex *A = gsl_matrix_complex_calloc(res->size, s->size);                                    
    gettimeofday(&e, NULL);                                                                                   
    t[0]+= (e.tv_sec - st.tv_sec) * 1000.0;                                                                   
    t[0]+= (e.tv_usec - st.tv_usec) / 1000.0;                                                                 
                                                                                                              
    for(i=0;i<s->size;i++)                                                                                    
    {                                                                                                         
        gettimeofday(&st, NULL);                                                                              
        if(s->s[i]->full == NULL)                                                                             
        {                                                                                                     
            low_rank_mult2(s->s[i]->rk, x, res, s->s[i]->hfi, s->s[i]->hfo);              
        }                                                                                                     
        else                                                                                                  
        {                                                                                                     
            full_mult2(s->s[i]->full, x, res, s->s[i]->hfi, s->s[i]->hfo);                
        }                                                                                                     
        gettimeofday(&e, NULL);                                                                               
        t[1]+= (e.tv_sec - st.tv_sec) * 1000.0;                                                               
        t[1]+= (e.tv_usec - st.tv_usec) / 1000.0;                                                             
        gettimeofday(&st, NULL);                                                                                                                                                      
        t[0]+= (e.tv_sec - st.tv_sec) * 1000.0;                                                               
        t[0]+= (e.tv_usec - st.tv_usec) / 1000.0;                                                             
    }                                                                                                         
                                                                                                              
    /*gettimeofday(&st, NULL);                                                                                  
    gsl_blas_zgemv (CblasNoTrans, gsl_complex_rect(1,0), A, ones, gsl_complex_rect(0,0), res);                
    gettimeofday(&e, NULL);*/                                                                                   
    t[0]+= (e.tv_sec - st.tv_sec) * 1000.0;                                                                   
    t[0]+= (e.tv_usec - st.tv_usec) / 1000.0;                                                                 
    //gsl_vector_complex_free(ones);                                                                            
    //gsl_matrix_complex_free(A);                                                                               
}  
void mat_mat(psupermataux A, gsl_matrix_complex *B, double *t)
{
    gsl_vector_complex_view a1;
    gsl_vector_complex *a2 = gsl_vector_complex_calloc(B->size1);
    int i;

    for(i=0;i<B->size2;i++)
    {
        a1 = gsl_matrix_complex_column(B,i);
        mat_vec(A, &a1.vector,a2, t);
        gsl_matrix_complex_set_col(B,i,a2);
    }
    gsl_vector_complex_free(a2); 

}
