#include "functions.h"

// LU decomoposition of a general matrix
void dgetrf_(int* M, int *N, double* A, int* lda, int* IPIV, int* INFO);
// generate inverse of a matrix given its LU decomposition
void dgetri_(int* N, double* A, int* lda, int* IPIV, double* WORK, int* lwork, int* INFO);

void inverse(int N, double *A)
{

    int IPIV[N+1];
    int LWORK = N*N;
    double WORK[LWORK];
    int INFO;

    dgetrf_(&N,&N,A,&N,IPIV,&INFO);
    dgetri_(&N,A,&N,IPIV,WORK,&LWORK,&INFO);

}

int check_if_equal( int *v1, int *v2 )
{
    int i, j;
    int t = 0;
    for(i=0; i<3; i++)
    {
        for(j=0; j<3; j++)
        {
            if( v1[i]==v2[j] )
            {
                t++;
            }
        }
    }
    return t;
/*
    if(v1[0]==v2[0]||v1[0]==v2[1]||v1[0]==v2[2]||v1[1]==v2[0]||v1[1]==v2[1]||v1[1]==v2[2]||v1[2]==v2[0]||v1[2]==v2[1]||v1[2]==v2[2])
    {
        return 1;
    }
    else
    {
        return 0;
    }*/
}

double norm(double *x)
{
    return sqrt(pow(x[0],2)+pow(x[1],2)+pow(x[2],2));
}

void vec_prod(int sz,double *a, double *b,double *c, double *res)
{
    int i;
    
    for (i = 0; i < sz; i++)
    {
        res[i] = a[i]*b[i]*c[i];
    }
    
}

void vec_prod_maxwell(int sz,double *a, double *b,double *phi, double *res)                                             
{                                                                                                             
    int i;                                                                                                    
                                                                                                              
    for (i = 0; i < sz; i++)                                                                                  
    {                                                                                                         
        res[3*i] = a[i]*b[i]*phi[3*i];
        res[3*i+1] = a[i]*b[i]*phi[3*i+1]; 
        res[3*i+2] = a[i]*b[i]*phi[3*i+2]; 

    }                                                                                                         
                                                                                                              
}

gsl_vector_complex *vec_prod2(int sz,double *a, double *b,double *c)                                             
{                                                                                                             
    int i;
    gsl_vector_complex *res = gsl_vector_complex_alloc(sz);                   
                                                                                                              
    for (i = 0; i < sz; i++)                                                                                  
    {                                                                                                         
        gsl_vector_complex_set(res,i,gsl_complex_rect(a[i]*b[i]*c[i],0));                                                                              
    }                                                    
    return res;                                                                                                                                                    
}

void matrix_prod(int s_r_A, int s_c_A, int s_r_B,int s_c_B, double *A, double *B)
{
    
    int i, j, k, index=0;
    double res[s_r_A*s_c_B], sum;
    
    for(i=0; i < s_r_A; i++)
    {
        for(j=0; j < s_c_B; j++)
        {
            sum = 0.0;
            
            for(k = 0; k < s_c_A; k++)
            {
                sum+=A[i*s_c_A+k]*B[k*s_c_B+j];
            }
            
            res[index] = sum;
            index++;
        }
    }
    
    for(i=0; i<s_r_A*s_c_B; i++)
    {
        B[i] = res[i];
    }
    
}

double dot_prod(int sz, int idx1, int idx2, double *a, double *b)
{
    int  i;
    double res = 0;
    
    
    for (i=0; i<sz; i++)
    {
        res += a[i+idx1]*b[3*idx2+i];
        
    }
    
    return res;
}

double dot_prod2(int sz, double *a, double *b)                                             
{                                                                                                             
    int  i;                                                                                                   
    double res = 0;                                                                                           
                                                                                                              
                                                                                                              
    for (i=0; i<sz; i++)                                                                                      
    {                                                                                                         
        res += a[3*i]*b[3*i]+a[3*i+1]*b[3*i+1]+a[3*i+2]*b[3*i+2];                                                                         
                                                                                                              
    }                                                                                                         
                                                                                                              
    return res;                                                                                               
} 

void mat_vec_prod(int s_r_A, int s_c_A,int idx, double *A, double *b,double *res)
{
    int i;
    double aux[s_r_A];
    
    
    for (i=0; i< s_r_A; i++)
    {
        aux[i] = dot_prod(s_c_A,i*s_c_A,idx,A,b);
    }
    
    for (i=0; i< s_r_A*s_c_A; i++)
    {
        res[i]=aux[i];
    }
    
}


void transpose( int row, int col, double *src, double *dst)
{
    int i, j, index = 0;
    double aux[col*row];
    
    for (i=0; i<col; i++)
    {
        for (j=0; j<row; j++)
        {
            aux[index]=src[j*col+i];
            index++;
        }
    }

   for (i=0; i<col*row; i++)
    {
        dst[i]=aux[i];
    }    
}


void sust_add_vec(double *a, double *b, int opt, int index1, int index2, double *res)
{
    int size_b = 3, i;
    
    if (opt==0)
    {
        for (i=0; i<size_b ; i++)
        {
            res[i]=a[3*index1+i]+b[3*index2+i];
        }
    }
    else
    {
        for (i=0; i<size_b ; i++)
        {
            res[i] = a[3*index1+i]-b[3*index2+i];
        }
    }
}

void sust_add_vec_mat(int col_A, int opt,double *A, double *b, int row_A, double *res)
{
    int          i, j;
    double aux[col_A];
    
    for (i = 0; i < row_A ; i++)
    {
        sust_add_vec(A, b, opt, i, 0, aux);
        
        for(j = 0; j < col_A; j++)
        {
            res[i*col_A+j] = aux[j];
        }
    }
}

double sum_v(double *v, int sz)
{
    double sum = 0.0;
    int i;

    for(i = 0; i < sz; i++)
    {
        sum+=v[i];
    }
    
    return sum;
}

int is_elmt(int e,int *array, int sz)
{
    int i;

    for(i = 0; i<sz; i++)
    {
        if(e==array[i])
        {
            return 1;
        }
    }
    return 0;
}
