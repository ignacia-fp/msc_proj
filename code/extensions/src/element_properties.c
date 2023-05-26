#include "element_properties.h"
#define PI 3.14159265358979323846
//

int compare_vert(double *v1, double *v2)
{
    int i;
    
    for( i = 0; i < 3; i++)
    {
        if( v1[i] != v2[i] )
        {
            return 1;
        }
    }
    return 0;
}

double len_edge(double *v1, double *v2)
{
    int i;
    double z[3];

    for ( i = 0; i < 3; i++ )
    {
        z[i] = v1[i] - v2[i];
    }
    
    return norm(z);
}

double deriv_t( double x1, double x2, double x3, int opt )
{
    if (opt == 0)
    {
        return x2-x1;
        
    }
    else
    {
        return x3-x1;
    }
}

double deriv_q( double x1, double x2, double x3,double x4,double eta, double chi, int opt )
{
    if (opt == 0)
    {
        return 0.25*((x1-x2)*(eta-1)+(x3-x4)*(1+eta));
        
    }
    else
    {
        return 0.25*((x1-x4)*(chi-1)+(x3-x2)*(1+chi));
    }
}

double area( double *M)
{
    double dxeta, dxchi, dyeta, dychi, dzeta, dzchi, point[3];

    dxeta = deriv_t( M[0], M[3], M[6], 1 );
    dxchi = deriv_t( M[0], M[3], M[6], 0 );
    dyeta = deriv_t( M[1], M[4], M[7], 1 );
    dychi = deriv_t( M[1], M[4], M[7], 0 );
    dzeta = deriv_t( M[2], M[5], M[8], 1 );
    dzchi = deriv_t( M[2], M[5], M[8], 0 );
    point[0] = dychi*dzeta-dyeta*dzchi;
    point[1] = -(dxchi*dzeta-dxeta*dzchi);
    point[2] = dxchi*dyeta-dxeta*dychi;
    
    return norm(point);
}

void evaluate_kernel( double *x, double *y, double kappa, int npoints_x, int npoints_y,double *di, double *dj, double *res )
{
    int            i,j,m;
    double          z[3];
    double          nrm;
    double complex res_j;
    
    res[0] = 0.0;
    res[1] = 0.0;
    
    for( j = 0; j < npoints_y; j++ )
    {
        res_j = 0.0 + 0.0*I;
        for( i = 0; i < npoints_x; i++ )
        {
            for(m = 0; m < 3; m++)
            {
                z[m] = x[m+3*i] - y[m+3*j];
            }
            nrm = norm(z);
            res_j += cexp((1*I)*kappa*nrm)/(nrm*4.0*PI)*di[i];
        }
        res[0] += creal(res_j*dj[j]);
        res[1] += cimag(res_j*dj[j]);

    }
}

void evaluate_kernel_maxwell( double *x, double *y, double kappa, int npoints_x, int npoints_y,double *di, double *dj, double *res )
{                                                                                                             
    int            i,j,m;                                                                                     
    double          z[3];                                                                                     
    double          nrm;                                                                                      
    double complex res_j;   
    double aux;                                                                                  
                                                                                                              
    res[0] = 0.0;                                                                                             
    res[1] = 0.0;                                                                                             
                                                                                                              
    for( j = 0; j < npoints_y; j++ )                                                                          
    {                                                                                                         
        res_j = 0.0 + 0.0*I;                                                                                  
        for( i = 0; i < npoints_x; i++ )                                                                      
        {                                                                                                     
            for(m = 0; m < 3; m++)                                                                            
            {                                                                                                 
                z[m] = x[m+3*i] - y[m+3*j];                                                                   
            }                                                                                                 
            nrm = norm(z);
            aux = di[3*i]*dj[3*j]+di[3*i+1]*dj[3*j+1]+di[3*i+2]*dj[3*j+2];                                                                              

            res_j += cexp((1*I)*kappa*nrm)/(nrm*4.0*PI)*aux; 
        }                                                                                                     
        res[0] += creal(res_j);                                                                         
        res[1] += cimag(res_j);    



    }    

}

void evaluate_kernel_maxwell_T( double *x, double *y, double kappa, int npoints_x, int npoints_y,double *di, double *dj, double *di2, double *dj2, double *cres1, double *cres2 )
{                                                                                                             
    int            i,j,m;                                                                                     
    double          z[3];                                                                                     
    double          nrm;                                                                                      
    double complex res_j, res_j2;                                                                                     
    double aux;                                                                                               
    double complex cone = 0.0+1.0*I;                   
    double complex res1, res2, res3;
    res1 = 0.0+ 0.0*I;                                                                                             
    res2 = 0.0+ 0.0*I;                                                                                  

    for( j = 0; j < npoints_y; j++ )                                                                          
    {                                                                                                         
        res_j = 0.0 + 0.0*I;   
        res_j2 = 0.0 + 0.0*I; 
        for( i = 0; i < npoints_x; i++ )                                                                      
        {                                                                                                     
            for(m = 0; m < 3; m++)                                                                            
            {                                                                                                 
                z[m] = x[m+3*i] - y[m+3*j];                                                                   
            }                                                                                                 
            nrm = norm(z);                                                                                    
            aux = di2[3*i]*dj2[3*j]+di2[3*i+1]*dj2[3*j+1]+di2[3*i+2]*dj2[3*j+2];                                    
            res_j += cexp((1*I)*kappa*nrm)/(nrm*4.0*PI)*di[i];                                                                      
            res_j2 += cexp((1*I)*kappa*nrm)/(nrm*4.0*PI)*aux;                                                  
        }                                                                                                     
        res1 += res_j*dj[j];                                                                            
        res2 += res_j2;


    }                                                                                                         
   cres1[0] = creal(res1);
   cres1[1] = cimag(res1);
   cres2[0] = creal(res2);                                                                                    
   cres2[1] = cimag(res2);
}  
gsl_complex HP_kern( double kappa, double *weights_x, double *weights_y, double *phii, double *phij, double *points_x, double *points_y, double *IE_x, double *IE_y, gsl_vector *n_x, gsl_vector *n_y, gsl_vector *curl_x, gsl_vector *curl_y, int sz_x, int sz_y )
{
    int            i,j,m;
    double          z[3];
    double          nrm;
    double resn, resc;
    double complex res, resj;
    
    
    res = 0.0 + 0.0*I;
    gsl_blas_ddot(n_x, n_y, &resn);                                                                   
    gsl_blas_ddot(curl_x, curl_y, &resc);  
    for( j = 0; j < sz_y; j++ )
    {
        resj = 0.0 + 0.0*I;
        for( i = 0; i < sz_x; i++ )
        {
            for(m = 0; m < 3; m++)
            {
                z[m] = points_x[m+3*i] - points_y[m+3*j];
            }
            nrm = norm(z);
            resj += cexp((1*I)*kappa*nrm)/(nrm*4.0*PI)*(weights_x[i]*IE_x[i]*(resc-kappa*kappa*resn*phii[i]*phij[j]));
          
        }
       
        res+= resj*weights_y[j]*IE_y[j];   
    
    }
    return gsl_complex_rect(creal(res), cimag(res));
}

gsl_complex HP_kern2( double kappa, double *weights_x, double *weights_y, double *phii, double *phij, double *points_x, double *points_y, double *IE_x, double *IE_y, gsl_vector *n_x, gsl_vector *n_y, gsl_matrix *curl_x, gsl_matrix *curl_y, int sz_x, int sz_y )
{                                                                                                             
    int            i,j,m;                                                                                     
    double          z[3];                                                                                     
    double          nrm;                                                                                      
    double resn, resc;                                                                                        
    double complex res, resj;                                                                                 
                                                                                                              
    res = 0.0 + 0.0*I;                                
    gsl_blas_ddot(n_x, n_y, &resn);

    for( j = 0; j < sz_y; j++ )                                                                               
    {                                                                                                         
        resj = 0.0 + 0.0*I; 
        gsl_vector_view curl1 = gsl_matrix_row(curl_y,j);  
        for( i = 0; i < sz_x; i++ )                                                                           
        {                    
            gsl_vector_view curl2 = gsl_matrix_row(curl_x,i); 
            gsl_blas_ddot(&curl1.vector, &curl2.vector, &resc);                                                             
            for(m = 0; m < 3; m++)                                                                            
            {                                                                                                 
                z[m] = points_x[m+3*i] - points_y[m+3*j];                                                     
            }                                                                                                 
            nrm = norm(z);                                                                                    
            resj += cexp((1*I)*kappa*nrm)/(nrm*4.0*PI)*(weights_x[i]*IE_x[i]*(resc-kappa*kappa*resn*phii[i]*phij[j]));
      
     
       } 
    
        res+= resj*weights_y[j]*IE_y[j];                                                                                                                                                          
    }                                     

    return gsl_complex_rect(creal(res), cimag(res));                                                          
} 

gsl_complex HP_kern3( double kappa, double *weights_x, double *weights_y, double *phii, double *phij, double *points_x, double *points_y, double *IE_x, double *IE_y, gsl_vector *n_x, gsl_vector *n_y, gsl_vector *curl_x, gsl_vector *curl_y, int sz_x, int sz_y )
{                                                                                                             
    int            i,j,m;                                                                                     
    double          z[3];                                                                                     
    double          nrm;                                                                                      
    double resn, resc;                                                                                        
    double complex res, resj;                                                                                 
                                                                                                              
    res = 0.0 + 0.0*I;                                                                                        
    gsl_blas_ddot(n_x, n_y, &resn);                                                                           
                                                                                                              
    for( j = 0; j < sz_y; j++ )                                                                               
    {                                                                                                         
        resj = 0.0 + 0.0*I;                                                                                   
        gsl_vector_view curl1 = gsl_matrix_row(curl_y,j);                                                     
        gsl_blas_ddot(&curl1.vector, curl_x, &resc);                                                 
        for( i = 0; i < sz_x; i++ )                                                                           
        {                                                                                                                                                    
            for(m = 0; m < 3; m++)                                                                            
            {                                                                                                 
                z[m] = points_x[m+3*i] - points_y[m+3*j];                                                     
            }                                                                                                 
            nrm = norm(z);                                                                                    
            resj += cexp((1*I)*kappa*nrm)/(nrm*4.0*PI)*(weights_x[i]*IE_x[i]*(resc-kappa*kappa*resn*phii[i]*phij[j]));
                                                                                                              
                                                                                                              
       }                                                                                                      
                                                                                                              
        res+= resj*weights_y[j]*IE_y[j];                                                                      
    }                                                                                                         
                                                                                                              
    return gsl_complex_rect(creal(res), cimag(res));                                                          
} 

gsl_complex HP_kern4( double kappa, double *weights_x, double *weights_y, double *phii, double *phij, double *points_x, double *points_y, double *IE_x, double *IE_y, gsl_vector *n_x, gsl_vector *n_y, gsl_matrix *curl_x, gsl_vector *curl_y, int sz_x, int sz_y )
{                                                                                                             
    int            i,j,m;                                                                                     
    double          z[3];                                                                                     
    double          nrm;                                                                                      
    double resn, resc;                                                                                        
    double complex res, resj;                                                                                 
                                                                                                              
    res = 0.0 + 0.0*I;                                                                                        
    gsl_blas_ddot(n_x, n_y, &resn);                                                                           
                                                                                                              
    for( j = 0; j < sz_y; j++ )                                                                               
    {                                                                                                         
        resj = 0.0 + 0.0*I;                                                                                   
        for( i = 0; i < sz_x; i++ )                                                                           
        {                                                                                                     
            gsl_vector_view curl2 = gsl_matrix_row(curl_x,i);                                                 
            gsl_blas_ddot(curl_y, &curl2.vector, &resc);                                               
            for(m = 0; m < 3; m++)                                                                            
            {                                                                                                 
                z[m] = points_x[m+3*i] - points_y[m+3*j];                                                     
            }                                                                                                 
            nrm = norm(z);                                                                                    
            resj += cexp((1*I)*kappa*nrm)/(nrm*4.0*PI)*(weights_x[i]*IE_x[i]*(resc-kappa*kappa*resn*phii[i]*phij[j]));
                                                                                                              
                                                                                                              
       }                                                                                                      
                                                                                                              
        res+= resj*weights_y[j]*IE_y[j];                                                                      
    }                                                                                                         
                                                                                                              
    return gsl_complex_rect(creal(res), cimag(res));                                                          
}   
void cross_product(const gsl_vector *u, const gsl_vector *v, gsl_vector *product)
{
    double p1 = gsl_vector_get(u, 1)*gsl_vector_get(v, 2)
    - gsl_vector_get(u, 2)*gsl_vector_get(v, 1);
    
    double p2 = gsl_vector_get(u, 2)*gsl_vector_get(v, 0)
    - gsl_vector_get(u, 0)*gsl_vector_get(v, 2);
    
    double p3 = gsl_vector_get(u, 0)*gsl_vector_get(v, 1)
    - gsl_vector_get(u, 1)*gsl_vector_get(v, 0);
    
    gsl_vector_set(product, 0, p1);
    gsl_vector_set(product, 1, p2);
    gsl_vector_set(product, 2, p3);
}

void cross_product2( double *u, gsl_vector *v, int sz )
{
    int i;
    double aux[3];
    for(i=0;i<sz;i++)
    {
        aux[0] = u[3*i+1]*gsl_vector_get(v, 2)-u[3*i+2]*gsl_vector_get(v, 1);
        aux[1] = u[3*i+2]*gsl_vector_get(v, 0)-u[3*i]*gsl_vector_get(v, 2);
        aux[2] = u[3*i]*gsl_vector_get(v, 1)-u[3*i+1]*gsl_vector_get(v, 0);
        u[3*i] = aux[0];
        u[3*i+1] = aux[1];
        u[3*i+2] = aux[2];
    }
}

gsl_vector *get_normal(gsl_vector *p1, gsl_vector *p2, gsl_vector *p3)
{
    gsl_vector *v1 = gsl_vector_alloc(3);
    gsl_vector *v2 = gsl_vector_alloc(3);
    gsl_vector *normal = gsl_vector_alloc(3);
    double nrm;
    gsl_blas_dcopy(p2, v1);
    gsl_blas_dcopy(p3, v2);
    gsl_blas_daxpy(-1.0,p1,v1);
    gsl_blas_daxpy(-1.0,p1,v2);
    cross_product(v1,v2,normal);
    nrm = gsl_blas_dnrm2 (normal);
    gsl_blas_dscal (1.0/nrm, normal);
    return normal;
}

gsl_vector *normal(double *points, double *v)
{
    int i;
    gsl_vector *point1 = gsl_vector_complex_alloc(3);
    gsl_vector *point2 = gsl_vector_complex_alloc(3);
    gsl_vector *point3 = gsl_vector_complex_alloc(3);
    gsl_vector *n = gsl_vector_alloc(3);
    gsl_vector_set(point1,0,v[0]);
    gsl_vector_set(point1,1,v[1]);
    gsl_vector_set(point1,2,v[2]);
    gsl_vector_set(point2,0,v[3]);
    gsl_vector_set(point2,1,v[4]);
    gsl_vector_set(point2,2,v[5]);
    gsl_vector_set(point3,0,v[6]); 
    gsl_vector_set(point3,1,v[7]); 
    gsl_vector_set(point3,2,v[8]); 
    //gsl_vector_set(point3,0,points[0]);
    //gsl_vector_set(point3,1,points[1]);
    //gsl_vector_set(point3,2,points[2]);
    gsl_vector_memcpy(n,get_normal(point3, point1, point2));

    return n;
}


gsl_matrix *invert_a_matrix(gsl_matrix *matrix)
{
    gsl_permutation *p = gsl_permutation_alloc(matrix->size1);
    int s;

    // Compute the LU decomposition of this matrix
    gsl_linalg_LU_decomp(matrix, p, &s);

    // Compute the  inverse of the LU decomposition
    gsl_matrix *inv = gsl_matrix_alloc(matrix->size1, matrix->size1);
    gsl_linalg_LU_invert(matrix, p, inv);

    gsl_permutation_free(p);

    return inv;
}

void divergence(double a,double len, double *v,char *tp, int opt, double *div, int sz)
{
    int i;

    if(strcmp(tp,"RWG")==0)
    {
        if(opt==0)
        {
            for(i=0;i<sz;i++)
            {    
                div[i] = len/a;
            }
        }
        else
        {
            for(i=0;i<sz;i++)
            {
                div[i] = -len/a;
            }
        }
    }
}

gsl_vector *curl_bf(double *v, int dof, gsl_vector *n)
{
    int i;
    gsl_vector *grad = gsl_vector_alloc(2);
    gsl_vector *curl = gsl_vector_alloc(3);
    gsl_vector *aux = gsl_vector_alloc(2);
    gsl_vector *aux2 = gsl_vector_alloc(3);
    gsl_matrix *J = gsl_matrix_alloc(3, 2);
    gsl_matrix *G = gsl_matrix_alloc(2, 2);
    gsl_matrix *Inv = gsl_matrix_alloc(2, 2);
    double det;
    if(dof==0)
    {
        gsl_vector_set(grad,0,-1);
        gsl_vector_set(grad,1,-1);
    }
    else if(dof==1)
    {
        gsl_vector_set(grad,0, 1);
        gsl_vector_set(grad,1, 0);
    }
    else if(dof==2)
    {
        gsl_vector_set(grad,0, 0);
        gsl_vector_set(grad,1, 1);
    }
    gsl_matrix_set(J,0,0, v[3]-v[0]);
    gsl_matrix_set(J,0,1, v[6]-v[0]);
    gsl_matrix_set(J,1,0, v[4]-v[1]);
    gsl_matrix_set(J,1,1, v[7]-v[1]);
    gsl_matrix_set(J,2,0, v[5]-v[2]);
    gsl_matrix_set(J,2,1, v[8]-v[2]);
    
    gsl_blas_dgemm(CblasTrans,CblasNoTrans, 1.0, J, J, 0.0, G);
    det = fabs(gsl_matrix_get(G,0,0)*gsl_matrix_get(G,1,1)- gsl_matrix_get(G,0,1)*gsl_matrix_get(G,1,0));
    gsl_matrix_set(Inv,0,0,gsl_matrix_get(G,1,1)/det);
    gsl_matrix_set(Inv,0,1,-gsl_matrix_get(G,1,0)/det);
    gsl_matrix_set(Inv,1,0,-gsl_matrix_get(G,0,1)/det);
    gsl_matrix_set(Inv,1,1,gsl_matrix_get(G,0,0)/det);
    gsl_blas_dgemv(CblasNoTrans,1,Inv,grad,0,aux);
    gsl_blas_dgemv(CblasNoTrans,1,J,aux,0,aux2);
    cross_product(n,aux2,curl);
    gsl_vector_free(grad);
    gsl_matrix_free(Inv);
    gsl_vector_free(aux);
    gsl_vector_free(aux2);
    gsl_matrix_free(J);
    gsl_matrix_free(G);

    return curl;
}


gsl_vector *curl_bf2(double *v, int dof, gsl_vector *n, double eta, double chi)                                                        
{                                                                                                             
    int i;                                                                                                    
    gsl_vector *grad = gsl_vector_alloc(2);                                                                   
    gsl_vector *curl = gsl_vector_alloc(3);                                                                   
    gsl_vector *aux = gsl_vector_alloc(2);                                                                    
    gsl_vector *aux2 = gsl_vector_alloc(3);                                                                   
    gsl_matrix *J = gsl_matrix_alloc(3, 2);                                                                   
    gsl_matrix *G = gsl_matrix_alloc(2, 2);                                                                   
    gsl_matrix *Inv = gsl_matrix_alloc(2, 2);                                                                 
    double det;                                                                                               
    if(eta>chi)                                                                                                
    {                                                                                                         
        gsl_vector_set(grad,0,2);                                                                            
        gsl_vector_set(grad,1,0);                                                                            
    }                          
    else
    {
        gsl_vector_set(grad,0,0);                                                                             
        gsl_vector_set(grad,1,2);
    }                                                                               
    gsl_matrix_set(J,0,0, v[3]-v[0]);                                                                         
    gsl_matrix_set(J,0,1, v[6]-v[0]);                                                                         
    gsl_matrix_set(J,1,0, v[4]-v[1]);                                                                         
    gsl_matrix_set(J,1,1, v[7]-v[1]);                                                                         
    gsl_matrix_set(J,2,0, v[5]-v[2]);                                                                         
    gsl_matrix_set(J,2,1, v[8]-v[2]);                                                                         
                                                                                                             
    gsl_blas_dgemm(CblasTrans,CblasNoTrans, 1.0, J, J, 0.0, G);                                               
    det = fabs(gsl_matrix_get(G,0,0)*gsl_matrix_get(G,1,1)- gsl_matrix_get(G,0,1)*gsl_matrix_get(G,1,0));     
   // printf("det %f %f %f %f %f %f %f %f %f %f\n",v[0],v[1],v[2],v[3],v[4],v[5],v[6],v[7],v[8],det);
    gsl_matrix_set(Inv,0,0,gsl_matrix_get(G,1,1)/det);                                                        
    gsl_matrix_set(Inv,0,1,-gsl_matrix_get(G,1,0)/det);                                                       
    gsl_matrix_set(Inv,1,0,-gsl_matrix_get(G,0,1)/det);                                                       
    gsl_matrix_set(Inv,1,1,gsl_matrix_get(G,0,0)/det);                                                        
    gsl_blas_dgemv(CblasNoTrans,1,Inv,grad,0,aux);                                                            
    gsl_blas_dgemv(CblasNoTrans,1,J,aux,0,aux2);                                                              
    cross_product(n,aux2,curl);                                                                               
    gsl_vector_free(grad);                                                                                    
    gsl_matrix_free(Inv);                                                                                     
    gsl_vector_free(aux);                                                                                     
    gsl_vector_free(aux2);                                                                                    
    gsl_matrix_free(J);                                                                                       
    gsl_matrix_free(G);                                                                                       
                                                                                                              
    return curl;                                                                                              
}
gsl_matrix *curl_bfq(double *v, double *eta, double *chi, int sz, int dof, gsl_vector *n)
{
    int i, j, id;
    gsl_vector *grad = gsl_vector_alloc(2);
    gsl_matrix *curl = gsl_matrix_alloc(sz*sz,3);
    gsl_vector *aux = gsl_vector_alloc(2);
    gsl_vector *aux2 = gsl_vector_alloc(3);
    gsl_matrix *J = gsl_matrix_alloc(3, 2);
    gsl_matrix *G = gsl_matrix_alloc(2, 2);
    gsl_matrix *Inv = gsl_matrix_alloc(2, 2);
    double vert1[9];
    double vert2[9];
    double det;
    id = 0;
   for(i = 0; i<sz; i++)
    {

        for(j = 0; j<sz; j++)
        {
            gsl_vector_view row = gsl_matrix_row(curl,id);
            
            if((chi[i]<0 && eta[j]==chi[i])||eta[j]<chi[i])                                               
            {                                                                                             
                gsl_vector_set(grad,0,-0.5);                                                                   
                gsl_vector_set(grad,1,0.0);                                                                 
            }                                                                                             
            else if((chi[i]>0 && eta[j]==chi[i])||eta[j]>chi[i])                                          
            { 
                gsl_vector_set(grad,0,0.0);
                gsl_vector_set(grad,1,-0.5); 
            }
            //gsl_vector_set(grad,0,0.25*(eta[j]-1));
            gsl_matrix_set(J,0,0, deriv_q( v[0], v[3], v[6], v[9], chi[i], eta[j], 1 ));
            gsl_matrix_set(J,0,1, deriv_q( v[0], v[3], v[6], v[9], chi[i], eta[j], 0 ));
            gsl_matrix_set(J,1,0, deriv_q( v[1], v[4], v[7], v[10], chi[i], eta[j], 1 ));
            gsl_matrix_set(J,1,1, deriv_q( v[1], v[4], v[7], v[10], chi[i], eta[j], 0 ));
            gsl_matrix_set(J,2,0, deriv_q( v[2], v[5], v[8], v[11], chi[i], eta[j], 1 ));
            gsl_matrix_set(J,2,1, deriv_q( v[2], v[5], v[8], v[11], chi[i], eta[j], 0 ));   
            gsl_blas_dgemm(CblasTrans,CblasNoTrans, 1.0, J, J, 0.0, G);
            det = fabs(gsl_matrix_get(G,0,0)*gsl_matrix_get(G,1,1)- gsl_matrix_get(G,0,1)*gsl_matrix_get(G,1,0));
            gsl_matrix_set(Inv,0,0,gsl_matrix_get(G,1,1)/det);
            gsl_matrix_set(Inv,0,1,-gsl_matrix_get(G,1,0)/det);
            gsl_matrix_set(Inv,1,0,-gsl_matrix_get(G,0,1)/det);
            gsl_matrix_set(Inv,1,1,gsl_matrix_get(G,0,0)/det);
            gsl_blas_dgemv(CblasNoTrans,1,Inv,grad,0,aux);
            gsl_blas_dgemv(CblasNoTrans,1,J,aux,0,aux2);
            cross_product(n,aux2,&row.vector);
           // printf("curl: %f %f %f \n",gsl_vector_get(&row.vector,0),gsl_vector_get(&row.vector,1), gsl_vector_get(&row.vector,2));
            id++;
            
        }
    }
    gsl_vector_free(grad);
    gsl_matrix_free(Inv);
    gsl_vector_free(aux);
    gsl_vector_free(aux2);
    gsl_matrix_free(J);
    gsl_matrix_free(G);
    
    return curl;
    
}

gsl_matrix *curl_bfq2(double *v, double *eta, double *chi, int sz, int dof, gsl_vector *n)                     
{                                                                                                             
    int i, j, id;                                                                                             
    gsl_vector *grad = gsl_vector_alloc(2);                                                                   
    gsl_matrix *curl = gsl_matrix_alloc(sz*sz,3);                                                             
    gsl_vector *aux = gsl_vector_alloc(2);                                                                    
    gsl_vector *aux2 = gsl_vector_alloc(3);                                                                   
    gsl_matrix *J = gsl_matrix_alloc(3, 2);                                                                   
    gsl_matrix *G = gsl_matrix_alloc(2, 2);                                                                   
    gsl_matrix *Inv = gsl_matrix_alloc(2, 2);                                                                 
    double vert1[9];                                                                                          
    double vert2[9];                                                                                          
    double det;                                                                                               
    id = 0;                                                                                                   
   for(i = 0; i<sz; i++)                                                                                      
    {                                                                                                         
                                                                                                              
        for(j = 0; j<sz; j++)                                                                                 
        {                                                                                                     
            gsl_vector_view row = gsl_matrix_row(curl,id);                                                    
                                                                                                              
            if((chi[i]<0 && eta[j]==chi[i])||eta[j]<chi[i])                                                   
            {                                                                                                 
                gsl_vector_set(grad,0,0.5);                                                                  
                gsl_vector_set(grad,1,0.0);                                                                   
            }                                                                                                 
            else if((chi[i]>0 && eta[j]==chi[i])||eta[j]>chi[i])                                              
            {                                                                                                 
                gsl_vector_set(grad,0,0.0);                                                                   
                gsl_vector_set(grad,1,0.5);                                                                  
            }                                                                                                 
            //gsl_vector_set(grad,0,0.25*(eta[j]-1));                                                           
            gsl_matrix_set(J,0,0, deriv_q( v[0], v[3], v[6], v[9], chi[i], eta[j], 1 ));                      
            gsl_matrix_set(J,0,1, deriv_q( v[0], v[3], v[6], v[9], chi[i], eta[j], 0 ));                      
            gsl_matrix_set(J,1,0, deriv_q( v[1], v[4], v[7], v[10], chi[i], eta[j], 1 ));                     
            gsl_matrix_set(J,1,1, deriv_q( v[1], v[4], v[7], v[10], chi[i], eta[j], 0 ));                     
            gsl_matrix_set(J,2,0, deriv_q( v[2], v[5], v[8], v[11], chi[i], eta[j], 1 ));                     
            gsl_matrix_set(J,2,1, deriv_q( v[2], v[5], v[8], v[11], chi[i], eta[j], 0 ));                     
            gsl_blas_dgemm(CblasTrans,CblasNoTrans, 1.0, J, J, 0.0, G);                                       
            det = fabs(gsl_matrix_get(G,0,0)*gsl_matrix_get(G,1,1)- gsl_matrix_get(G,0,1)*gsl_matrix_get(G,1,0));
            gsl_matrix_set(Inv,0,0,gsl_matrix_get(G,1,1)/det);                                                
            gsl_matrix_set(Inv,0,1,-gsl_matrix_get(G,1,0)/det);                                               
            gsl_matrix_set(Inv,1,0,-gsl_matrix_get(G,0,1)/det);                                               
            gsl_matrix_set(Inv,1,1,gsl_matrix_get(G,0,0)/det);                                                
            gsl_blas_dgemv(CblasNoTrans,1,Inv,grad,0,aux);                                                    
            gsl_blas_dgemv(CblasNoTrans,1,J,aux,0,aux2);                                                      
            cross_product(n,aux2,&row.vector);                                                                
           // printf("curl: %f %f %f \n",gsl_vector_get(&row.vector,0),gsl_vector_get(&row.vector,1), gsl_vector_get(&row.vector,2));
            id++;                                                                                             
                                                                                                              
        }                                                                                                     
    }                                                                                                         
    gsl_vector_free(grad);                                                                                    
    gsl_matrix_free(Inv);                                                                                     
    gsl_vector_free(aux);                                                                                     
    gsl_vector_free(aux2);                                                                                    
    gsl_matrix_free(J);                                                                                       
    gsl_matrix_free(G);                                                                                       
                                                                                                              
    return curl;
}
void get_ie_t( int sz_eta, double *M, double *IE)
{

    int i;
    double A = area(M);

    for(i = 0; i < sz_eta; i++)
    {
        IE[i] = A;
    }
}

void get_ie_q( int sz_eta,int sz_chi,double *M,double *eta, double *chi, double *IE )
{
    double dxeta, dxchi, dyeta, dychi, dzeta, dzchi, point[3];
    int  index=0, i, j;

    for( i = 0; i < sz_eta; i++ )
    {
        for( j=0; j < sz_chi; j++ )
        {
            dxeta = deriv_q( M[0], M[3], M[6], M[9], chi[i], eta[j], 1 );
            dxchi = deriv_q( M[0], M[3], M[6], M[9], chi[i], eta[j], 0 );
            dyeta = deriv_q( M[1], M[4], M[7], M[10], chi[i], eta[j], 1 );
            dychi = deriv_q( M[1], M[4], M[7], M[10], chi[i], eta[j], 0 );
            dzeta = deriv_q( M[2], M[5], M[8], M[11], chi[i], eta[j], 1 );
            dzchi = deriv_q( M[2], M[5], M[8], M[11], chi[i], eta[j], 0 );
            point[0]=dychi*dzeta-dyeta*dzchi;
            point[1]=-(dxchi*dzeta-dxeta*dzchi);
            point[2]=dxchi*dyeta-dxeta*dychi;
            IE[index] = norm(point);
            index++;
        }
    }
}

void points_t( double *M, int sz_eta, double *eta, double *chi, double *points_aux)
{
    int i;
    
    for ( i = 0; i<sz_eta; i++ )
    {
        points_aux[3*i]   = nodal_function_t( M[0], M[3], M[6], eta[i], chi[i] );
        points_aux[3*i+1] = nodal_function_t( M[1], M[4], M[7], eta[i], chi[i] );
        points_aux[3*i+2] = nodal_function_t( M[2], M[5], M[8], eta[i], chi[i] );
        
    }
}

void points_q( double *M, int sz_eta, int sz_chi, double *eta, double *chi, double *points_aux)
{
    int index = 0, i, j;
    
    for ( i = 0; i<sz_eta; i++ )
    {
        for( j = 0; j<sz_chi; j++ )
        {
            points_aux[3*index]   = nodal_function_q( M[0], M[3], M[6], M[9], eta[i], chi[j] );
            points_aux[3*index+1] = nodal_function_q( M[1], M[4], M[7], M[10], eta[i], chi[j] );
            points_aux[3*index+2] = nodal_function_q( M[2], M[5], M[8], M[11], eta[i], chi[j] );
            index++;
        }
    }
}

double nodal_function_t( double x1, double x2, double x3, double eta, double chi )
{
    double N1 = 1-eta-chi;
    double N2 = chi;
    double N3 = eta;
    
    return x1*N1+x2*N2+x3*N3;
}

double nodal_function_q( double x1, double x2, double x3, double x4, double eta, double chi )
{
    double N1 = 0.25*(1-chi)*(1-eta);
    double N2 = 0.25*(1+chi)*(1-eta);
    double N3 = 0.25*(1+chi)*(1+eta);
    double N4 = 0.25*(1-chi)*(1+eta);
    
    return x1*N1+x2*N2+x3*N3+x4*N4;
}


void quad_rules_PRIMAL( double *eta, double *chi, int *sz, int order )
{
    int i;
    double v1[3] = {1/6.,2./3, 1./6};
    double v2[3] = {1./6,1./6, 2./3};
    double v3[4] = {1./3, 1./5, 1./5, 3./5};
    double v4[4] = {1./3, 1./5, 3./5, 1./5};
    
    for( i = 0; i < 4; i++ )
    {
        eta[i]=1.0/3.0;
        chi[i]=1.0/3.0;
    }
    
    if( order == 1 )
    {
        sz[0] = 1;
        sz[1] = 1;
        eta[0] = 1.0/3;
        chi[0] = 1.0/3;
    }
    else if( order == 2 )
    {
        sz[0] = 3;
        sz[1] = 3;
        
        for( i = 0; i < sz[0] ; i++ )
        {
            eta[i] = v2[i];
            chi[i] = v1[i];
        }
    }
    else
    {
        sz[0] = 4;
        sz[1] = 4;
        
        for( i = 0; i < sz[0]; i++  )
        {
            eta[i] = v4[i];
            chi[i] = v3[i];
        }
    }
}

void quad_rules_DUAL( double *eta, double *chi, int *sz, int order )
{
    int i;
    double v1[2] = {-1./sqrt(3),1./sqrt(3)};
    double v2[3] = {-sqrt(3./5),0,sqrt(3./5)};
    double v3[4] = {-sqrt((3/7)-(2/7)*sqrt(6/5)),sqrt((3/7)-(2/7)*sqrt(6/5)),-sqrt((3/7)+(2/7)*sqrt(6/5)),sqrt((3/7)+(2/7)*sqrt(6/5))};
    double v4[5] = { 0, 0.538469,-0.538469,0.90618,-0.90618};
    for( i = 0; i < 4; i++ )
    {
        eta[i]=0.0;
        chi[i]=0.0;
    }
    
    if( order == 1 )
    {
        sz[0] = 1;
        sz[1] = 1;
        eta[0] = 0.0;
        chi[0] = 0.0;
    }
    else if( order == 2 )
    {
        sz[0] = 2;
        sz[1] = 2;
        
        for( i = 0; i < sz[0] ; i++ )
        {
            eta[i] = v1[i];
            chi[i] = v1[i];
        }
    }
    else if ( order == 3 )
    {
        sz[0] = 3;
        sz[1] = 3;
        
        for( i = 0; i < sz[0]; i++  )
        {
            eta[i] = v2[i];
            chi[i] = v2[i];
        }
    }
    else if ( order == 4 )
    {
        sz[0] = 4;
        sz[1] = 4;
        
        for( i = 0; i < sz[0]; i++ )
        {
            eta[i] = v3[i];
            chi[i] = v3[i];
        }
    }
    else if ( order == 5 )                                                                                                   
    {                                                                                                         
        sz[0] = 5;                                                                                            
        sz[1] = 5;                                                                                            
                                                                                                              
        for( i = 0; i < sz[0]; i++ )                                                                          
        {                                                                                                     
            eta[i] = v4[i];                                                                                   
            chi[i] = v4[i];                                                                                   
        }                                                                                                     
    }
}

void weights_int( double *weights, int order )
{
    if( order == 1 )
    {
        weights[0] = 2;
    }
    else if( order == 2 )
    {
        weights[0] = 1.;
        weights[1] = 1.;
    }
    else if( order == 3 )
    {
        weights[0] = 5./9;
        weights[1] = 8./9;
        weights[2] = 5./9;
    }
    else if( order == 4 )
    {
        weights[0] = (18+sqrt(30))/36;
        weights[1] = (18+sqrt(30))/36;
        weights[2] = (18-sqrt(30))/36;
        weights[3] = (18-sqrt(30))/36;
    }
    else if( order == 5 )                                                                                      
    {
        weights[0] = 0.568889;
        weights[1] = 0.478629;
        weights[2] = 0.478629;
        weights[3] = 0.236927;
        weights[4] = 0.236927;
    }
}

void quad_rules_ab( double a, double b, double *eta, double *chi, double *weights1, double *weights2, int *sz, int order)
{
    
    int i;
    quad_rules_DUAL( eta, chi, sz, order );
    weights_int( weights1, order );
    weights_int( weights2, order );
    for(i=0; i< sz[0]; i++)
    {
        eta[i] = (b-a)*0.5*eta[i] + (a+b)*0.5;
        chi[i] = (b-a)*0.5*chi[i] + (a+b)*0.5;
        weights1[i] = (b-a)*0.5*weights1[i];
        weights2[i] = (b-a)*0.5*weights2[i];
    }
}


void weights_DUAL( double *weights, int order, int sz )
{
    double v0[1] = {2.};
    double v1[2] = {1., 1.};
    double v2[3] = {5./9, 8./9, 5./9};
    double v3[4] = {(18+sqrt(30))/36,(18+sqrt(30))/36,(18-sqrt(30))/36, (18-sqrt(30))/36};
    double v4[5] = {0.568889,0.478629,0.478629,0.236927,0.236927};
    int i, j, index=0;
    
    if( order==1 )
    {
        for( i = 0; i < sz; i++ )
        {
            for( j = 0;j <sz; j++ )
            {
                weights[index] = v0[i]*v0[j];
                index++;
            }
            
        }
    }
    else if ( order == 2 )
    {
        for( i = 0; i < sz; i++ )
        {
            for( j = 0;j <sz; j++ )
            {
                weights[index] = v1[i]*v1[j];
                index++;
            }
            
        }
    }
    else if ( order == 3 )
    {
        for( i = 0; i < sz; i++ )
        {
            for( j = 0;j <sz; j++ )
            {
                weights[index] = v2[i]*v2[j];
                index++;
            }
            
        }
    }
    else if ( order == 4 )
    {
        for( i = 0; i < sz; i++ )
        {
            for( j = 0;j <sz; j++ )
            {
                weights[index] = v3[i]*v3[j];
                index++;
            }
            
        }
    }
    else if ( order == 5 )                                                                                    
    {                                                                                                         
        for( i = 0; i < sz; i++ )                                                                             
        {                                                                                                     
            for( j = 0;j <sz; j++ )                                                                           
            {                                                                                                 
                weights[index] = v4[i]*v4[j];                                                                 
                index++;                                                                                      
            }                                                                                                 
                                                                                                              
        }                                                                                                     
    } 
}

void weights_PRIMAL( double *weights, int order, int sz )
{
    double v0[1] = {1./2};
    double v1[3] = {1./6, 1./6, 1./6};
    double v2[4] = {-27./96, 25./96, 25./96, 25./96};
    
    int i;
    
    if( order==1 )
    {
        for( i = 0; i < sz; i++ )
        {
            weights[i] = v0[i];
        }
    }
    else if ( order == 2 )
    {
        for( i = 0; i < sz; i++ )
        {
            weights[i] = v1[i];
        }
    }
    else if ( order == 3 )
    {
        for( i = 0; i < sz; i++ )
        {
            weights[i] = v2[i];
        }
    }
}

void ev_basis_function( int sz, int dof, double *phi, double *eta, double *chi, char *tp )
{
    int i,j,id, msz;
    
    if( sz % 2 == 0 )
    {
        msz = sz/2;
    }
    else
    {
        msz = (sz+1)/2;
    }
    
    if( strcmp(tp,"P0") == 0 || strcmp(tp,"P0d") == 0 )
    {
    
        for( i = 0; i < sz; i++ )
        {
            phi[i] = 1.0;
        }
        
    }
    else if( strcmp(tp,"P1di") == 0 || strcmp(tp,"P1bi") == 0)                                                  
    {
        
        id = 0;
        /*for(i = 0; i < msz; i++)                                                                              
        {                                                                                                     
            for( j = 0; j< msz; j++)                                                                          
            {        
                phi[id] = 0.25*(1-chi[i])*(1-eta[j]); 
                //phi[id] = 0.25*(1+chi[i])*(1+eta[j]);
                id++;
            }                                                                                                 
        } */
        for(i = 0; i < msz; i++)                                                                              
        {                                                                                                     
            for( j = 0; j< msz; j++)                                                                          
            {
                if((chi[i]<0 && eta[j]==chi[i])||eta[j]<chi[i])                                                   
                {                         
                    phi[id] = 0.5*(1-chi[i]);                                                                                                  
                }                                                                                                 
                else if((chi[i]>0 && eta[j]==chi[i])||eta[j]>chi[i])                                              
                {                                                                                                 
                    phi[id] = 0.5*(1-eta[j]);                          
                }
                id++;
            }
        }
                                                                                                              
    }  
    else if( strcmp(tp,"P1dii") == 0 )                                                
    {
        id = 0;
        for(i = 0; i < msz; i++)                                                                              
        {                                                                                                     
            for( j = 0; j< msz; j++)                                                                          
            {                                                                                                 
                if((chi[i]<0 && eta[j]==chi[i])||eta[j]<chi[i])                                               
                {                                                                                             
                    phi[id] = 0.5*(chi[i]);                                                                 
                }                                                                                             
                else if((chi[i]>0 && eta[j]==chi[i])||eta[j]>chi[i])                                          
                {                                                                                             
                    phi[id] = 0.5*(eta[j]);                                                                 
                }                                                                                             
                id++;                                                                                         
            }                                                                                                 
        }                                                                                                     
                                                                                                              
    }
    else if( strcmp(tp,"P1") == 0 || strcmp(tp,"P1b") == 0 || strcmp(tp,"P1d") == 0 )
    {
    
        if(dof==0)
        {
            
            for(i = 0; i < sz; i++)
            {
                phi[i] = 1-chi[i]-eta[i];
            }
        }
        else if(dof==1)
        {
        
            for(i = 0; i < sz; i++)
            {
                phi[i] = chi[i];
            }
        }
        else if(dof==2)
        {
            
            for(i = 0; i < sz; i++)
            {
                phi[i] = eta[i];
            }
        }
        
    }
    else if( strcmp(tp,"P1ii") == 0  )                         
    {
        for(i = 0; i < sz; i++)                                                                           
        {
            if(eta[i]>chi[i])
            {
                phi[i]=chi[i]/0.5;
            }
            else
            {
                phi[i]=eta[i]/0.5; 
            }
        }

    }
    
}

void ev_basis_function_maxwell(double *points, double *phi, double *p0, double A, double len, gsl_vector *n, int sz, int opt, char *tp)
{
    int i;

    if(opt==0)
    {
        for(i=0;i<sz;i++)
        {

            phi[3*i] = (points[3*i]-p0[0])*len/(2*A);
            phi[3*i+1] = (points[3*i+1]-p0[1])*len/(2*A); 
            phi[3*i+2] = (points[3*i+2]-p0[2])*len/(2*A); 
        }
    }
    else
    {
        for(i=0;i<sz;i++)                                                                                     
        {                                                                                                     
            phi[3*i] = -(points[3*i]-p0[0])*len/(2*A);                                                         
            phi[3*i+1] = -(points[3*i+1]-p0[1])*len/(2*A);                                                   
            phi[3*i+2] = -(points[3*i+2]-p0[2])*len/(2*A); 
        } 
    }
    if( strcmp(tp,"NC") == 0 )
    {
        cross_product2(phi,n,sz);
    }
}

void evaluate_sing_kern(double *y,  double *M, double *V, int order, double *res, int dof, double kappa, char *tp)
{
    int       i, j, index;
    double       point[3];
    double           z[3];
    double    weights[16];
    double eta[4], chi[4];
    double chi2[1],eta2[1];
    double       nrm,nrm2;
    double  complex res_j;
    double              J;
    double           T[9];
    double jac[3], phi[1];
    int             sz[2];
    double  points_aux[3];
    
    jac[0] = (M[1]-y[1])*(M[5]-M[2]) - (M[2]-y[2])*(M[4]-M[1]);
    jac[1] = (M[0]-y[0])*(M[5]-M[2]) - (M[2]-y[2])*(M[3]-M[0]);
    jac[2] = (M[0]-y[0])*(M[4]-M[1]) - (M[1]-y[1])*(M[3]-M[0]);
    J = norm(jac);
    quad_rules_DUAL( eta, chi, sz, order );
    weights_DUAL( weights, order, sz[0] );
    
    index = 0;
    phi[0] = 1.0;
    res[0] = 0.0;
    res[1] = 0.0;

    for(i = 0; i<sz[0]; i++)
    {
        eta[i] = (eta[i]+1)*0.5;
        chi[i] = (chi[i]+1)*0.5;
    }
    for( i = 0; i < sz[0]; i++ )
    {
        for( j=0; j < sz[1]; j++ )
        { 
            z[0] = -y[0]+(1-eta[j])*M[0]+eta[j]*M[3];
            z[1] = -y[1]+(1-eta[j])*M[1]+eta[j]*M[4];
            z[2] = -y[2]+(1-eta[j])*M[2]+eta[j]*M[5];
            
            point[0] = chi[i]*z[0];
            point[1] = chi[i]*z[1];
            point[2] = chi[i]*z[2];
            
            points_aux[0] = point[0]+y[0];
            points_aux[1] = point[1]+y[1];
            points_aux[2] = point[2]+y[2];
            
            nrm = norm(z);
            nrm2 = norm(point);
            
            if(strcmp(tp,"P1")==0  || strcmp(tp,"P1b")==0 )
            {
                transpose( 3, 3, V, T );
                if(T[0] ==0.0 && T[1] == 0.0 && T[2] ==0.0)                                                               
                {                                                                                                         
                    double T2[4], points_aux2[2];                                                           
                    T2[0] = T[4]-T[3];                                                                                    
                    T2[1] = T[5]-T[3];                                                                                    
                    T2[2] = T[7]-T[6];                                                                                    
                    T2[3] = T[8]-T[6];                                                                                    
                    inverse(2, T2);                                                                                                                                                                                                                 
                    points_aux2[0] = points_aux[1]-T[3];                                                        
                    points_aux2[1] = points_aux[2]-T[6];                                                                                                             
                    matrix_prod( 2, 2, 2,1, T2, points_aux2 );                                               
                    eta2[0] = points_aux2[0];                                                                      
                    chi2[0] = points_aux2[1];                                                                                          
                } 
                else if(T[3] ==0.0 && T[4] == 0.0 && T[5] ==0.0)                                                   
                {                                                                                             
                    double T2[4], points_aux2[2];                                               
                    T2[0] = T[1]-T[0];                                                                                  
                    T2[1] = T[2]-T[0];                                                                        
                    T2[2] = T[7]-T[6];                                                                        
                    T2[3] = T[8]-T[6];                                                                        
                    inverse(2, T2);                                                                           
                    points_aux2[0] = points_aux[0]-T[0];                                            
                    points_aux2[1] = points_aux[2]-T[6];                                                                                                                                   
                    matrix_prod( 2, 2, 2,1, T2, points_aux2 );                                                  
                    eta2[0] = points_aux2[0];                                                                 
                    chi2[0] = points_aux2[1];                                                                 
                }
                else if(T[6] ==0.0 && T[7] == 0.0 && T[8] ==0.0)                                              
                {                                                                                             
                    double T2[4], points_aux2[2];                                               
                    T2[0] = T[1]-T[0];                                                                        
                    T2[1] = T[2]-T[0];                                                                        
                    T2[2] = T[4]-T[3];                                                                        
                    T2[3] = T[5]-T[3];                                                                        
                    inverse(2, T2);                                                                           
                    points_aux2[0] = points_aux[0]-T[0];                                                    
                    points_aux2[1] = points_aux[1]-T[3];
                    matrix_prod( 2, 2, 2, 1, T2, points_aux2 );                                                  
                    eta2[0] = points_aux2[0];                                                                 
                    chi2[0] = points_aux2[1];                                                                 
                }
                else
                {
                    inverse(3, T);
                    matrix_prod( 3, 3, 3, 1, T, points_aux );
                    eta2[0] = points_aux[1];
                    chi2[0] = points_aux[2];
                }
                ev_basis_function( 1, dof, phi, eta2, chi2, tp );               
            }
            res_j =0.25*phi[0]*cexp((1*I)*kappa*nrm2)*J*weights[index]/(nrm*4.0*PI);
            res[0] += creal(res_j);
            res[1] += cimag(res_j);
            index++;
        }
    }

}

void evaluate_sing_kern_maxwell_weakly(double *y,  double *M, double *p0, double A, double len, gsl_vector *n, int opt, int order, double *res, double kappa)
{                                                                                                             
    int       i, j, index;                                                                                    
    double       point[3];                                                                                    
    double           z[3];                                                                                    
    double    weights[16];                                                                                    
    double eta[4], chi[4];                                                                                    
    double chi2[1],eta2[1];                                                                                   
    double       nrm,nrm2;                                                                                    
    double  complex res_j[3];                                                                                    
    double              J;                                                                                                                                         
    double jac[3], phi[3];                                                                                    
    int             sz[2];                                                                                    
    double  points_aux[3];                                                                                    
                                                                                                              
    jac[0] = (M[1]-y[1])*(M[5]-M[2]) - (M[2]-y[2])*(M[4]-M[1]);                                               
    jac[1] = (M[0]-y[0])*(M[5]-M[2]) - (M[2]-y[2])*(M[3]-M[0]);                                               
    jac[2] = (M[0]-y[0])*(M[4]-M[1]) - (M[1]-y[1])*(M[3]-M[0]);                                               
    J = norm(jac);                                                                                            
    quad_rules_DUAL( eta, chi, sz, order );                                                                   
    weights_DUAL( weights, order, sz[0] );                                                                    
                                                                                                              
    index = 0;                                                                                                
    phi[0] = 1.0;                                                                                             
    res[0] = 0.0;                                                                                             
    res[1] = 0.0;                                                                                             
                                                                                                              
    for(i = 0; i<sz[0]; i++)                                                                                  
    {                                                                                                         
        eta[i] = (eta[i]+1)*0.5;                                                                              
        chi[i] = (chi[i]+1)*0.5;                                                                              
    }                                                                                                         
    for( i = 0; i < sz[0]; i++ )                                                                              
    {                                                                                                         
        for( j=0; j < sz[1]; j++ )                                                                            
        {                                                                                                     
            z[0] = -y[0]+(1-eta[j])*M[0]+eta[j]*M[3];                                                         
            z[1] = -y[1]+(1-eta[j])*M[1]+eta[j]*M[4];                                                         
            z[2] = -y[2]+(1-eta[j])*M[2]+eta[j]*M[5];                                                         
                                                                                                              
            point[0] = chi[i]*z[0];                                                                           
            point[1] = chi[i]*z[1];                                                                           
            point[2] = chi[i]*z[2];                                                                           
                                                                                                              
            points_aux[0] = point[0]+y[0];                                                                    
            points_aux[1] = point[1]+y[1];                                                                    
            points_aux[2] = point[2]+y[2];                                                                    
                                                                                                              
            nrm = norm(z);                                                                                    
            nrm2 = norm(point);
            ev_basis_function_maxwell(points_aux, phi, p0, A, len, n, 1, opt, "RWG");
            res_j[0] =0.25*phi[0]*cexp((1*I)*kappa*nrm2)*J*weights[index]/(nrm*4.0*PI);                          
            res_j[1] =0.25*phi[1]*cexp((1*I)*kappa*nrm2)*J*weights[index]/(nrm*4.0*PI);
            res_j[2] =0.25*phi[2]*cexp((1*I)*kappa*nrm2)*J*weights[index]/(nrm*4.0*PI);
            res[0] += creal(res_j[0]);                                                                           
            res[1] += creal(res_j[1]); 
            res[2] += creal(res_j[2]);  
            res[3] += cimag(res_j[0]); 
            res[4] += cimag(res_j[1]);   
            res[5] += cimag(res_j[2]);                                                 
            index++;
         }                                                                                                     
    }                                                                                                         
                                                                                                              
}     

void evaluate_sing_kern_maxwell_hypersingular(double *y,  double *M, double *V, double a, double len, int order, int opt, double *res, double kappa, char *tp)
{                                                                                                             
    int       i, j, index;                                                                                    
    double       point[3];                                                                                    
    double           z[3];                                                                                    
    double    weights[16];                                                                                    
    double eta[4], chi[4];                                                                                    
    double chi2[1],eta2[1];                                                                                   
    double       nrm,nrm2;                                                                                    
    double  complex res_j;                                                                                    
    double              J;                                                                                    
    double           T[9];                                                                                    
    double jac[3], phi[1];                                                                                    
    int             sz[2];                                                                                    
    double  points_aux[3];                                                                                    
                                                                                                              
    jac[0] = (M[1]-y[1])*(M[5]-M[2]) - (M[2]-y[2])*(M[4]-M[1]);                                               
    jac[1] = (M[0]-y[0])*(M[5]-M[2]) - (M[2]-y[2])*(M[3]-M[0]);                                               
    jac[2] = (M[0]-y[0])*(M[4]-M[1]) - (M[1]-y[1])*(M[3]-M[0]);                                               
    J = norm(jac);                                                                                            
    quad_rules_DUAL( eta, chi, sz, order );                                                                   
    weights_DUAL( weights, order, sz[0] );                                                                    
                                                                                                              
    index = 0;                                                                                                
    phi[0] = 1.0;                                                                                             
    res[0] = 0.0;                                                                                             
    res[1] = 0.0;                                                                                             
                                                                                                              
    for(i = 0; i<sz[0]; i++)                                                                                  
    {                                                                                                         
        eta[i] = (eta[i]+1)*0.5;                                                                              
        chi[i] = (chi[i]+1)*0.5;                                                                              
    }                                                                                                         
    for( i = 0; i < sz[0]; i++ )                                                                              
    {                                                                                                         
        for( j=0; j < sz[1]; j++ )                                                                            
        {                                                                                                     
            z[0] = -y[0]+(1-eta[j])*M[0]+eta[j]*M[3];                                                         
            z[1] = -y[1]+(1-eta[j])*M[1]+eta[j]*M[4];                                                         
            z[2] = -y[2]+(1-eta[j])*M[2]+eta[j]*M[5];                                                         
                                                                                                              
            point[0] = chi[i]*z[0];                                                                           
            point[1] = chi[i]*z[1];                                                                           
            point[2] = chi[i]*z[2];                                                                           
                                                                                                              
            points_aux[0] = point[0]+y[0];                                                                    
            points_aux[1] = point[1]+y[1];                                                                    
            points_aux[2] = point[2]+y[2];                                                                    
                                                                                                              
            nrm = norm(z);                                                                                    
            nrm2 = norm(point); 
            divergence(a,len, NULL,"RWG", opt, phi, 1);
            res_j =0.25*phi[0]*cexp((1*I)*kappa*nrm2)*J*weights[index]/(nrm*4.0*PI);                          
            res[0] += creal(res_j);                                                                           
            res[1] += cimag(res_j);                                                                           
            index++;                                                                                          
        }                                                                                                     
    }
}
void singular(double *P0, double *V, double *weights, double *res, double *IE, int *sz,  int order, int dof, double kappa, char *tp, double phi)
{
    double y[3], M[6], A[9], a;
    double cres[2];
    int k;
    
    cres[0] = 0.0;
    cres[1] = 0.0;
    res[0] = 0.0;
    res[1] = 0.0;
    for( k = 0; k < sz[0]; k++)
    {
            y[0] = P0[3*k];                                                                                       
            y[1] = P0[3*k+1];                                                                                     
            y[2] = P0[3*k+2];
            A[0] = P0[3*k];
            A[1] = P0[3*k+1];                                                                                 
            A[2] = P0[3*k+2];
            M[0] = V[0];
            M[1] = V[1];
            M[2] = V[2];
            M[3] = V[3];
            M[4] = V[4];
            M[5] = V[5];
            A[3] = V[0];                                                                                      
            A[4] = V[1];                                                                                      
            A[5] = V[2];                                                                                      
            A[6] = V[3];                                                                                      
            A[7] = V[4];                                                                                      
            A[8] = V[5];
            a = area(A);
            evaluate_sing_kern(y, M, V, order, cres, dof, kappa, tp);
            res[0]+= a*cres[0];
            res[1]+= a*cres[1];
            M[0] = V[3];
            M[1] = V[4];
            M[2] = V[5];
            M[3] = V[6];
            M[4] = V[7];
            M[5] = V[8];
            A[3] = V[3];                                                                                      
            A[4] = V[4];                                                                                      
            A[5] = V[5];                                                                                      
            A[6] = V[6];                                                                                      
            A[7] = V[7];                                                                                      
            A[8] = V[8];
            a = area(A);
            evaluate_sing_kern(y, M, V, order, cres, dof, kappa, tp);
            res[0]+= a*cres[0];
            res[1]+= a*cres[1];
            M[0] = V[6];                                                                                      
            M[1] = V[7];                                                                                      
            M[2] = V[8];                                                                                      
            M[3] = V[0];                                                                                      
            M[4] = V[1];                                                                                      
            M[5] = V[2];
            A[3] = V[6];                                                                                      
            A[4] = V[7];                                                                                      
            A[5] = V[8];                                                                                      
            A[6] = V[0];                                                                                      
            A[7] = V[1];                                                                                      
            A[8] = V[2];
            a = area(A);
            evaluate_sing_kern(y, M, V, order, cres, dof, kappa, tp);
            res[0]+= cres[0]*a;                                                                                 
            res[1]+= cres[1]*a;                                                                                 
            res[0] = res[0]*phi;                                                                          
            res[1] = res[1]*phi;
        
       
      
    }
}

void singular_maxwell_weakly(double *P0, double *V, double *p0, double *weights, double A_t, double len, gsl_vector *n, double *res, int *sz,  int order, int opt, double kappa,  double *phi)
{                                                                                                             
    double y[3], M[6], A[9], a;                                                                               
    double cres[6] = {0.0,0.0,0.0,0.0,0.0,0.0};  
    double res_aux[6] = {0.0,0.0,0.0,0.0,0.0,0.0};                                                                                         
    int k;                                                                             
    double ones[3] ={len/A_t,len/A_t,len/A_t};
                                                                                                                                                                                                       
    res[0] = 0.0;                                                                                             
    res[1] = 0.0;                                                                                             
    for( k = 0; k < sz[0]; k++)                                                                               
    {                                                                                                         
            y[0] = P0[3*k];                                                                                       
            y[1] = P0[3*k+1];                                                                                     
            y[2] = P0[3*k+2];                                                                                 
            A[0] = P0[3*k];                                                                                   
            A[1] = P0[3*k+1];                                                                                 
            A[2] = P0[3*k+2];                                                                                 
            M[0] = V[0];                                                                                      
            M[1] = V[1];                                                                                      
            M[2] = V[2];                                                                                      
            M[3] = V[3];                                                                                      
            M[4] = V[4];                                                                                      
            M[5] = V[5];                                                                                      
            A[3] = V[0];                                                                                      
            A[4] = V[1];                                                                                      
            A[5] = V[2];                                                                                      
            A[6] = V[3];                                                                                      
            A[7] = V[4];                                                                                      
            A[8] = V[5];                                                                                      
            a = area(A);
            evaluate_sing_kern_maxwell_weakly(y, M, p0, A_t, len, n, opt, order, cres, kappa);                                         
            res_aux[0]+= a*cres[0];                                                                               
            res_aux[1]+= a*cres[1];                                                      
            res_aux[2]+= a*cres[2];
            res_aux[3]+= a*cres[3];
            res_aux[4]+= a*cres[4];
            res_aux[5]+= a*cres[5];                         
            M[0] = V[3];                                                                                      
            M[1] = V[4];                                                                                      
            M[2] = V[5];                                                                                      
            M[3] = V[6];                                                                                      
            M[4] = V[7];                                                                                      
            M[5] = V[8];                                                                                      
            A[3] = V[3];                                                                                      
            A[4] = V[4];                                                                                      
            A[5] = V[5];                                                                                      
            A[6] = V[6];                                                                                      
            A[7] = V[7];                                                                                      
            A[8] = V[8];                                                                                      
            a = area(A);
            evaluate_sing_kern_maxwell_weakly(y, M, p0, A_t, len, n, opt, order, cres, kappa);                                                                     
            res_aux[0]+= a*cres[0];                                                                               
            res_aux[1]+= a*cres[1];                                                                               
            res_aux[2]+= a*cres[2];                                                                               
            res_aux[3]+= a*cres[3];                                                                               
            res_aux[4]+= a*cres[4];                                                                               
            res_aux[5]+= a*cres[5];                                                   
            M[0] = V[6];                                                                                      
            M[1] = V[7];                                                                                      
            M[2] = V[8];                                                                                      
            M[3] = V[0];                                                                                      
            M[4] = V[1];                                                                                      
            M[5] = V[2];                                                                                      
            A[3] = V[6];                                                                                      
            A[4] = V[7];                                                                                      
            A[5] = V[8];                                                                                      
            A[6] = V[0];                                                                                      
            A[7] = V[1];                                                                                      
            A[8] = V[2];                                                                                      
            a = area(A);
            evaluate_sing_kern_maxwell_weakly(y, M, p0, A_t, len, n, opt, order, cres, kappa);                                                                     
            res_aux[0]+= a*cres[0];                                                                               
            res_aux[1]+= a*cres[1];                                                                               
            res_aux[2]+= a*cres[2];                                                                               
            res_aux[3]+= a*cres[3];                                                                               
            res_aux[4]+= a*cres[4];                                                                               
            res_aux[5]+= a*cres[5];                                                                           
            res[0] = dot_prod(3, 3*k, 0, phi, res_aux);                                                                              
            res[1] = dot_prod(3, 3*k, 3, phi, res_aux);                                                                              
                                                                                                              
                                                                                                              
                                                                                                              
    }                                                                                                         
} 



void singular_maxwell_weakly2(double kappa, int *order, double *v, double *p01, double *p02, double *res, gsl_vector *n, double A_t, double len1, double len2, int opt1, int opt2)
{
    double eta1[4], eta2[4], eta3[4], chi[4], w1[4], w2[4], w3[4], w4[4];
    int sz[2], sz2[2];
    double z[2], x[2], y[2], z2[2], point1[3], point2[3], point3[3], point4[3];
    double vert1[3], vert2[3];
    int i,j,k,l;
    double n1, n2, Jac;
    double complex res_j = 0.0+0.0*I;
    double complex k1, k2;
    gsl_matrix *G = gsl_matrix_alloc(2, 2);                                                                   
    gsl_matrix *J = gsl_matrix_alloc(3, 2);                                                                   
    double det; 
    res[0] = 0.0;
    res[1] = 0.0;
    quad_rules_ab( 0.0, 1.0, eta1, eta2, w1, w2, sz, order[3]);
    quad_rules_ab( 0.0, 1.0, eta3, chi, w3, w4, sz2, order[3]);    
    double phii[3], phij[3];
    

    vert1[0] = v[0]-v[3];                                                                                     
    vert1[1] = v[1]-v[4];                                                                                     
    vert1[2] = v[2]-v[5];                                                                                     
    vert2[0] = v[3]-v[6];                                                                                     
    vert2[1] = v[4]-v[7];                                                                                     
    vert2[2] = v[5]-v[8];                                                                                     
                                                                                                              
    for(i=0; i< sz2[0]; i++)                                                                                  
    {                                                                                                         
        for(j=0; j< sz2[1]; j++)                                                                              
        {                                                                                                     
            for(k=0; k< sz[0]; k++)                                                                           
            {                                                                                                 
                for(l=0; l< sz[1]; l++)                                                                       
                {                                                                                             
                    point1[0] = vert1[0]+vert2[0]*eta1[l];                                                    
                    point1[1] = vert1[1]+vert2[1]*eta1[l];                                                    
                    point1[2] = vert1[2]+vert2[2]*eta1[l];                                                    
                    z[0] = chi[i];                                                                            
                    z[1] = chi[i]*eta1[l];                                                                    
                    y[0] = (1-z[0])*eta2[k];                                                                  
                    y[1] = y[0]*eta3[j];                                                                      
                    x[0] = z[0] + y[0];                                                                       
                    x[1] = z[1] + y[1];                                                                       
                    Jac = (1-z[0])*y[0];                                     
                    point3[0] = (1-x[0])*v[0]+(x[0]-x[1])*v[3]+x[1]*v[6];                         
                    point3[1] = (1-x[0])*v[1]+(x[0]-x[1])*v[4]+x[1]*v[7];                         
                    point3[2] = (1-x[0])*v[2]+(x[0]-x[1])*v[5]+x[1]*v[8];                         
                    point4[0] = (1-y[0])*v[0]+(y[0]-y[1])*v[3]+y[1]*v[6];                         
                    point4[1] = (1-y[0])*v[1]+(y[0]-y[1])*v[4]+y[1]*v[7];                         
                    point4[2] = (1-y[0])*v[2]+(y[0]-y[1])*v[5]+y[1]*v[8];                                  
                    n1 = sqrt( chi[i]*chi[i]*(point1[0]*point1[0]+point1[1]*point1[1]+point1[2]*point1[2]) ); 
                    n2 = sqrt( point1[0]*point1[0]+point1[1]*point1[1]+point1[2]*point1[2] );      
                    ev_basis_function_maxwell(point3, phii, p01, A_t, len1, NULL, 1, opt1, "RWG");             
                    ev_basis_function_maxwell(point4, phij, p02, A_t, len2, NULL, 1, opt2, "RWG");            
                    k1 = Jac*cexp((1*I)*kappa*n1)*dot_prod(3,0,0,phii,phij)/(n2*4.0*PI);
                    ev_basis_function_maxwell(point4, phii, p01, A_t, len1, NULL, 1, opt1, "RWG");            
                    ev_basis_function_maxwell(point3, phij, p02, A_t, len2, NULL, 1, opt2, "RWG");
                    k2 = Jac*cexp((1*I)*kappa*n1)*dot_prod(3,0,0,phii,phij)/(n2*4.0*PI);         
                    res_j+= w1[l]*w2[k]*w3[j]*w4[i]*(k1+k2);                                                                                         
                                                                                                              
                }                                                                                             
            }                                                                                                 
        }                                                                                                     
    }                                                                                                         
    res[0] += creal(res_j);                                                                                   
    res[1] += cimag(res_j);                                                                                   
    res_j = 0.0+0.0*I; 

    for(i=0; i< sz2[0]; i++)                                                                                  
    {                                                                                                         
        for(j=0; j< sz2[1]; j++)                                                                              
        {                                                                                                     
            for(k=0; k< sz[0]; k++)                                                                           
            {                                                                                                 
                for(l=0; l< sz[1]; l++)                                                                       
                {                                                                                             
                    point1[0] = eta1[l]*vert1[0]+vert2[0]*(eta1[l]-1);                                        
                    point1[1] = eta1[l]*vert1[1]+vert2[1]*(eta1[l]-1);                                        
                    point1[2] = eta1[l]*vert1[2]+vert2[2]*(eta1[l]-1);                                        
                    z[0] = chi[i]*eta1[l];                                                                    
                    z[1] = chi[i]*(eta1[l]-1);                                                                
                    y[0] = (1-z[0]+z[1])*eta2[k]-z[1];                                                        
                    y[1] = (z[1]+y[0])*eta3[j]-z[1];                                                          
                    x[0] = z[0] + y[0];                                                                       
                    x[1] = z[1] + y[1];                                                                       
                    Jac = (1-z[0]+z[1])*(y[0]+z[1]);                                  
                    point3[0] = (1-x[0])*v[0]+(x[0]-x[1])*v[3]+x[1]*v[6];                         
                    point3[1] = (1-x[0])*v[1]+(x[0]-x[1])*v[4]+x[1]*v[7];                         
                    point3[2] = (1-x[0])*v[2]+(x[0]-x[1])*v[5]+x[1]*v[8];                         
                    point4[0] = (1-y[0])*v[0]+(y[0]-y[1])*v[3]+y[1]*v[6];                         
                    point4[1] = (1-y[0])*v[1]+(y[0]-y[1])*v[4]+y[1]*v[7];                         
                    point4[2] = (1-y[0])*v[2]+(y[0]-y[1])*v[5]+y[1]*v[8];                         
                    n1 = sqrt( chi[i]*chi[i]*(point1[0]*point1[0]+point1[1]*point1[1]+point1[2]*point1[2]) ); 
                    n2 = sqrt( point1[0]*point1[0]+point1[1]*point1[1]+point1[2]*point1[2] );      
                    ev_basis_function_maxwell(point3, phii, p01, A_t, len1, NULL, 1, opt1, "RWG");             
                    ev_basis_function_maxwell(point4, phij, p02, A_t, len2, NULL, 1, opt2, "RWG"); 
                    k1 = Jac*cexp((1*I)*kappa*n1)*dot_prod(3,0,0,phii,phij)/(n2*4.0*PI);                      
                    ev_basis_function_maxwell(point4, phii, p01, A_t, len1, NULL, 1, opt1, "RWG");            
                    ev_basis_function_maxwell(point3, phij, p02, A_t, len2, NULL, 1, opt2, "RWG");            
                    k2 = Jac*cexp((1*I)*kappa*n1)*dot_prod(3,0,0,phii,phij)/(n2*4.0*PI);                      
                    res_j+= w1[l]*w2[k]*w3[j]*w4[i]*(k1+k2);              
                                                                                                              
                }                                                                                             
            }                                                                                                 
        }                                                                                                     
    }                                                                                                         
    res[0] += creal(res_j);                                                                                   
    res[1] += cimag(res_j);                                                                                   
    res_j = 0.0+0.0*I;                                                                                        
    for(i=0; i< sz2[0]; i++)                                                                                  
    {                                                                                                         
        for(j=0; j< sz2[1]; j++)                                                                              
        {                                                                                                     
            for(k=0; k< sz[0]; k++)                                                                           
            {                                                                                                 
                for(l=0; l< sz[1]; l++)                                                                       
                {                                                                                             
                    point1[0] = eta1[l]*vert1[0]+vert2[0];                                                    
                    point1[1] = eta1[l]*vert1[1]+vert2[1];                                                    
                    point1[2] = eta1[l]*vert1[2]+vert2[2];                                                    
                    z[0] = chi[i]*eta1[l];                                                                    
                    z[1] = chi[i];                                                                            
                    y[0] = (1-z[1])*eta2[k]+z[1]-z[0];                                                        
                    y[1] = (y[0]-z[1]+z[0])*eta3[j];                                                          
                    x[0] = z[0] + y[0];                                                                       
                    x[1] = z[1] + y[1];                                                                       
                    Jac = (1-z[1])*(y[0]-z[1]+z[0]);                                  
                    point3[0] = (1-x[0])*v[0]+(x[0]-x[1])*v[3]+x[1]*v[6];                         
                    point3[1] = (1-x[0])*v[1]+(x[0]-x[1])*v[4]+x[1]*v[7];                         
                    point3[2] = (1-x[0])*v[2]+(x[0]-x[1])*v[5]+x[1]*v[8];                         
                    point4[0] = (1-y[0])*v[0]+(y[0]-y[1])*v[3]+y[1]*v[6];                         
                    point4[1] = (1-y[0])*v[1]+(y[0]-y[1])*v[4]+y[1]*v[7];                         
                    point4[2] = (1-y[0])*v[2]+(y[0]-y[1])*v[5]+y[1]*v[8];                         

                    n1 = sqrt( chi[i]*chi[i]*(point1[0]*point1[0]+point1[1]*point1[1]+point1[2]*point1[2]) ); 
                    n2 = sqrt( point1[0]*point1[0]+point1[1]*point1[1]+point1[2]*point1[2] );  
                    ev_basis_function_maxwell(point3, phii, p01, A_t, len1, NULL, 1, opt1, "RWG");             
                    ev_basis_function_maxwell(point4, phij, p02, A_t, len2, NULL, 1, opt2, "RWG");                
                    k1 = Jac*cexp((1*I)*kappa*n1)*dot_prod(3,0,0,phii,phij)/(n2*4.0*PI);                      
                    ev_basis_function_maxwell(point4, phii, p01, A_t, len1, NULL, 1, opt1, "RWG");            
                    ev_basis_function_maxwell(point3, phij, p02, A_t, len2, NULL, 1, opt2, "RWG");            
                    k2 = Jac*cexp((1*I)*kappa*n1)*dot_prod(3,0,0,phii,phij)/(n2*4.0*PI);                      
                    res_j+= w1[l]*w2[k]*w3[j]*w4[i]*(k1+k2);    
                                                                                                              
                                                                                                              
                }                                                                                             
            }                                                                                                 
        }                                                                                                     
    }
    /*gsl_matrix_set(J,0,0, v[3]-v[0]);                                                                         
    gsl_matrix_set(J,0,1, v[6]-v[3]);                                                                         
    gsl_matrix_set(J,1,0, v[4]-v[1]);                                                                         
    gsl_matrix_set(J,1,1, v[7]-v[4]);                                                                         
    gsl_matrix_set(J,2,0, v[5]-v[2]);                                                                         
    gsl_matrix_set(J,2,1, v[8]-v[5]);                                                                         
    gsl_blas_dgemm(CblasTrans,CblasNoTrans, 1.0, J, J, 0.0, G);                                               
    det = fabs(gsl_matrix_get(G,0,0)*gsl_matrix_get(G,1,1)- gsl_matrix_get(G,0,1)*gsl_matrix_get(G,1,0));     
    printf("%f\n",det);    
    for(i=0; i< sz2[0]; i++)
    {
        for(j=0; j< sz2[1]; j++)
        {
            for(k=0; k< sz[0]; k++)
            {
                for(l=0; l< sz[1]; l++)
                {
                    x[0] = 1;
                    x[1] = 1-eta1[l]+eta2[k]*eta1[l];
                    z[0] = eta1[l]*eta2[k]*eta3[j];
                    z[1] = eta2[k]*eta1[l];
                    z2[0] = eta3[j];
                    z2[1] = 1;
                    y[0] = x[0]-z[0];
                    y[1] = x[1]-z[1];
                    //point1[0] = (1-chi[i]*x[0]-chi[i]*x[1])*v[0] + x[0]*chi[i]*v[3] + chi[i]*x[1]*v[6];
                    //point1[1] = (1-chi[i]*x[0]-chi[i]*x[1])*v[1] + x[0]*chi[i]*v[4] + chi[i]*x[1]*v[7]; 
                    //point1[2] = (1-chi[i]*x[0]-chi[i]*x[1])*v[2] + x[0]*chi[i]*v[5] + chi[i]*x[1]*v[8]; 
                    //point2[0] = (1-chi[i]*y[0]-chi[i]*y[1])*v[0] + y[0]*chi[i]*v[3] + chi[i]*y[1]*v[6];      
                    //point2[1] = (1-chi[i]*y[0]-chi[i]*y[1])*v[1] + y[0]*chi[i]*v[4] + chi[i]*y[1]*v[7];      
                    //point2[2] = (1-chi[i]*y[0]-chi[i]*y[1])*v[2] + y[0]*chi[i]*v[5] + chi[i]*y[1]*v[8];
                    
                    point1[0] = v[0] + (v[3]-v[0])*x[0]*chi[i] + x[1]*chi[i]*(v[6]-v[3]);
                    point1[1] = v[1] + (v[4]-v[1])*x[0]*chi[i] + x[1]*chi[i]*(v[7]-v[4]);
                    point1[2] = v[2] + (v[5]-v[2])*x[0]*chi[i] + x[1]*chi[i]*(v[8]-v[5]);
                    point2[0] = v[0] + (v[3]-v[0])*y[0]*chi[i] + y[1]*chi[i]*(v[6]-v[3]);                     
                    point2[1] = v[1] + (v[4]-v[1])*y[0]*chi[i] + y[1]*chi[i]*(v[7]-v[4]);                     
                    point2[2] = v[2] + (v[5]-v[2])*y[0]*chi[i] + y[1]*chi[i]*(v[8]-v[5]); 
                    point3[0] = (v[3]-v[0])*chi[i]*z[0]+(v[6]-v[3])*chi[i]*z[1];                              
                    point3[1] = (v[4]-v[1])*chi[i]*z[0]+(v[7]-v[4])*chi[i]*z[1];                              
                    point3[2] = (v[5]-v[2])*chi[i]*z[0]+(v[8]-v[5])*chi[i]*z[1];                              
                    point4[0] = (v[3]-v[0])*z2[0]+(v[6]-v[3])*z2[1];                                          
                    point4[1] = (v[4]-v[1])*z2[0]+(v[7]-v[4])*z2[1];                                          
                    point4[2] = (v[5]-v[2])*z2[0]+(v[8]-v[5])*z2[1];
                    ev_basis_function_maxwell(point1, phii, p01, A_t, len1, n, 1, opt1, "RWG");
                    ev_basis_function_maxwell(point2, phij, p02, A_t, len2, n, 1, opt2, "RWG"); 
                    //n1 = sqrt( (z[0]*z[0]+z[1]*z[1])*chi[i]*chi[i] );
                    //n2 = sqrt( z2[0]*z2[0]+z2[1]*z2[1] );
                    n1 = sqrt((point3[0]*point3[0]+point3[1]*point3[1]+point3[2]*point3[2])*chi[i]*chi[i]);
                    n2 = sqrt((point4[0]*point4[0]+point4[1]*point4[1]+point4[2]*point4[2]));
                    k1 = chi[i]*chi[i]*eta1[l]*cexp((1*I)*kappa*n1)*dot_prod(3,0,0,phii,phij)/(n2*4.0*PI);
                    ev_basis_function_maxwell(point2, phii, p01, A_t, len2, n, 1, opt1, "RWG");                
                    ev_basis_function_maxwell(point1, phij, p02, A_t, len1, n, 1, opt2, "RWG"); 
                    k2 = chi[i]*chi[i]*eta1[l]*cexp((1*I)*kappa*n1)*dot_prod(3,0,0,phii,phij)/(n2*4.0*PI);
                    res_j+= w1[l]*w2[k]*w3[j]*w4[i]*(k1+k2);
                    
                }
            }
        }
    }
    res[0] += creal(res_j);
    res[1] += cimag(res_j);
    printf("%f %f\n", res[0], res[1]);
    res_j = 0.0+0.0*I; 
    for(i=0; i< sz2[0]; i++)                                                                                  
    {                                                                                                         
        for(j=0; j< sz2[1]; j++)                                                                              
        {                                                                                                     
            for(k=0; k< sz[0]; k++)                                                                           
            {                                                                                                 
                for(l=0; l< sz[1]; l++)                                                                       
                {                                                                                             
                    x[0] = 1;                                                                                 
                    x[1] = eta1[l]-eta1[l]*eta2[k]+eta1[l]*eta2[k]*eta3[j];                                                         
                    z[0] = eta1[l]*eta2[k];                                                    
                    z[1] = eta2[k]*eta1[l]*eta3[j];                                                                   
                    z2[0] = 1;                                                                
                    z2[1] = eta3[j];                                                                                
                    y[0] = x[0]-z[0];                                                                         
                    y[1] = x[1]-z[1];                                                                         
                    //point1[0] = (1-chi[i]*x[0]-chi[i]*x[1])*v[0] + x[0]*chi[i]*v[3] + chi[i]*x[1]*v[6];      
                    //point1[1] = (1-chi[i]*x[0]-chi[i]*x[1])*v[1] + x[0]*chi[i]*v[4] + chi[i]*x[1]*v[7];      
                    //point1[2] = (1-chi[i]*x[0]-chi[i]*x[1])*v[2] + x[0]*chi[i]*v[5] + chi[i]*x[1]*v[8];      
                    //point2[0] = (1-chi[i]*y[0]-chi[i]*y[1])*v[0] + y[0]*chi[i]*v[3] + chi[i]*y[1]*v[6];      
                    //point2[1] = (1-chi[i]*y[0]-chi[i]*y[1])*v[1] + y[0]*chi[i]*v[4] + chi[i]*y[1]*v[7];      
                    //point2[2] = (1-chi[i]*y[0]-chi[i]*y[1])*v[2] + y[0]*chi[i]*v[5] + chi[i]*y[1]*v[8]; 
                    
                    point1[0] = v[0] + (v[3]-v[0])*x[0]*chi[i] + x[1]*chi[i]*(v[6]-v[3]);                     
                    point1[1] = v[1] + (v[4]-v[1])*x[0]*chi[i] + x[1]*chi[i]*(v[7]-v[4]);                     
                    point1[2] = v[2] + (v[5]-v[2])*x[0]*chi[i] + x[1]*chi[i]*(v[8]-v[5]);                     
                    point2[0] = v[0] + (v[3]-v[0])*y[0]*chi[i] + y[1]*chi[i]*(v[6]-v[3]);                     
                    point2[1] = v[1] + (v[4]-v[1])*y[0]*chi[i] + y[1]*chi[i]*(v[7]-v[4]);                     
                    point2[2] = v[2] + (v[5]-v[2])*y[0]*chi[i] + y[1]*chi[i]*(v[8]-v[5]);
                    point3[0] = (v[3]-v[0])*chi[i]*z[0]+(v[6]-v[3])*chi[i]*z[1];                              
                    point3[1] = (v[4]-v[1])*chi[i]*z[0]+(v[7]-v[4])*chi[i]*z[1];                              
                    point3[2] = (v[5]-v[2])*chi[i]*z[0]+(v[8]-v[5])*chi[i]*z[1];                              
                    point4[0] = (v[3]-v[0])*z2[0]+(v[6]-v[3])*z2[1];                                          
                    point4[1] = (v[4]-v[1])*z2[0]+(v[7]-v[4])*z2[1];                                          
                    point4[2] = (v[5]-v[2])*z2[0]+(v[8]-v[5])*z2[1];      
                    ev_basis_function_maxwell(point1, phii, p01, A_t, len1, NULL, 1, opt1, "RWG");             
                    ev_basis_function_maxwell(point2, phij, p02, A_t, len2, NULL, 1, opt2, "RWG");             
                    //n1 = sqrt( (z[0]*z[0]+z[1]*z[1])*chi[i]*chi[i] );                                                         
                    //n2 = sqrt( z2[0]*z2[0]+z2[1]*z2[1]);                                                     
                    n1 = sqrt((point3[0]*point3[0]+point3[1]*point3[1]+point3[2]*point3[2])*chi[i]*chi[i]);   
                    n2 = sqrt((point4[0]*point4[0]+point4[1]*point4[1]+point4[2]*point4[2]));
                    k1 = chi[i]*chi[i]*eta1[l]*cexp((1*I)*kappa*n1)*dot_prod(3,0,0,phii,phij)/(n2*4.0*PI);    
                    ev_basis_function_maxwell(point2, phii, p01, A_t, len2, n, 1, opt1, "RWG");                
                    ev_basis_function_maxwell(point1, phij, p02, A_t, len1, n, 1, opt2, "RWG");                 
                    k2 = chi[i]*chi[i]*eta1[l]*cexp((1*I)*kappa*n1)*dot_prod(3,0,0,phii,phij)/(n2*4.0*PI);    
                    res_j+= w1[l]*w2[k]*w3[j]*w4[i]*(k1+k2);   
                                                                                                              
                }                                                                                             
            }                                                                                                 
        }                                                                                                     
    }
    res[0] += creal(res_j);                                                                                   
    res[1] += cimag(res_j);         
    printf("%f %f\n", res[0], res[1]); 
    res_j = 0.0+0.0*I;
    for(i=0; i< sz2[0]; i++)                                                                                  
    {                                                                                                         
        for(j=0; j< sz2[1]; j++)                                                                              
        {                                                                                                     
            for(k=0; k< sz[0]; k++)                                                                           
            {                                                                                                 
                for(l=0; l< sz[1]; l++)                                                                       
                {                                                                                             
                    x[0] = 1-eta1[l]*eta2[k]*eta3[j];                                                                                 
                    x[1] = eta1[l]-eta1[l]*eta2[k]*eta3[j];                                   
                    z[0] = -eta1[l]*eta2[k]*eta3[j];                                                                   
                    z[1] = eta1[l]*eta2[k]-eta2[k]*eta1[l]*eta3[j];                                                           
                    z2[0] = -eta3[j];                                                                                
                    z2[1] = 1-eta3[j];                                                                          
                    y[0] = x[0]-z[0];                                                                         
                    y[1] = x[1]-z[1];                                                                         
                    //point1[0] = (1-chi[i]*x[0]-chi[i]*x[1])*v[0] + x[0]*chi[i]*v[3] + chi[i]*x[1]*v[6];      
                    //point1[1] = (1-chi[i]*x[0]-chi[i]*x[1])*v[1] + x[0]*chi[i]*v[4] + chi[i]*x[1]*v[7];      
                    //point1[2] = (1-chi[i]*x[0]-chi[i]*x[1])*v[2] + x[0]*chi[i]*v[5] + chi[i]*x[1]*v[8];      
                    //point2[0] = (1-chi[i]*y[0]-chi[i]*y[1])*v[0] + y[0]*chi[i]*v[3] + chi[i]*y[1]*v[6];      
                    //point2[1] = (1-chi[i]*y[0]-chi[i]*y[1])*v[1] + y[0]*chi[i]*v[4] + chi[i]*y[1]*v[7];      
                    //point2[2] = (1-chi[i]*y[0]-chi[i]*y[1])*v[2] + y[0]*chi[i]*v[5] + chi[i]*y[1]*v[8];
                    
                    point1[0] = v[0] + (v[3]-v[0])*x[0]*chi[i] + x[1]*chi[i]*(v[6]-v[3]);                     
                    point1[1] = v[1] + (v[4]-v[1])*x[0]*chi[i] + x[1]*chi[i]*(v[7]-v[4]);                     
                    point1[2] = v[2] + (v[5]-v[2])*x[0]*chi[i] + x[1]*chi[i]*(v[8]-v[5]);                     
                    point2[0] = v[0] + (v[3]-v[0])*y[0]*chi[i] + y[1]*chi[i]*(v[6]-v[3]);                     
                    point2[1] = v[1] + (v[4]-v[1])*y[0]*chi[i] + y[1]*chi[i]*(v[7]-v[4]);                     
                    point2[2] = v[2] + (v[5]-v[2])*y[0]*chi[i] + y[1]*chi[i]*(v[8]-v[5]);
                    point3[0] = (v[3]-v[0])*chi[i]*z[0]+(v[6]-v[3])*chi[i]*z[1];                              
                    point3[1] = (v[4]-v[1])*chi[i]*z[0]+(v[7]-v[4])*chi[i]*z[1];                              
                    point3[2] = (v[5]-v[2])*chi[i]*z[0]+(v[8]-v[5])*chi[i]*z[1];                              
                    point4[0] = (v[3]-v[0])*z2[0]+(v[6]-v[3])*z2[1];                                          
                    point4[1] = (v[4]-v[1])*z2[0]+(v[7]-v[4])*z2[1];                                          
                    point4[2] = (v[5]-v[2])*z2[0]+(v[8]-v[5])*z2[1];       
                    ev_basis_function_maxwell(point1, phii, p01, A_t, len1, NULL, 1, opt1, "RWG");             
                    ev_basis_function_maxwell(point2, phij, p02, A_t, len2, NULL, 1, opt2, "RWG");             
                    //n1 = sqrt( (z[0]*z[0]+z[1]*z[1])*chi[i]*chi[i] );                                                         
                    //n2 = sqrt( z2[0]*z2[0]+z2[1]*z2[1] );      
                    n1 = sqrt((point3[0]*point3[0]+point3[1]*point3[1]+point3[2]*point3[2])*chi[i]*chi[i]);   
                    n2 = sqrt((point4[0]*point4[0]+point4[1]*point4[1]+point4[2]*point4[2]));
                    k1 = chi[i]*chi[i]*eta1[l]*cexp((1*I)*kappa*n1)*dot_prod(3,0,0,phii,phij)/(n2*4.0*PI);    
                    ev_basis_function_maxwell(point2, phii, p01, A_t, len2, n, 1, opt1, "RWG");                
                    ev_basis_function_maxwell(point1, phij, p02, A_t, len1, n, 1, opt2, "RWG");                 
                    k2 = chi[i]*chi[i]*eta1[l]*cexp((1*I)*kappa*n1)*dot_prod(3,0,0,phii,phij)/(n2*4.0*PI);    
                    res_j+= w1[l]*w2[k]*w3[j]*w4[i]*(k1+k2);                                                  
                                                                                                              
                }                                                                                             
            }                                                                                                 
        }                                                                                                     
    } */
    res[0] += creal(res_j);                                                                                   
    res[1] += cimag(res_j);
    res[0] = res[0]*4*A_t*A_t;                                                                              
    res[1] = res[1]*4*A_t*A_t;

}

void singular_maxwell_hypersingular2(double kappa, int *order, double *v, double *p01, double *p02, double *res, gsl_vector *n, double A_t, double len1, double len2, int opt1, int opt2)
{                                                                                                             
    double eta1[4], eta2[4], eta3[4], chi[4], w1[4], w2[4], w3[4], w4[4];             
    int sz[2], sz2[2];                                                                                        
    double z[2], x[2], y[2], z2[2], point1[3], point2[3], point3[3], point4[3], vert1[3], vert2[3];                                                     
    int i,j,k,l;                                                                                              
    double n1, n2, k1, k2, Jac;                                                                                    
    double complex res_j = 0.0+0.0*I;
    gsl_matrix *G = gsl_matrix_alloc(2, 2);
    gsl_matrix *J = gsl_matrix_alloc(3, 2);
    double det;                                                                         
    res[0] = 0.0;                                                                                             
    res[1] = 0.0;                                                                                             
    quad_rules_ab( 0.0, 1.0, eta1, eta2, w1, w2, sz, order[3]);                                   
    quad_rules_ab( 0.0, 1.0, eta3, chi, w3, w4, sz2, order[3]);                                   
    double phii[sz[0]], phij[sz[0]];                                                    
    divergence(A_t,len1, NULL,"RWG", opt1, phii, 1);                                          
    divergence(A_t,len2, NULL,"RWG", opt2, phij, 1);


    vert1[0] = v[0]-v[3];
    vert1[1] = v[1]-v[4];
    vert1[2] = v[2]-v[5];
    vert2[0] = v[3]-v[6];                                                                                     
    vert2[1] = v[4]-v[7];                                                                                     
    vert2[2] = v[5]-v[8];

    for(i=0; i< sz2[0]; i++)                                                                                  
    {                                                                                                         
        for(j=0; j< sz2[1]; j++)                                                                              
        {                                                                                                     
            for(k=0; k< sz[0]; k++)                                                                           
            {                                                                                                 
                for(l=0; l< sz[1]; l++)                                                                       
                {                                                                                             
                    point1[0] = vert1[0]+vert2[0]*eta1[l];
                    point1[1] = vert1[1]+vert2[1]*eta1[l];
                    point1[2] = vert1[2]+vert2[2]*eta1[l];
                    z[0] = chi[i];
                    z[1] = chi[i]*eta1[l];
                    y[0] = (1-z[0])*eta2[k];
                    y[1] = y[0]*eta3[j]; 
                    x[0] = z[0] + y[0];
                    x[1] = z[1] + y[1];
                    Jac = (1-z[0])*y[0];
                    point3[0] = (1-x[0])*v[0]+(x[0]-x[1])*v[3]+x[1]*v[6];                                     
                    point3[1] = (1-x[0])*v[1]+(x[0]-x[1])*v[4]+x[1]*v[7];                                     
                    point3[2] = (1-x[0])*v[2]+(x[0]-x[1])*v[5]+x[1]*v[8];                                     
                    point4[0] = (1-y[0])*v[0]+(y[0]-y[1])*v[3]+y[1]*v[6];                                     
                    point4[1] = (1-y[0])*v[1]+(y[0]-y[1])*v[4]+y[1]*v[7];                                     
                    point4[2] = (1-y[0])*v[2]+(y[0]-y[1])*v[5]+y[1]*v[8];
                    n1 = sqrt( chi[i]*chi[i]*(point1[0]*point1[0]+point1[1]*point1[1]+point1[2]*point1[2]) );                 
                    n2 = sqrt( point1[0]*point1[0]+point1[1]*point1[1]+point1[2]*point1[2] );                 
                    res_j+= w1[l]*w2[k]*w3[j]*w4[i]*2*Jac*cexp((1*I)*kappa*n1)*phii[0]*phij[0]/(n2*4.0*PI);
                                                                                                              
                                                                                                              
                }                                                                                             
            }                                                                                                 
        }                                                                                                     
    }
    res[0] += creal(res_j);                                                                                   
    res[1] += cimag(res_j);                                                                                   
    res_j = 0.0+0.0*I; 
    
    for(i=0; i< sz2[0]; i++)                                                                                  
    {                                                                                                         
        for(j=0; j< sz2[1]; j++)                                                                              
        {                                                                                                     
            for(k=0; k< sz[0]; k++)                                                                           
            {                                                                                                 
                for(l=0; l< sz[1]; l++)                                                                       
                {                                                                                             
                    point1[0] = eta1[l]*vert1[0]+vert2[0]*(eta1[l]-1);                                                    
                    point1[1] = eta1[l]*vert1[1]+vert2[1]*(eta1[l]-1);                                                    
                    point1[2] = eta1[l]*vert1[2]+vert2[2]*(eta1[l]-1);                                                    
                    z[0] = chi[i]*eta1[l];                                                                            
                    z[1] = chi[i]*(eta1[l]-1);                                                                    
                    y[0] = (1-z[0]+z[1])*eta2[k]-z[1];                                                                  
                    y[1] = (z[1]+y[0])*eta3[j]-z[1];                                                                      
                    x[0] = z[0] + y[0];                                                                       
                    x[1] = z[1] + y[1];                                                                       
                    Jac = (1-z[0]+z[1])*(y[0]+z[1]);                     
                    point3[0] = (1-x[0])*v[0]+(x[0]-x[1])*v[3]+x[1]*v[6];                                     
                    point3[1] = (1-x[0])*v[1]+(x[0]-x[1])*v[4]+x[1]*v[7];                                     
                    point3[2] = (1-x[0])*v[2]+(x[0]-x[1])*v[5]+x[1]*v[8];                                     
                    point4[0] = (1-y[0])*v[0]+(y[0]-y[1])*v[3]+y[1]*v[6];                                     
                    point4[1] = (1-y[0])*v[1]+(y[0]-y[1])*v[4]+y[1]*v[7];                                     
                    point4[2] = (1-y[0])*v[2]+(y[0]-y[1])*v[5]+y[1]*v[8];               
                    n1 = sqrt( chi[i]*chi[i]*(point1[0]*point1[0]+point1[1]*point1[1]+point1[2]*point1[2]) );  
                    n2 = sqrt( point1[0]*point1[0]+point1[1]*point1[1]+point1[2]*point1[2] );                 
                    res_j+= w1[l]*w2[k]*w3[j]*w4[i]*2*Jac*cexp((1*I)*kappa*n1)*phii[0]*phij[0]/(n2*4.0*PI);     
                                                                                                              
                                                                                                              
                }                                                                                             
            }                                                                                                 
        }                                                                                                     
    }
    res[0] += creal(res_j);                                                                                   
    res[1] += cimag(res_j);                                                                                   
    res_j = 0.0+0.0*I; 
    for(i=0; i< sz2[0]; i++)                                                                                  
    {                                                                                                         
        for(j=0; j< sz2[1]; j++)                                                                              
        {                                                                                                     
            for(k=0; k< sz[0]; k++)                                                                           
            {                                                                                                 
                for(l=0; l< sz[1]; l++)                                                                       
                {                                                                                             
                    point1[0] = eta1[l]*vert1[0]+vert2[0];                                        
                    point1[1] = eta1[l]*vert1[1]+vert2[1];                                        
                    point1[2] = eta1[l]*vert1[2]+vert2[2];                                        
                    z[0] = chi[i]*eta1[l];                                                                    
                    z[1] = chi[i];                                                      
                    y[0] = (1-z[1])*eta2[k]+z[1]-z[0];                                                        
                    y[1] = (y[0]-z[1]+z[0])*eta3[j];                                                          
                    x[0] = z[0] + y[0];                                                                       
                    x[1] = z[1] + y[1];                                                                       
                    Jac = (1-z[1])*(y[0]-z[1]+z[0]);  
                    point3[0] = (1-x[0])*v[0]+(x[0]-x[1])*v[3]+x[1]*v[6];                                     
                    point3[1] = (1-x[0])*v[1]+(x[0]-x[1])*v[4]+x[1]*v[7];                                     
                    point3[2] = (1-x[0])*v[2]+(x[0]-x[1])*v[5]+x[1]*v[8];                                     
                    point4[0] = (1-y[0])*v[0]+(y[0]-y[1])*v[3]+y[1]*v[6];                                     
                    point4[1] = (1-y[0])*v[1]+(y[0]-y[1])*v[4]+y[1]*v[7];                                     
                    point4[2] = (1-y[0])*v[2]+(y[0]-y[1])*v[5]+y[1]*v[8]; 
                    n1 = sqrt( chi[i]*chi[i]*(point1[0]*point1[0]+point1[1]*point1[1]+point1[2]*point1[2]) );  
                    n2 = sqrt( point1[0]*point1[0]+point1[1]*point1[1]+point1[2]*point1[2] );                 
                    res_j+= w1[l]*w2[k]*w3[j]*w4[i]*2*Jac*cexp((1*I)*kappa*n1)*phii[0]*phij[0]/(n2*4.0*PI);     
                                                                                                              
                                                                                                              
                }                                                                                             
            }                                                                                                 
        }                                                                                                     
    } 
    /*gsl_matrix_set(J,0,0, v[3]-v[0]);                                                                         
    gsl_matrix_set(J,0,1, v[6]-v[3]);                                                                         
    gsl_matrix_set(J,1,0, v[4]-v[1]);                                                                         
    gsl_matrix_set(J,1,1, v[7]-v[4]);                                                                         
    gsl_matrix_set(J,2,0, v[5]-v[2]);                                                                         
    gsl_matrix_set(J,2,1, v[8]-v[5]);
    gsl_blas_dgemm(CblasTrans,CblasNoTrans, 1.0, J, J, 0.0, G); 
    det = fabs(gsl_matrix_get(G,0,0)*gsl_matrix_get(G,1,1)- gsl_matrix_get(G,0,1)*gsl_matrix_get(G,1,0));

    for(i=0; i< sz2[0]; i++)                                                                                  
    {                                                                                                         
        for(j=0; j< sz2[1]; j++)                                                                              
        {                                                                                                     
            for(k=0; k< sz[0]; k++)                                                                           
            {                                                                                                 
                for(l=0; l< sz[1]; l++)                                                                       
                {                                                                                             
                    x[0] = 1;                                                                                 
                    x[1] = 1-eta1[l]+eta2[k]*eta1[l];                                                         
                    z[0] = eta1[l]*eta2[k]*eta3[j];                                                           
                    z[1] = eta2[k]*eta1[l];                                                                   
                    z2[0] = eta3[j];                                                                          
                    z2[1] = 1;
                    point1[0] = (v[3]-v[0])*chi[i]*z[0]+(v[6]-v[3])*chi[i]*z[1];
                    point1[1] = (v[4]-v[1])*chi[i]*z[0]+(v[7]-v[4])*chi[i]*z[1];       
                    point1[2] = (v[5]-v[2])*chi[i]*z[0]+(v[8]-v[5])*chi[i]*z[1];     
                    point2[0] = (v[3]-v[0])*z2[0]+(v[6]-v[3])*z2[1];      
                    point2[1] = (v[4]-v[1])*z2[0]+(v[7]-v[4])*z2[1];       
                    point2[2] = (v[5]-v[2])*z2[0]+(v[8]-v[5])*z2[1];                                                                                
                    y[0] = x[0]-z[0];                                                                         
                    y[1] = x[1]-z[1]; 

                    n1 = sqrt( point1[0]*point1[0]+point1[1]*point1[1]+point1[2]*point1[2] );
                    n2 = sqrt( point2[0]*point2[0]+point2[1]*point2[1]+point2[2]*point2[2] ); 
                    res_j+= w1[l]*w2[k]*w3[j]*w4[i]*2*chi[i]*chi[i]*eta1[l]*cexp((1*I)*kappa*n1)*phii[0]*phij[0]/(n2*4.0*PI);


                }                                                                                             
            }                                                                                                 
        }                                                                                                     
    }
    res[0] += creal(res_j);                                                                                   
    res[1] += cimag(res_j);                                                                                   
    res_j = 0.0+0.0*I;                                                                                        
    for(i=0; i< sz2[0]; i++)                                                                                  
    {                                                                                                         
        for(j=0; j< sz2[1]; j++)                                                                              
        {                                                                                                     
            for(k=0; k< sz[0]; k++)                                                                           
            {                                                                                                 
                for(l=0; l< sz[1]; l++)                                                                       
                {                                                                                             
                    x[0] = 1;                                                                                 
                    x[1] = eta1[l]-eta1[l]*eta2[k]+eta1[l]*eta2[k]*eta3[j];                                   
                    z[0] = eta1[l]*eta2[k];                                                                   
                    z[1] = eta2[k]*eta1[l]*eta3[j];                                                           
                    z2[0] = 1;                                                                                
                    z2[1] = eta3[j];
                    point1[0] = (v[3]-v[0])*chi[i]*z[0]+(v[6]-v[3])*chi[i]*z[1];                              
                    point1[1] = (v[4]-v[1])*chi[i]*z[0]+(v[7]-v[4])*chi[i]*z[1];                                                                         
                    point1[2] = (v[5]-v[2])*chi[i]*z[0]+(v[8]-v[5])*chi[i]*z[1];                                                                         
                    point2[0] = (v[3]-v[0])*z2[0]+(v[6]-v[3])*z2[1];                                                                                    
                    point2[1] = (v[4]-v[1])*z2[0]+(v[7]-v[4])*z2[1];                                          
                    point2[2] = (v[5]-v[2])*z2[0]+(v[8]-v[5])*z2[1];                                                                          
                    y[0] = x[0]-z[0];                                                                         
                    y[1] = x[1]-z[1];                                                                         
                    n1 = sqrt( point1[0]*point1[0]+point1[1]*point1[1]+point1[2]*point1[2] );                 
                    n2 = sqrt( point2[0]*point2[0]+point2[1]*point2[1]+point2[2]*point2[2] );                    
                    res_j+= w1[l]*w2[k]*w3[j]*w4[i]*2*chi[i]*chi[i]*eta1[l]*cexp((1*I)*kappa*n1)*phii[0]*phij[0]/(n2*4.0*PI);
                                                                                                              
                }                                                                                             
            }                                                                                                 
        }                                                                                                     
    } 
    res[0] += creal(res_j);                                                                                   
    res[1] += cimag(res_j);                                                                                   
    res_j = 0.0+0.0*I;                                                                                        
    for(i=0; i< sz2[0]; i++)                                                                                  
    {                                                                                                         
        for(j=0; j< sz2[1]; j++)                                                                              
        {                                                                                                     
            for(k=0; k< sz[0]; k++)                                                                           
            {                                                                                                 
                for(l=0; l< sz[1]; l++)                                                                       
                {                                                                                             
                    x[0] = 1-eta1[l]*eta2[k]*eta3[j];                                                         
                    x[1] = eta1[l]-eta1[l]*eta2[k]*eta3[j];                                                   
                    z[0] = -eta1[l]*eta2[k]*eta3[j];                                                          
                    z[1] = eta1[l]*eta2[k]-eta2[k]*eta1[l]*eta3[j];                                           
                    z2[0] = -eta3[j];                                                                         
                    z2[1] = 1-eta3[j];                                                                        
                    point1[0] = (v[3]-v[0])*chi[i]*z[0]+(v[6]-v[3])*chi[i]*z[1];                              
                    point1[1] = (v[4]-v[1])*chi[i]*z[0]+(v[7]-v[4])*chi[i]*z[1];                              
                    point1[2] = (v[5]-v[2])*chi[i]*z[0]+(v[8]-v[5])*chi[i]*z[1];                              
                    point2[0] = (v[3]-v[0])*z2[0]+(v[6]-v[3])*z2[1];                                          
                    point2[1] = (v[4]-v[1])*z2[0]+(v[7]-v[4])*z2[1];                                          
                    point2[2] = (v[5]-v[2])*z2[0]+(v[8]-v[5])*z2[1];
                    y[0] = x[0]-z[0];                                                                         
                    y[1] = x[1]-z[1];                                                                         
                    n1 = sqrt( point1[0]*point1[0]+point1[1]*point1[1]+point1[2]*point1[2] );                 
                    n2 = sqrt( point2[0]*point2[0]+point2[1]*point2[1]+point2[2]*point2[2] ); 
                    res_j+= w1[l]*w2[k]*w3[j]*w4[i]*2*chi[i]*chi[i]*eta1[l]*cexp((1*I)*kappa*n1)*phii[0]*phij[0]/(n2*4.0*PI);
                                                                                                              
                }                                                                                             
            }                                                                                                 
        }                                                                                                     
    } */                                                                                                        
    res[0] += creal(res_j);                                                                                   
    res[1] += cimag(res_j);                                                                         
    res[0] = res[0]*4*A_t*A_t;                                                                              
    res[1] = res[1]*4*A_t*A_t;
                                                                                                              
}

void singular_maxwell_weakly_ce(double kappa, int *order, double *v1, double *v2,  double *p01, double *p02, double *res, gsl_vector *n, double len1, double len2, double A_t1, double A_t2, int *edges1, int *edges2, int opt1, int opt2)
{                                                                                                             
    double eta1[4], eta2[4], eta3[4], chi[4], w1[4], w2[4], w3[4], w4[4];       
    double vert1[9] = {v1[0], v1[1], v1[2], v1[3], v1[4], v1[5], v1[6], v1[7], v1[8]};                        
    double vert2[9] = {v2[0], v2[1], v2[2], v2[3], v2[4], v2[5], v2[6], v2[7], v2[8]};                                
    double  vert3[3], vert4[3], vert5[3];                                                                                
    int sz[2], sz2[2];                                                                                        
    double z[2], x[2], y[2], x2[2], y2[2], z2[3], point1[3], point2[3], point3[3], point4[3];                                        
    int i,j,k,l;                                                                                              
    double n1, n2, k1, k2;                                                                                    
    double det1, det2;                                                                                        
    gsl_matrix *G1 = gsl_matrix_alloc(2, 2);                                                                  
    gsl_matrix *J1 = gsl_matrix_alloc(3, 2);                                                                  
    gsl_matrix *G2 = gsl_matrix_alloc(2, 2);                                                                  
    gsl_matrix *J2 = gsl_matrix_alloc(3, 2);                                                                  
    double complex res_j = 0.0+0.0*I;                                                                         
    res[0] = 0.0;                                                                                             
    res[1] = 0.0;                                                                                             
    quad_rules_ab( 0.0, 1.0, eta1, eta2, w1, w2, sz, order[2]);                                               
    quad_rules_ab( 0.0, 1.0, eta3, chi, w3, w4, sz2, order[2]);                                               
    double phii[3], phij[3];                                                                                                                                  
                                                                                                              
    if(edges1[0]==edges2[0])                                                                                  
    {                                                                                                         
        vert1[0] = v1[0];                                                                                     
        vert1[1] = v1[1];                                                                                     
        vert1[2] = v1[2];                                                                                     
        vert1[3] = v1[3];                                                                                     
        vert1[4] = v1[4];                                                                                     
        vert1[5] = v1[5];                                                                                     
        vert1[6] = v1[6];                                                                                     
        vert1[7] = v1[7];                                                                                     
        vert1[8] = v1[8];                                                                                     
        vert2[0] = v1[0];                                                                                     
        vert2[1] = v1[1];                                                                                     
        vert2[2] = v1[2];                                                                                     
        vert2[3] = v1[3];                                                                                     
        vert2[4] = v1[4];                                                                                     
        vert2[5] = v1[5];                                                                                     
        vert2[6] = v2[6];                                                                                     
        vert2[7] = v2[7];                                                                                     
        vert2[8] = v2[8];                                                                                     
    }                                                                                                         
    else if(edges1[0]==edges2[1])                                                                             
    {                                                                                                         
        vert1[0] = v1[0];                                                                                     
        vert1[1] = v1[1];                                                                                     
        vert1[2] = v1[2];                                                                                     
        vert1[3] = v1[3];                                                                                     
        vert1[4] = v1[4];                                                                                     
        vert1[5] = v1[5];                                                                                     
        vert1[6] = v1[6];                                                                                     
        vert1[7] = v1[7];                                                                                     
        vert1[8] = v1[8];                                                                                     
        vert2[0] = v1[0];                                                                                     
        vert2[1] = v1[1];                                                                                     
        vert2[2] = v1[2];                                                                                     
        vert2[3] = v1[3];                                                                                     
        vert2[4] = v1[4]; 
        vert2[5] = v1[5];                                                                                     
        vert2[6] = v2[0];                                                                                     
        vert2[7] = v2[1];                                                                                     
        vert2[8] = v2[2];                                                                                     
    }                                                                                                         
    else if(edges1[0]==edges2[2])                                                                             
    {                                                                                                         
        vert1[0] = v1[0];                                                                                     
        vert1[1] = v1[1];                                                                                     
        vert1[2] = v1[2];                                                                                     
        vert1[3] = v1[3];                                                                                     
        vert1[4] = v1[4];                                                                                     
        vert1[5] = v1[5];                                                                                     
        vert1[6] = v1[6];                                                                                     
        vert1[7] = v1[7];                                                                                     
        vert1[8] = v1[8];                                                                                     
        vert2[0] = v1[0];                                                                                     
        vert2[1] = v1[1];                                                                                     
        vert2[2] = v1[2];                                                                                     
        vert2[3] = v1[3];                                                                                     
        vert2[4] = v1[4];                                                                                     
        vert2[5] = v1[5];                                                                                     
        vert2[6] = v2[3];                                                                                     
        vert2[7] = v2[4];                                                                                     
        vert2[8] = v2[5];                                                                                     
    }                                                                                                         
    else if(edges1[1]==edges2[0])                                                                             
    {                                                                                                         
                                                                                                              
        vert1[0] = v1[3];                                                                                     
        vert1[1] = v1[4];                                                                                     
        vert1[2] = v1[5];                                                                                     
        vert1[3] = v1[6];                                                                                     
        vert1[4] = v1[7];                                                                                     
        vert1[5] = v1[8];                                                                                     
        vert1[6] = v1[0];                                                                                     
        vert1[7] = v1[1];                                                                                     
        vert1[8] = v1[2];                                                                                     
        vert2[0] = v1[3];                                                                                     
        vert2[1] = v1[4];                                                                                     
        vert2[2] = v1[5];                                                                                     
        vert2[3] = v1[6];                                                                                     
        vert2[4] = v1[7];                                                                                     
        vert2[5] = v1[8];                                                                                     
        vert2[6] = v2[6];                                                                                     
        vert2[7] = v2[7];                                                                                     
        vert2[8] = v2[8];                                                                                     
                                                                                                              
    }
    else if(edges1[1]==edges2[1])                                                                             
    {                                                                                                         
        vert1[0] = v1[3];                                                                                     
        vert1[1] = v1[4];                                                                                     
        vert1[2] = v1[5];                                                                                     
        vert1[3] = v1[6];                                                                                     
        vert1[4] = v1[7];                                                                                     
        vert1[5] = v1[8];                                                                                     
        vert1[6] = v1[0];                                                                                     
        vert1[7] = v1[1];                                                                                     
        vert1[8] = v1[2];                                                                                     
        vert2[0] = v1[3];                                                                                     
        vert2[1] = v1[4];                                                                                     
        vert2[2] = v1[5];                                                                                     
        vert2[3] = v1[6];                                                                                     
        vert2[4] = v1[7];                                                                                     
        vert2[5] = v1[8];                                                                                     
        vert2[6] = v2[0];                                                                                     
        vert2[7] = v2[1];                                                                                     
        vert2[8] = v2[2];                                                                                     
                                                                                                              
    }                                                                                                         
    else if(edges1[1]==edges2[2])                                                                             
    {                                                                                                         
        vert1[0] = v1[3];                                                                                     
        vert1[1] = v1[4];                                                                                     
        vert1[2] = v1[5];                                                                                     
        vert1[3] = v1[6];                                                                                     
        vert1[4] = v1[7];                                                                                     
        vert1[5] = v1[8];                                                                                     
        vert1[6] = v1[0];                                                                                     
        vert1[7] = v1[1];                                                                                     
        vert1[8] = v1[2];                                                                                     
        vert2[0] = v1[3];                                                                                     
        vert2[1] = v1[4];                                                                                     
        vert2[2] = v1[5];                                                                                     
        vert2[3] = v1[6];                                                                                     
        vert2[4] = v1[7];                                                                                     
        vert2[5] = v1[8];                                                                                     
        vert2[6] = v2[3];                                                                                     
        vert2[7] = v2[4];                                                                                     
        vert2[8] = v2[5];                                                                                     
                                                                                                              
    } 
    else if(edges1[2]==edges2[0])                                                                             
    {                                                                                                         
        vert1[0] = v1[6];                                                                                     
        vert1[1] = v1[7];                                                                                     
        vert1[2] = v1[8];                                                                                     
        vert1[3] = v1[0];                                                                                     
        vert1[4] = v1[1];                                                                                     
        vert1[5] = v1[2];                                                                                     
        vert1[6] = v1[3];                                                                                     
        vert1[7] = v1[4];                                                                                     
        vert1[8] = v1[5];                                                                                     
        vert2[0] = v1[6];                                                                                     
        vert2[1] = v1[7];                                                                                     
        vert2[2] = v1[8];                                                                                     
        vert2[3] = v1[0];                                                                                     
        vert2[4] = v1[1];                                                                                     
        vert2[5] = v1[2];                                                                                     
        vert2[6] = v2[6];                                                                                     
        vert2[7] = v2[7];                                                                                     
        vert2[8] = v2[8];                                                                                     
                                                                                                              
    }                                                                                                         
    else if(edges1[2]==edges2[1])                                                                             
    {                                                                                                         
        vert1[0] = v1[6];                                                                                     
        vert1[1] = v1[7];                                                                                     
        vert1[2] = v1[8];                                                                                     
        vert1[3] = v1[0];                                                                                     
        vert1[4] = v1[1];                                                                                     
        vert1[5] = v1[2];                                                                                     
        vert1[6] = v1[3];                                                                                     
        vert1[7] = v1[4];                                                                                     
        vert1[8] = v1[5];                                                                                     
        vert2[0] = v1[6];                                                                                     
        vert2[1] = v1[7];                                                                                     
        vert2[2] = v1[8];                                                                                     
        vert2[3] = v1[0];                                                                                     
        vert2[4] = v1[1];                                                                                     
        vert2[5] = v1[2];                                                                                     
        vert2[6] = v2[0];                                                                                     
        vert2[7] = v2[1];                                                                                     
        vert2[8] = v2[2];                                                                                     
    }                                                                                                         
    else if(edges1[2]==edges2[2])                                                                             
    {                                                                                                         
        vert1[0] = v1[6];                                                                                     
        vert1[1] = v1[7];                                                                                     
        vert1[2] = v1[8];                                                                                     
        vert1[3] = v1[0];                                                                                     
        vert1[4] = v1[1];                                                                                     
        vert1[5] = v1[2];                                                                                     
        vert1[6] = v1[3];                                                                                     
        vert1[7] = v1[4];                                                                                     
        vert1[8] = v1[5];                                                                                     
        vert2[0] = v1[6];                                                                                     
        vert2[1] = v1[7];                                                                                     
        vert2[2] = v1[8];                                                                                     
        vert2[3] = v1[0]; 
        vert2[4] = v1[1];                                                                                     
        vert2[5] = v1[2];                                                                                     
        vert2[6] = v2[3];                                                                                     
        vert2[7] = v2[4];                                                                                     
        vert2[8] = v2[5];                                                                                     
    }                                                                                                         
     
    vert3[0] = vert1[0]-vert1[3];                                                                             
    vert3[1] = vert1[1]-vert1[4];                                                                             
    vert3[2] = vert1[2]-vert1[5];                                                                             
    vert4[0] = vert1[3]-vert2[6];                                                                             
    vert4[1] = vert1[4]-vert2[7];                                                                             
    vert4[2] = vert1[5]-vert2[8];                                                                             
    vert5[0] = vert1[6]-vert2[6];                                                                             
    vert5[1] = vert1[7]-vert2[7];                                                                             
    vert5[2] = vert1[8]-vert2[8];                                                                             
                                                                                                              
    for(i=0; i< sz2[0]; i++)                                                                                  
    {                                                                                                         
        for(j=0; j< sz2[1]; j++)                                                                              
        {                                                                                                     
            for(k=0; k< sz[0]; k++)                                                                           
            {                                                                                                 
                for(l=0; l< sz[1]; l++)                                                                       
                {   
                    y[0] = (1-chi[i])*eta3[j]+chi[i];
                    y[1] = chi[i]*(1-eta1[l]+eta1[l]*eta2[k]);
                    z[0] = -chi[i]*eta1[l];
                    z[1] = -chi[i]*eta1[l]*eta2[k];
                    x[0] = y[0]+z[0];
                    x[1] = y[1]+z[1];                                                                                        
                    point1[0] = -eta1[l]*vert3[0]-eta1[l]*eta2[k]*vert4[0]+(1-eta1[l]+eta1[l]*eta2[k])*vert5[0];
                    point1[1] = -eta1[l]*vert3[1]-eta1[l]*eta2[k]*vert4[1]+(1-eta1[l]+eta1[l]*eta2[k])*vert5[1];
                    point1[2] = -eta1[l]*vert3[2]-eta1[l]*eta2[k]*vert4[2]+(1-eta1[l]+eta1[l]*eta2[k])*vert5[2];
                    n1 = sqrt( (point1[0]*point1[0]+point1[1]*point1[1]+point1[2]*point1[2])*chi[i]*chi[i] ); 
                    n2 = sqrt( point1[0]*point1[0]+point1[1]*point1[1]+point1[2]*point1[2] );  
                    point3[0] = (1-x[0])*vert1[0]+(x[0]-x[1])*vert1[3]+x[1]*vert1[6];                         
                    point3[1] = (1-x[0])*vert1[1]+(x[0]-x[1])*vert1[4]+x[1]*vert1[7];                         
                    point3[2] = (1-x[0])*vert1[2]+(x[0]-x[1])*vert1[5]+x[1]*vert1[8];                         
                    point4[0] = (1-y[0])*vert2[0]+(y[0]-y[1])*vert2[3]+y[1]*vert2[6];                         
                    point4[1] = (1-y[0])*vert2[1]+(y[0]-y[1])*vert2[4]+y[1]*vert2[7];                         
                    point4[2] = (1-y[0])*vert2[2]+(y[0]-y[1])*vert2[5]+y[1]*vert2[8];                
                    ev_basis_function_maxwell(point3, phii, p01, A_t1, len1, NULL, 1, opt1, "RWG");            
                    ev_basis_function_maxwell(point4, phij, p02, A_t2, len2, NULL, 1, opt2, "RWG");  
                    res_j+= w1[l]*w2[k]*w3[j]*w4[i]*chi[i]*eta1[l]*(1-chi[i])*cexp((1*I)*kappa*n1)*dot_prod(3,0,0,phii,phij)/(n2*4.0*PI);
                 }                                                                                            
            }                                                                                                 
        }                                                                                                     
    }                                                                                                         
    res[0] += creal(res_j);                                                                                   
    res[1] += cimag(res_j);                                                                                   
    res_j = 0.0+0.0*I;                                                                                        
                                                                                                              
    for(i=0; i< sz2[0]; i++)                                                                                  
    {                                                                                                         
        for(j=0; j< sz2[1]; j++)                                                                              
        {                                                                                                     
            for(k=0; k< sz[0]; k++)                                                                           
            {                                                                                                 
                for(l=0; l< sz[1]; l++)                                                                       
                {                                                                                             
                    y[0] = (1-chi[i])*eta3[j]+chi[i]*(1-eta1[l]);                                                         
                    y[1] = chi[i]*(1-eta1[l]);                                                
                    z[0] = chi[i]*eta1[l];                                                                   
                    z[1] = chi[i]*eta1[l]*eta2[k];                                                           
                    x[0] = y[0]+z[0];                                                                         
                    x[1] = y[1]+z[1]; 
                    point1[0] = eta1[l]*vert3[0]+eta1[l]*eta2[k]*vert4[0]+(1-eta1[l])*vert5[0];               
                    point1[1] = eta1[l]*vert3[1]+eta1[l]*eta2[k]*vert4[1]+(1-eta1[l])*vert5[1];               
                    point1[2] = eta1[l]*vert3[2]+eta1[l]*eta2[k]*vert4[2]+(1-eta1[l])*vert5[2];               
                    n1 = sqrt( (point1[0]*point1[0]+point1[1]*point1[1]+point1[2]*point1[2])*chi[i]*chi[i] ); 
                    n2 = sqrt( point1[0]*point1[0]+point1[1]*point1[1]+point1[2]*point1[2] );
                    point3[0] = (1-x[0])*vert1[0]+(x[0]-x[1])*vert1[3]+x[1]*vert1[6];                         
                    point3[1] = (1-x[0])*vert1[1]+(x[0]-x[1])*vert1[4]+x[1]*vert1[7];                         
                    point3[2] = (1-x[0])*vert1[2]+(x[0]-x[1])*vert1[5]+x[1]*vert1[8];                         
                    point4[0] = (1-y[0])*vert2[0]+(y[0]-y[1])*vert2[3]+y[1]*vert2[6];                         
                    point4[1] = (1-y[0])*vert2[1]+(y[0]-y[1])*vert2[4]+y[1]*vert2[7];                         
                    point4[2] = (1-y[0])*vert2[2]+(y[0]-y[1])*vert2[5]+y[1]*vert2[8];         
                    ev_basis_function_maxwell(point3, phii, p01, A_t1, len1, NULL, 1, opt1, "RWG");            
                    ev_basis_function_maxwell(point4, phij, p02, A_t2, len2, NULL, 1, opt2, "RWG");           
                    res_j+= w1[l]*w2[k]*w3[j]*w4[i]*chi[i]*eta1[l]*(1-chi[i])*cexp((1*I)*kappa*n1)*dot_prod(3,0,0,phii,phij)/(n2*4.0*PI);
                 }                                                                                            
            }                                                                                                 
        }                                                                                                     
    } 
    res[0] += creal(res_j);                                                                                   
    res[1] += cimag(res_j);                                                                                   
    res_j = 0.0+0.0*I;                                                                                        
    for(i=0; i< sz2[0]; i++)                                                                                  
    {                                                                                                         
        for(j=0; j< sz2[1]; j++)                                                                              
        {                                                                                                     
            for(k=0; k< sz[0]; k++)                                                                           
            {                                                                                                 
                for(l=0; l< sz[1]; l++)                                                                       
                {         
                    y[0] = (1-chi[i])*eta3[j]+chi[i];                                              
                    y[1] = chi[i]*(1-eta1[l]);                                                                
                    z[0] = -chi[i]*eta1[l]*eta2[k];                                                                    
                    z[1] = chi[i]*eta1[l]*(1-eta2[k]);                                                            
                    x[0] = y[0]+z[0];                                                                         
                    x[1] = y[1]+z[1];                                                                                     
                    point1[0] = eta1[l]*eta2[k]*vert3[0]+eta1[l]*(1-eta2[k])*vert4[0]+(1-eta1[l])*vert5[0];   
                    point1[1] = eta1[l]*eta2[k]*vert3[1]+eta1[l]*(1-eta2[k])*vert4[1]+(1-eta1[l])*vert5[1];   
                    point1[2] = eta1[l]*eta2[k]*vert3[2]+eta1[l]*(1-eta2[k])*vert4[2]+(1-eta1[l])*vert5[2];   
                    n1 = sqrt( (point1[0]*point1[0]+point1[1]*point1[1]+point1[2]*point1[2])*chi[i]*chi[i] ); 
                    n2 = sqrt( point1[0]*point1[0]+point1[1]*point1[1]+point1[2]*point1[2] );   
                    point3[0] = (1-x[0])*vert1[0]+(x[0]-x[1])*vert1[3]+x[1]*vert1[6];                         
                    point3[1] = (1-x[0])*vert1[1]+(x[0]-x[1])*vert1[4]+x[1]*vert1[7];                         
                    point3[2] = (1-x[0])*vert1[2]+(x[0]-x[1])*vert1[5]+x[1]*vert1[8];                         
                    point4[0] = (1-y[0])*vert2[0]+(y[0]-y[1])*vert2[3]+y[1]*vert2[6];                         
                    point4[1] = (1-y[0])*vert2[1]+(y[0]-y[1])*vert2[4]+y[1]*vert2[7];                         
                    point4[2] = (1-y[0])*vert2[2]+(y[0]-y[1])*vert2[5]+y[1]*vert2[8];      
                    ev_basis_function_maxwell(point3, phii, p01, A_t1, len1, NULL, 1, opt1, "RWG");            
                    ev_basis_function_maxwell(point4, phij, p02, A_t2, len2, NULL, 1, opt2, "RWG");           
                    res_j+= w1[l]*w2[k]*w3[j]*w4[i]*chi[i]*eta1[l]*(1-chi[i])*cexp((1*I)*kappa*n1)*dot_prod(3,0,0,phii,phij)/(n2*4.0*PI);
                 }                                                                                            
            }                                                                                                 
        }                                                                                                     
    }                                                                                                         
    res[0] += creal(res_j);                                                                                   
    res[1] += cimag(res_j);                                                                                   
    res_j = 0.0+0.0*I;                                                                                        
    for(i=0; i< sz2[0]; i++)                                                                                  
    {                                                                                                         
        for(j=0; j< sz2[1]; j++)                                                                              
        {                                                                                                     
            for(k=0; k< sz[0]; k++)                                                                           
            {                                                                                                 
                for(l=0; l< sz[1]; l++)                                                                       
                {                                    
                    y[0] = (1-chi[i])*eta3[j]+chi[i]*(1-eta1[l]*eta2[k]);                                                         
                    y[1] = chi[i]*(1-eta1[l]*eta2[k]);                                                                
                    z[0] = chi[i]*eta1[l]*eta2[k];                                                           
                    z[1] = -chi[i]*eta1[l]*(1-eta2[k]);                                                        
                    x[0] = y[0]+z[0];                                                                         
                    x[1] = y[1]+z[1];                                                         
                    point1[0] = eta1[l]*eta2[k]*vert3[0]-eta1[l]*(1-eta2[k])*vert4[0]+(1-eta1[l])*vert5[0];   
                    point1[1] = eta1[l]*eta2[k]*vert3[1]-eta1[l]*(1-eta2[k])*vert4[1]+(1-eta1[l])*vert5[1];   
                    point1[2] = eta1[l]*eta2[k]*vert3[2]-eta1[l]*(1-eta2[k])*vert4[2]+(1-eta1[l])*vert5[2];   
                    n1 = sqrt( (point1[0]*point1[0]+point1[1]*point1[1]+point1[2]*point1[2])*chi[i]*chi[i] ); 
                    n2 = sqrt( point1[0]*point1[0]+point1[1]*point1[1]+point1[2]*point1[2] ); 
                    point3[0] = (1-x[0])*vert1[0]+(x[0]-x[1])*vert1[3]+x[1]*vert1[6];                         
                    point3[1] = (1-x[0])*vert1[1]+(x[0]-x[1])*vert1[4]+x[1]*vert1[7];                         
                    point3[2] = (1-x[0])*vert1[2]+(x[0]-x[1])*vert1[5]+x[1]*vert1[8];                         
                    point4[0] = (1-y[0])*vert2[0]+(y[0]-y[1])*vert2[3]+y[1]*vert2[6];                         
                    point4[1] = (1-y[0])*vert2[1]+(y[0]-y[1])*vert2[4]+y[1]*vert2[7];                         
                    point4[2] = (1-y[0])*vert2[2]+(y[0]-y[1])*vert2[5]+y[1]*vert2[8];        
                    ev_basis_function_maxwell(point3, phii, p01, A_t1, len1, NULL, 1, opt1, "RWG");            
                    ev_basis_function_maxwell(point4, phij, p02, A_t2, len2, NULL, 1, opt2, "RWG");           
                    res_j+= w1[l]*w2[k]*w3[j]*w4[i]*chi[i]*eta1[l]*(1-chi[i])*cexp((1*I)*kappa*n1)*dot_prod(3,0,0,phii,phij)/(n2*4.0*PI);
                 }                                                                                            
            }                                                                                                 
        }                                                                                                     
    }                                                                                                         
    res[0] += creal(res_j);                                                                                   
    res[1] += cimag(res_j);                                                                                   
    res_j = 0.0+0.0*I; 
    for(i=0; i< sz2[0]; i++)                                                                                  
    {                                                                                                         
        for(j=0; j< sz2[1]; j++)                                                                              
        {                                                                                                     
            for(k=0; k< sz[0]; k++)                                                                           
            {                                                                                                 
                for(l=0; l< sz[1]; l++)                                                                       
                {           
                    y[0] = (1-chi[i])*eta3[j]+chi[i];                                      
                    y[1] = chi[i];                                                        
                    z[0] = -chi[i]*eta1[l]*eta2[k];                                                            
                    z[1] = -chi[i]*eta1[l];                                                                                   
                    x[0] = y[0]+z[0];                                                                         
                    x[1] = y[1]+z[1];
                    point1[0] = -eta1[l]*eta2[k]*vert3[0]-eta1[l]*vert4[0]+vert5[0];                          
                    point1[1] = -eta1[l]*eta2[k]*vert3[1]-eta1[l]*vert4[1]+vert5[1];                          
                    point1[2] = -eta1[l]*eta2[k]*vert3[2]-eta1[l]*vert4[2]+vert5[2];                          
                    n1 = sqrt( (point1[0]*point1[0]+point1[1]*point1[1]+point1[2]*point1[2])*chi[i]*chi[i] ); 
                    n2 = sqrt( point1[0]*point1[0]+point1[1]*point1[1]+point1[2]*point1[2] );
                    point3[0] = (1-x[0])*vert1[0]+(x[0]-x[1])*vert1[3]+x[1]*vert1[6];                         
                    point3[1] = (1-x[0])*vert1[1]+(x[0]-x[1])*vert1[4]+x[1]*vert1[7];                         
                    point3[2] = (1-x[0])*vert1[2]+(x[0]-x[1])*vert1[5]+x[1]*vert1[8];                         
                    point4[0] = (1-y[0])*vert2[0]+(y[0]-y[1])*vert2[3]+y[1]*vert2[6];                         
                    point4[1] = (1-y[0])*vert2[1]+(y[0]-y[1])*vert2[4]+y[1]*vert2[7];                         
                    point4[2] = (1-y[0])*vert2[2]+(y[0]-y[1])*vert2[5]+y[1]*vert2[8];         
                    ev_basis_function_maxwell(point3, phii, p01, A_t1, len1, NULL, 1, opt1, "RWG");            
                    ev_basis_function_maxwell(point4, phij, p02, A_t2, len2, NULL, 1, opt2, "RWG");           
                    res_j+= w1[l]*w2[k]*w3[j]*w4[i]*chi[i]*eta1[l]*(1-chi[i])*cexp((1*I)*kappa*n1)*dot_prod(3,0,0,phii,phij)/(n2*4.0*PI);
                 }                                                                                            
            }                                                                                                 
        }                                                                                                     
    }                                                                                                         
    res[0] += creal(res_j);                                                                                   
    res[1] += cimag(res_j);                                                                                   
    res_j = 0.0+0.0*I;                                                                                        
    for(i=0; i< sz2[0]; i++)                                                                                  
    {                                                                                                         
        for(j=0; j< sz2[1]; j++)                                                                              
        {                                                                                                     
            for(k=0; k< sz[0]; k++)                                                                           
            {                                                                                                 
                for(l=0; l< sz[1]; l++)                                                                       
                {                    
                    y[0] = (1-chi[i])*eta3[j]+chi[i]*(1-eta1[l]*eta2[k]);                                                         
                    y[1] = chi[i]*(1-eta1[l]);                                                                            
                    z[0] = chi[i]*eta1[l]*eta2[k];                                                           
                    z[1] = chi[i]*eta1[l];                                                                   
                    x[0] = y[0]+z[0];                                                                         
                    x[1] = y[1]+z[1];                                                                          
                    point1[0] = eta1[l]*eta2[k]*vert3[0]+eta1[l]*vert4[0]+(1-eta1[l])*vert5[0];               
                    point1[1] = eta1[l]*eta2[k]*vert3[1]+eta1[l]*vert4[1]+(1-eta1[l])*vert5[1];               
                    point1[2] = eta1[l]*eta2[k]*vert3[2]+eta1[l]*vert4[2]+(1-eta1[l])*vert5[2];               
                    n1 = sqrt( (point1[0]*point1[0]+point1[1]*point1[1]+point1[2]*point1[2])*chi[i]*chi[i] ); 
                    n2 = sqrt( point1[0]*point1[0]+point1[1]*point1[1]+point1[2]*point1[2] );        
                    point3[0] = (1-x[0])*vert1[0]+(x[0]-x[1])*vert1[3]+x[1]*vert1[6];                         
                    point3[1] = (1-x[0])*vert1[1]+(x[0]-x[1])*vert1[4]+x[1]*vert1[7];                         
                    point3[2] = (1-x[0])*vert1[2]+(x[0]-x[1])*vert1[5]+x[1]*vert1[8];                         
                    point4[0] = (1-y[0])*vert2[0]+(y[0]-y[1])*vert2[3]+y[1]*vert2[6];                         
                    point4[1] = (1-y[0])*vert2[1]+(y[0]-y[1])*vert2[4]+y[1]*vert2[7];                         
                    point4[2] = (1-y[0])*vert2[2]+(y[0]-y[1])*vert2[5]+y[1]*vert2[8]; 
                    ev_basis_function_maxwell(point3, phii, p01, A_t1, len1, NULL, 1, opt1, "RWG");            
                    ev_basis_function_maxwell(point4, phij, p02, A_t2, len2, NULL, 1, opt2, "RWG");           
                    res_j+= w1[l]*w2[k]*w3[j]*w4[i]*chi[i]*eta1[l]*(1-chi[i])*cexp((1*I)*kappa*n1)*dot_prod(3,0,0,phii,phij)/(n2*4.0*PI);
                 }                                                                                            
            }                                                                                                 
        }                                                                                                     
    }                                                                                                          
    /*gsl_matrix_set(J1,0,0, vert1[3]-vert1[0]);                                                                
    gsl_matrix_set(J1,0,1, vert1[6]-vert1[3]);                                                                
    gsl_matrix_set(J1,1,0, vert1[4]-vert1[1]);                                                                
    gsl_matrix_set(J1,1,1, vert1[7]-vert1[4]);                                                                
    gsl_matrix_set(J1,2,0, vert1[5]-vert1[2]);                                                                
    gsl_matrix_set(J1,2,1, vert1[8]-vert1[5]);                                                                
    gsl_blas_dgemm(CblasTrans,CblasNoTrans, 1.0, J1, J1, 0.0, G1);                                            
    det1 = sqrt(fabs(gsl_matrix_get(G1,0,0)*gsl_matrix_get(G1,1,1)- gsl_matrix_get(G1,0,1)*gsl_matrix_get(G1,1,0)));
    gsl_matrix_set(J2,0,0, vert2[3]-vert2[0]);                                                                
    gsl_matrix_set(J2,0,1, vert2[6]-vert2[3]);                                                                
    gsl_matrix_set(J2,1,0, vert2[4]-vert2[1]);                                                                
    gsl_matrix_set(J2,1,1, vert2[7]-vert2[4]);                                                                
    gsl_matrix_set(J2,2,0, vert2[5]-vert2[2]);                                                                
    gsl_matrix_set(J2,2,1, vert2[8]-vert2[5]);                                                                
    gsl_blas_dgemm(CblasTrans,CblasNoTrans, 1.0, J2, J2, 0.0, G2);                                            
    det2 = sqrt(fabs(gsl_matrix_get(G2,0,0)*gsl_matrix_get(G2,1,1)- gsl_matrix_get(G2,0,1)*gsl_matrix_get(G2,1,0)));

   for(i=0; i< sz2[0]; i++)                                                                                  
    {                                                                                                         
        for(j=0; j< sz2[1]; j++)                                                                              
        {                                                                                                     
            for(k=0; k< sz[0]; k++)                                                                           
            {                                                                                                 
                for(l=0; l< sz[1]; l++)                                                                       
                {         
                    x[0] = 1;
                    z[0] = -eta1[l]*eta2[k];                                                                  
                    z[1] = eta1[l]*(1-eta2[k]);                                                               
                    z[2] = eta1[l]*eta3[j];           
                    x[1] = z[2];
                    y[0] = x[0]+z[0];
                    y[1] = z[1];
                    z2[0] = -eta2[k];                                                                         
                    z2[1] = (1-eta2[k]);                                                                      
                    z2[2] = eta3[j];                                                                          
                    point1[0] = (vert1[3]-vert1[0])*x[0]+(vert1[6]-vert1[3])*x[1]-((vert2[3]-vert2[0])*y[0]+(vert2[6]-vert2[3])*y[1]);
                    point1[1] = (vert1[4]-vert1[1])*x[1]+(vert1[7]-vert1[4])*x[1]-((vert2[4]-vert2[1])*y[0]+(vert2[7]-vert2[4])*y[1]);
                    point1[2] = (vert1[5]-vert1[2])*x[2]+(vert1[8]-vert1[5])*x[1]-((vert2[5]-vert2[2])*y[0]+(vert2[8]-vert2[5])*y[1]);   
                    point1[1] = (vert1[4]-vert1[1])*z[0]+(vert1[7]-vert1[4])*z[2]+(vert2[7]-vert2[4])*z[1];   
                    point1[2] = (vert1[5]-vert1[2])*z[0]+(vert1[8]-vert1[5])*z[2]+(vert2[8]-vert2[5])*z[1];   
                    point2[0] = (vert1[3]-vert1[0])*z2[0]+(vert1[6]-vert1[3])*z2[2]+(vert2[6]-vert2[3])*z2[1];
                    point2[1] = (vert1[4]-vert1[1])*z2[0]+(vert1[7]-vert1[4])*z2[2]+(vert2[7]-vert2[4])*z2[1];
                    point2[2] = (vert1[5]-vert1[2])*z2[0]+(vert1[8]-vert1[5])*z2[2]+(vert2[8]-vert2[5])*z2[1];
                    x[0] = chi[i]*x[0];                                                                       
                    x[1] = chi[i]*x[1];                                                                       
                    y[0] = chi[i]*y[0];                                                                       
                    y[1] = chi[i]*y[1];
                    point3[0] = vert1[0]+x[0]*(vert1[3]-vert1[0])+x[1]*(vert1[6]-vert1[3]);
                    point3[1] = vert1[1]+x[0]*(vert1[4]-vert1[1])+x[1]*(vert1[7]-vert1[4]); 
                    point3[2] = vert1[2]+x[0]*(vert1[5]-vert1[2])+x[1]*(vert1[8]-vert1[5]); 
                    point4[0] = vert2[0]+y[0]*(vert2[3]-vert2[0])+y[1]*(vert2[6]-vert2[3]);                   
                    point4[1] = vert2[1]+y[0]*(vert2[4]-vert2[1])+y[1]*(vert2[7]-vert2[4]);                   
                    point4[2] = vert2[2]+y[0]*(vert2[5]-vert2[2])+y[1]*(vert2[8]-vert2[5]);
                    ev_basis_function_maxwell(point3, phii, p01, A_t1, len1, NULL, 1, opt1, "RWG");            
                    ev_basis_function_maxwell(point4, phij, p02, A_t2, len2, NULL, 1, opt2, "RWG");
                    n1 = sqrt( (point1[0]*point1[0]+point1[1]*point1[1]+point1[2]*point1[2])*chi[i]*chi[i] ); 
                    n2 = sqrt( (point1[0]*point1[0]+point1[1]*point1[1]+point1[2]*point1[2]) ); 
                    //n2 = sqrt( point2[0]*point2[0]+point2[1]*point2[1]+point2[2]*point2[2] );                 
                    res_j+= w1[l]*w2[k]*w3[j]*w4[i]*chi[i]*chi[i]*eta1[l]*eta1[l]*cexp((1*I)*kappa*n1)*dot_prod(3,0,0,phii,phij)/(n2*4.0*PI);
                                                                                                              
                }                                                                                             
            }                                                                                                 
        }                                                                                                     
    }
    res[0] += creal(res_j);                                                                                   
    res[1] += cimag(res_j);                                                                                   
    res_j = 0.0+0.0*I;                                                                                        
    for(i=0; i< sz2[0]; i++)                                                                                  
    {                                                                                                         
        for(j=0; j< sz2[1]; j++)                                                                              
        {                                                                                                     
            for(k=0; k< sz[0]; k++)                                                                           
            {                                                                                                 
                for(l=0; l< sz[1]; l++)                                                                       
                {            
                    x[0] = 1;
                    z[0] = -eta1[l]*eta2[k]*eta3[j];                                                          
                    z[1] = eta2[k]*eta1[l]*(1-eta3[j]);                                                       
                    z[2] = eta1[l];                                                                           
                    z2[0] = -eta2[k]*eta3[j];                                                                 
                    z2[1] = eta2[k]*(1-eta3[j]);                                                              
                    z2[2] = 1; 
                    x[1] = z[2];                                                                              
                    y[0] = x[0]+z[0];                                                                         
                    y[1] = z[1];           
                    //x[0] = chi[i]*x[0];                                                                       
                    //x[1] = chi[i]*x[1];                                                                       
                    //y[0] = chi[i]*y[0];                                                                       
                    //y[1] = chi[i]*y[1];                                                                      
                    point1[0] = (vert1[3]-vert1[0])*z[0]+(vert1[6]-vert1[3])*z[2]+(vert2[6]-vert2[3])*z[1];   
                    point1[1] = (vert1[4]-vert1[1])*z[0]+(vert1[7]-vert1[4])*z[2]+(vert2[7]-vert2[4])*z[1];   
                    point1[2] = (vert1[5]-vert1[2])*z[0]+(vert1[8]-vert1[5])*z[2]+(vert2[8]-vert2[5])*z[1];   
                    point2[0] = (vert1[3]-vert1[0])*z2[0]+(vert1[6]-vert1[3])*z2[2]+(vert2[6]-vert2[3])*z2[1];
                    point2[1] = (vert1[4]-vert1[1])*z2[0]+(vert1[7]-vert1[4])*z2[2]+(vert2[7]-vert2[4])*z2[1];
                    point2[2] = (vert1[5]-vert1[2])*z2[0]+(vert1[8]-vert1[5])*z2[2]+(vert2[8]-vert2[5])*z2[1];
                    
                    point1[0] = (vert1[3]-vert1[0])*x[0]+(vert1[6]-vert1[3])*x[1]-((vert2[3]-vert2[0])*y[0]+(vert2[6]-vert2[3])*y[1]);
                    point1[1] = (vert1[4]-vert1[1])*x[1]+(vert1[7]-vert1[4])*x[1]-((vert2[4]-vert2[1])*y[0]+(vert2[7]-vert2[4])*y[1]);
                    point1[2] = (vert1[5]-vert1[2])*x[2]+(vert1[8]-vert1[5])*x[1]-((vert2[5]-vert2[2])*y[0]+(vert2[8]-vert2[5])*y[1]);
                    x[0] = chi[i]*x[0];                                                                       
                    x[1] = chi[i]*x[1];                                                                       
                    y[0] = chi[i]*y[0];                                                                       
                    y[1] = chi[i]*y[1];     
                    point3[0] = vert1[0]+x[0]*(vert1[3]-vert1[0])+x[1]*(vert1[6]-vert1[3]);                   
                    point3[1] = vert1[1]+x[0]*(vert1[4]-vert1[1])+x[1]*(vert1[7]-vert1[4]);                   
                    point3[2] = vert1[2]+x[0]*(vert1[5]-vert1[2])+x[1]*(vert1[8]-vert1[5]);                   
                    point4[0] = vert2[0]+y[0]*(vert2[3]-vert2[0])+y[1]*(vert2[6]-vert2[3]);                   
                    point4[1] = vert2[1]+y[0]*(vert2[4]-vert2[1])+y[1]*(vert2[7]-vert2[4]);                   
                    point4[2] = vert2[2]+y[0]*(vert2[5]-vert2[2])+y[1]*(vert2[8]-vert2[5]);
                    ev_basis_function_maxwell(point3, phii, p01, A_t1, len1, NULL, 1, opt1, "RWG");           
                    ev_basis_function_maxwell(point4, phij, p02, A_t2, len2, NULL, 1, opt2, "RWG");
                    n1 = sqrt( (point1[0]*point1[0]+point1[1]*point1[1]+point1[2]*point1[2])*chi[i]*chi[i] ); 
                    n2 = sqrt( (point1[0]*point1[0]+point1[1]*point1[1]+point1[2]*point1[2]) );
                    //n2 = sqrt( point2[0]*point2[0]+point2[1]*point2[1]+point2[2]*point2[2] );                 
                    res_j+= w1[l]*w2[k]*w3[j]*w4[i]*chi[i]*chi[i]*eta1[l]*eta1[l]*eta2[k]*cexp((1*I)*kappa*n1)*dot_prod(3,0,0,phii,phij)/(n2*4.0*PI);
                                                                                                              
                }                                                                                             
            }                                                                                                 
        }                                                                                                     
    }                                                                                                         
    res[0] += creal(res_j);                                                                                   
    res[1] += cimag(res_j);                                                                                   
    res_j = 0.0+0.0*I;     
    for(i=0; i< sz2[0]; i++)                                                                                  
    {                                                                                                         
        for(j=0; j< sz2[1]; j++)                                                                              
        {                                                                                                     
            for(k=0; k< sz[0]; k++)                                                                           
            {                                                                                                 
                for(l=0; l< sz[1]; l++)                                                                       
                {            
                    x[0] = 1-eta1[l]*eta2[k];
                    z[0] = eta1[l]*eta2[k];                                                                   
                    z[1] = eta2[k]*eta1[l]*eta3[j];                                                           
                    z[2] = eta1[l]*(1-eta2[k]);             
                    x[1] = z[2];                                                                              
                    y[0] = x[0]+z[0];                                                                         
                    y[1] = z[1];        
                    //x[0] = chi[i]*x[0];                                                                       
                    //x[1] = chi[i]*x[1];                                                                       
                    //y[0] = chi[i]*y[0];                                                                       
                    //y[1] = chi[i]*y[1];                                            
                    z2[0] = eta2[k];                                                                          
                    z2[1] = eta2[k]*eta3[j];                                                                  
                    z2[2] = (1-eta2[k]);                                                                      
                    point1[0] = (vert1[3]-vert1[0])*z[0]+(vert1[6]-vert1[3])*z[2]+(vert2[6]-vert2[3])*z[1];   
                    point1[1] = (vert1[4]-vert1[1])*z[0]+(vert1[7]-vert1[4])*z[2]+(vert2[7]-vert2[4])*z[1];   
                    point1[2] = (vert1[5]-vert1[2])*z[0]+(vert1[8]-vert1[5])*z[2]+(vert2[8]-vert2[5])*z[1];   
                    point2[0] = (vert1[3]-vert1[0])*z2[0]+(vert1[6]-vert1[3])*z2[2]+(vert2[6]-vert2[3])*z2[1];
                    point2[1] = (vert1[4]-vert1[1])*z2[0]+(vert1[7]-vert1[4])*z2[2]+(vert2[7]-vert2[4])*z2[1];
                    point2[2] = (vert1[5]-vert1[2])*z2[0]+(vert1[8]-vert1[5])*z2[2]+(vert2[8]-vert2[5])*z2[1];
                   
                    point1[0] = (vert1[3]-vert1[0])*x[0]+(vert1[6]-vert1[3])*x[1]-((vert2[3]-vert2[0])*y[0]+(vert2[6]-vert2[3])*y[1]);
                    point1[1] = (vert1[4]-vert1[1])*x[1]+(vert1[7]-vert1[4])*x[1]-((vert2[4]-vert2[1])*y[0]+(vert2[7]-vert2[4])*y[1]);
                    point1[2] = (vert1[5]-vert1[2])*x[2]+(vert1[8]-vert1[5])*x[1]-((vert2[5]-vert2[2])*y[0]+(vert2[8]-vert2[5])*y[1]); 
                    x[0] = chi[i]*x[0];                                                                       
                    x[1] = chi[i]*x[1];                                                                       
                    y[0] = chi[i]*y[0];                                                                       
                    y[1] = chi[i]*y[1];     
                    point3[0] = vert1[0]+x[0]*(vert1[3]-vert1[0])+x[1]*(vert1[6]-vert1[3]);                   
                    point3[1] = vert1[1]+x[0]*(vert1[4]-vert1[1])+x[1]*(vert1[7]-vert1[4]);                   
                    point3[2] = vert1[2]+x[0]*(vert1[5]-vert1[2])+x[1]*(vert1[8]-vert1[5]);                   
                    point4[0] = vert2[0]+y[0]*(vert2[3]-vert2[0])+y[1]*(vert2[6]-vert2[3]);                   
                    point4[1] = vert2[1]+y[0]*(vert2[4]-vert2[1])+y[1]*(vert2[7]-vert2[4]);                   
                    point4[2] = vert2[2]+y[0]*(vert2[5]-vert2[2])+y[1]*(vert2[8]-vert2[5]);
                    ev_basis_function_maxwell(point3, phii, p01, A_t1, len1, NULL, 1, opt1, "RWG");           
                    ev_basis_function_maxwell(point4, phij, p02, A_t2, len2, NULL, 1, opt2, "RWG");
                    n1 = sqrt( (point1[0]*point1[0]+point1[1]*point1[1]+point1[2]*point1[2])*chi[i]*chi[i] ); 
                    //n2 = sqrt( point2[0]*point2[0]+point2[1]*point2[1]+point2[2]*point2[2] );                 
                    n2 = sqrt( (point1[0]*point1[0]+point1[1]*point1[1]+point1[2]*point1[2]) );
                    res_j+= w1[l]*w2[k]*w3[j]*w4[i]*chi[i]*chi[i]*eta1[l]*eta1[l]*eta2[k]*cexp((1*I)*kappa*n1)*dot_prod(3,0,0,phii,phij)/(n2*4.0*PI);
                                                                                                              
                }                                                                                             
            }                                                                                                 
        }                                                                                                     
    }                                                                                                         
    res[0] += creal(res_j);                                                                                   
    res[1] += cimag(res_j);                                                                                   
    res_j = 0.0+0.0*I; 
    for(i=0; i< sz2[0]; i++)                                                                                  
    {                                                                                                         
        for(j=0; j< sz2[1]; j++)                                                                              
        {                                                                                                     
            for(k=0; k< sz[0]; k++)                                                                           
            {                                                                                                 
                for(l=0; l< sz[1]; l++)                                                                       
                {                            
                    x[0] = 1-eta1[l]*eta2[k]*eta3[j]; 
                    z[0] = eta1[l]*eta2[k]*eta3[j];                                                           
                    z[1] = eta1[l];                                                                           
                    z[2] = eta1[l]*eta2[k]*(1-eta3[j]);    
                    x[1] = z[2];                                                                              
                    y[0] = x[0]+z[0];                                                                         
                    y[1] = z[1];        
                    //x[0] = chi[i]*x[0];                                                                       
                    //x[1] = chi[i]*x[1];                                                                       
                    //y[0] = chi[i]*y[0];                                                                       
                    //y[1] = chi[i]*y[1];                                             
                    z2[0] = eta2[k]*eta3[j];                                                                  
                    z2[1] = 1;                                                                                
                    z2[2] = eta2[k]*(1-eta3[j]);                                                              
                    point1[0] = (vert1[3]-vert1[0])*z[0]+(vert1[6]-vert1[3])*z[2]+(vert2[6]-vert2[3])*z[1];   
                    point1[1] = (vert1[4]-vert1[1])*z[0]+(vert1[7]-vert1[4])*z[2]+(vert2[7]-vert2[4])*z[1];   
                    point1[2] = (vert1[5]-vert1[2])*z[0]+(vert1[8]-vert1[5])*z[2]+(vert2[8]-vert2[5])*z[1];   
                    point2[0] = (vert1[3]-vert1[0])*z2[0]+(vert1[6]-vert1[3])*z2[2]+(vert2[6]-vert2[3])*z2[1];
                    point2[1] = (vert1[4]-vert1[1])*z2[0]+(vert1[7]-vert1[4])*z2[2]+(vert2[7]-vert2[4])*z2[1];
                    point2[2] = (vert1[5]-vert1[2])*z2[0]+(vert1[8]-vert1[5])*z2[2]+(vert2[8]-vert2[5])*z2[1];
                    
                    point1[0] = (vert1[3]-vert1[0])*x[0]+(vert1[6]-vert1[3])*x[1]-((vert2[3]-vert2[0])*y[0]+(vert2[6]-vert2[3])*y[1]);
                    point1[1] = (vert1[4]-vert1[1])*x[1]+(vert1[7]-vert1[4])*x[1]-((vert2[4]-vert2[1])*y[0]+(vert2[7]-vert2[4])*y[1]);
                    point1[2] = (vert1[5]-vert1[2])*x[2]+(vert1[8]-vert1[5])*x[1]-((vert2[5]-vert2[2])*y[0]+(vert2[8]-vert2[5])*y[1]);
                    x[0] = chi[i]*x[0];                                                                       
                    x[1] = chi[i]*x[1];                                                                       
                    y[0] = chi[i]*y[0];                                                                       
                    y[1] = chi[i]*y[1];     
                    point3[0] = vert1[0]+x[0]*(vert1[3]-vert1[0])+x[1]*(vert1[6]-vert1[3]);                   
                    point3[1] = vert1[1]+x[0]*(vert1[4]-vert1[1])+x[1]*(vert1[7]-vert1[4]);                   
                    point3[2] = vert1[2]+x[0]*(vert1[5]-vert1[2])+x[1]*(vert1[8]-vert1[5]);                   
                    point4[0] = vert2[0]+y[0]*(vert2[3]-vert2[0])+y[1]*(vert2[6]-vert2[3]);                   
                    point4[1] = vert2[1]+y[0]*(vert2[4]-vert2[1])+y[1]*(vert2[7]-vert2[4]);                   
                    point4[2] = vert2[2]+y[0]*(vert2[5]-vert2[2])+y[1]*(vert2[8]-vert2[5]);
                    ev_basis_function_maxwell(point3, phii, p01, A_t1, len1, NULL, 1, opt1, "RWG");           
                    ev_basis_function_maxwell(point4, phij, p02, A_t2, len2, NULL, 1, opt2, "RWG");
                    n1 = sqrt( (point1[0]*point1[0]+point1[1]*point1[1]+point1[2]*point1[2])*chi[i]*chi[i] ); 
                    n2 = sqrt( (point1[0]*point1[0]+point1[1]*point1[1]+point1[2]*point1[2]) ); 
                    //n2 = sqrt( point2[0]*point2[0]+point2[1]*point2[1]+point2[2]*point2[2] );                 
                    res_j+= w1[l]*w2[k]*w3[j]*w4[i]*chi[i]*chi[i]*eta1[l]*eta1[l]*eta2[k]*cexp((1*I)*kappa*n1)*dot_prod(3,0,0,phii,phij)/(n2*4.0*PI);
                                                                                                              
                }                                                                                             
            }                                                                                                 
        }                                                                                                     
    }                                                                                                         
    res[0] += creal(res_j);                                                                                   
    res[1] += cimag(res_j);                                                                                   
    res_j = 0.0+0.0*I; 
    for(i=0; i< sz2[0]; i++)                                                                                  
    {                                                                                                         
        for(j=0; j< sz2[1]; j++)                                                                              
        {                                                                                                     
            for(k=0; k< sz[0]; k++)                                                                           
            {                                                                                                 
                for(l=0; l< sz[1]; l++)                                                                       
                {                                   
                    x[0] = 1-eta1[l]*eta2[k]*eta3[j]; 
                    z[0] = eta1[l]*eta2[k]*eta3[j];                                                           
                    z[1] = eta1[l]*eta2[k];                                                                   
                    z[2] = eta1[l]*(1-eta3[j]*eta2[k]);
                    x[1] = z[2];                                                                              
                    y[0] = x[0]+z[0];                                                                         
                    y[1] = z[1];        
                    //x[0] = chi[i]*x[0];                                                                       
                    //x[1] = chi[i]*x[1];                                                                       
                    //y[0] = chi[i]*y[0];                                                                       
                    //y[1] = chi[i]*y[1];                                                 
                    z2[0] = eta2[k]*eta3[j];                                                                  
                    z2[1] = eta2[k];                                                                          
                    z2[2] = 1-eta3[j]*eta2[k];                                                                
                    point1[0] = (vert1[3]-vert1[0])*z[0]+(vert1[6]-vert1[3])*z[2]+(vert2[6]-vert2[3])*z[1];   
                    point1[1] = (vert1[4]-vert1[1])*z[0]+(vert1[7]-vert1[4])*z[2]+(vert2[7]-vert2[4])*z[1];   
                    point1[2] = (vert1[5]-vert1[2])*z[0]+(vert1[8]-vert1[5])*z[2]+(vert2[8]-vert2[5])*z[1];   
                    point2[0] = (vert1[3]-vert1[0])*z2[0]+(vert1[6]-vert1[3])*z2[2]+(vert2[6]-vert2[3])*z2[1];
                    point2[1] = (vert1[4]-vert1[1])*z2[0]+(vert1[7]-vert1[4])*z2[2]+(vert2[7]-vert2[4])*z2[1];
                    point2[2] = (vert1[5]-vert1[2])*z2[0]+(vert1[8]-vert1[5])*z2[2]+(vert2[8]-vert2[5])*z2[1];
                    point1[0] = (vert1[3]-vert1[0])*x[0]+(vert1[6]-vert1[3])*x[1]-((vert2[3]-vert2[0])*y[0]+(vert2[6]-vert2[3])*y[1]);
                    point1[1] = (vert1[4]-vert1[1])*x[1]+(vert1[7]-vert1[4])*x[1]-((vert2[4]-vert2[1])*y[0]+(vert2[7]-vert2[4])*y[1]);
                    point1[2] = (vert1[5]-vert1[2])*x[2]+(vert1[8]-vert1[5])*x[1]-((vert2[5]-vert2[2])*y[0]+(vert2[8]-vert2[5])*y[1]);
                    x[0] = chi[i]*x[0];                                                                       
                    x[1] = chi[i]*x[1];                                                                       
                    y[0] = chi[i]*y[0];                                                                       
                    y[1] = chi[i]*y[1];     
                    point3[0] = vert1[0]+x[0]*(vert1[3]-vert1[0])+x[1]*(vert1[6]-vert1[3]);                   
                    point3[1] = vert1[1]+x[0]*(vert1[4]-vert1[1])+x[1]*(vert1[7]-vert1[4]);                   
                    point3[2] = vert1[2]+x[0]*(vert1[5]-vert1[2])+x[1]*(vert1[8]-vert1[5]);                   
                    point4[0] = vert2[0]+y[0]*(vert2[3]-vert2[0])+y[1]*(vert2[6]-vert2[3]);                   
                    point4[1] = vert2[1]+y[0]*(vert2[4]-vert2[1])+y[1]*(vert2[7]-vert2[4]);                   
                    point4[2] = vert2[2]+y[0]*(vert2[5]-vert2[2])+y[1]*(vert2[8]-vert2[5]);
                    ev_basis_function_maxwell(point3, phii, p01, A_t1, len1, NULL, 1, opt1, "RWG");           
                    ev_basis_function_maxwell(point4, phij, p02, A_t2, len2, NULL, 1, opt2, "RWG");
                    n1 = sqrt( (point1[0]*point1[0]+point1[1]*point1[1]+point1[2]*point1[2])*chi[i]*chi[i] ); 
                    n2 = sqrt( point1[0]*point1[0]+point1[1]*point1[1]+point1[2]*point1[2] );                 
                    //n2 = sqrt( point2[0]*point2[0]+point2[1]*point2[1]+point2[2]*point2[2] );
                    res_j+= w1[l]*w2[k]*w3[j]*w4[i]*chi[i]*chi[i]*eta1[l]*eta1[l]*eta2[k]*cexp((1*I)*kappa*n1)*dot_prod(3,0,0,phii,phij)/(n2*4.0*PI);
                                                                                                              
                }                                                                                             
            }                                                                                                 
        }                                                  
    }                                                      
    */
    res[0] += creal(res_j);                                                                                   
    res[1] += cimag(res_j);                                                                                   
    res[0] = res[0]*4*A_t1*A_t2;                                                                              
    res[1] = res[1]*4*A_t1*A_t2;
                                                                                                              
}                   
void singular_maxwell_hypersingular_ce(double kappa, int *order, double *v1, double *v2,  double *p01, double *p02, double *res, gsl_vector *n, double len1, double len2, double A_t1, double A_t2, int *edges1, int *edges2, int opt1, int opt2)
{                                                                                                             
    double eta1[4], eta2[4], eta3[4], chi[4], w1[4], w2[4], w3[4], w4[4];                                     
    double vert1[9] = {v1[0], v1[1], v1[2], v1[3], v1[4], v1[5], v1[6], v1[7], v1[8]};
    double vert2[9] = {v2[0], v2[1], v2[2], v2[3], v2[4], v2[5], v2[6], v2[7], v2[8]}; 
    double vert3[3], vert4[3], vert5[3];
    int sz[2], sz2[2];                                                                                        
    double  point1[3], point2[3];                                                     
    int i,j,k,l;                                                                                              
    double n1, n2;
    double complex res_j = 0.0+0.0*I;                                                                         
    res[0] = 0.0;                                                                                             
    res[1] = 0.0;                                                                                             
    quad_rules_ab( 0.0, 1.0, eta1, eta2, w1, w2, sz, order[2]);                                               
    quad_rules_ab( 0.0, 1.0, eta3, chi, w3, w4, sz2, order[2]);                                               
    double phii[sz[0]], phij[sz[0]];                                                                          
    divergence(A_t1,len1, NULL,"RWG", opt1, phii, 1);                                                          
    divergence(A_t2,len2, NULL,"RWG", opt2, phij, 1); 
    
    if(edges1[0]==edges2[0])
    {
        vert1[0] = v1[0];
        vert1[1] = v1[1];
        vert1[2] = v1[2];
        vert1[3] = v1[3];
        vert1[4] = v1[4];
        vert1[5] = v1[5];
        vert1[6] = v1[6];                                                                                
        vert1[7] = v1[7];                                                                                
        vert1[8] = v1[8];
        vert2[0] = v1[0];                                                                                     
        vert2[1] = v1[1];                                                                                     
        vert2[2] = v1[2];                                                                                     
        vert2[3] = v1[3];                                                                                     
        vert2[4] = v1[4];                                                                                     
        vert2[5] = v1[5];                                                                                     
        vert2[6] = v2[6];                                                                                     
        vert2[7] = v2[7];                                                                                     
        vert2[8] = v2[8]; 

    }
    else if(edges1[0]==edges2[1])                                                                            
    {                
        vert1[0] = v1[0];                                                                                     
        vert1[1] = v1[1];                                                                                     
        vert1[2] = v1[2];                                                                                     
        vert1[3] = v1[3];                                                                                     
        vert1[4] = v1[4];                                                                                     
        vert1[5] = v1[5];                                                                                     
        vert1[6] = v1[6];                                                                                     
        vert1[7] = v1[7];                                                                                     
        vert1[8] = v1[8];                                                                                     
        vert2[0] = v1[0];                                                                                     
        vert2[1] = v1[1];                                                                                     
        vert2[2] = v1[2];                                                                                     
        vert2[3] = v1[3];                                                                                     
        vert2[4] = v1[4];                                                                                     
        vert2[5] = v1[5];                                                                                     
        vert2[6] = v2[0];                                                                                     
        vert2[7] = v2[1];                                                                                     
        vert2[8] = v2[2];  
        
    } 
    else if(edges1[0]==edges2[2])                                                                            
    {                             
        vert1[0] = v1[0];                                                                                     
        vert1[1] = v1[1];                                                                                     
        vert1[2] = v1[2];                                                                                     
        vert1[3] = v1[3];                                                                                     
        vert1[4] = v1[4];                                                                                     
        vert1[5] = v1[5];                                                                                     
        vert1[6] = v1[6];                                                                                     
        vert1[7] = v1[7];                                                                                     
        vert1[8] = v1[8];             
        vert2[0] = v1[0];                                                                                     
        vert2[1] = v1[1];                                                                                     
        vert2[2] = v1[2];                                                                                     
        vert2[3] = v1[3];                                                                                     
        vert2[4] = v1[4];                                                                                     
        vert2[5] = v1[5];                                                                                     
        vert2[6] = v2[3];                                                                                     
        vert2[7] = v2[4];                                                                                     
        vert2[8] = v2[5];

    } 
    else if(edges1[1]==edges2[0])                                                                       
    {                  

        vert1[0] = v1[3];                                                                                     
        vert1[1] = v1[4];                                                                                     
        vert1[2] = v1[5];                                                                                     
        vert1[3] = v1[6];                                                                                     
        vert1[4] = v1[7];                                                                                     
        vert1[5] = v1[8];                                                                                     
        vert1[6] = v1[0];                                                                                     
        vert1[7] = v1[1];                                                                                     
        vert1[8] = v1[2];                                                                                     
        vert2[0] = v1[3];                                                                                     
        vert2[1] = v1[4];                                                                                     
        vert2[2] = v1[5];                                                                                     
        vert2[3] = v1[6];                                                                                     
        vert2[4] = v1[7];                                                                                     
        vert2[5] = v1[8];                                                                                     
        vert2[6] = v2[6];                                                                                     
        vert2[7] = v2[7];                                                                                     
        vert2[8] = v2[8]; 

                                                                                       
    }
    else if(edges1[1]==edges2[1])                                                                       
    {                  
        vert1[0] = v1[3];                                                                                     
        vert1[1] = v1[4];                                                                                     
        vert1[2] = v1[5];                                                                                     
        vert1[3] = v1[6];                                                                                     
        vert1[4] = v1[7];                                                                                     
        vert1[5] = v1[8];                                                                                     
        vert1[6] = v1[0];                                                                                     
        vert1[7] = v1[1];                                                                                     
        vert1[8] = v1[2];                                                                                     
        vert2[0] = v1[3];                                                                                     
        vert2[1] = v1[4];                                                                                     
        vert2[2] = v1[5];                                                                                     
        vert2[3] = v1[6];                                                                                     
        vert2[4] = v1[7];                                                                                     
        vert2[5] = v1[8];                                                                                     
        vert2[6] = v2[0];                                                                                     
        vert2[7] = v2[1];                                                                                     
        vert2[8] = v2[2];                                                                                       

    }
    else if(edges1[1]==edges2[2])                                                                       
    {                  
        vert1[0] = v1[3];                                                                                     
        vert1[1] = v1[4];                                                                                     
        vert1[2] = v1[5];                                                                                     
        vert1[3] = v1[6];                                                                                     
        vert1[4] = v1[7];                                                                                     
        vert1[5] = v1[8];                                                                                     
        vert1[6] = v1[0];                                                                                     
        vert1[7] = v1[1];                                                                                     
        vert1[8] = v1[2];                                                                                     
        vert2[0] = v1[3];                                                                                     
        vert2[1] = v1[4];                                                                                     
        vert2[2] = v1[5];                                                                                     
        vert2[3] = v1[6];                                                                                     
        vert2[4] = v1[7];                                                                                     
        vert2[5] = v1[8];                                                                                     
        vert2[6] = v2[3];                                                                                     
        vert2[7] = v2[4];                                                                                     
        vert2[8] = v2[5];                                                                                       

    }
    else if(edges1[2]==edges2[0])                                                                       
    {   
        vert1[0] = v1[6];                                                                                     
        vert1[1] = v1[7];                                                                                     
        vert1[2] = v1[8];                                                                                     
        vert1[3] = v1[0];                                                                                     
        vert1[4] = v1[1];                                                                                     
        vert1[5] = v1[2];                                                                                     
        vert1[6] = v1[3];                                                                                     
        vert1[7] = v1[4];                                                                                     
        vert1[8] = v1[5];                                                                                     
        vert2[0] = v1[6];                                                                                     
        vert2[1] = v1[7];                                                                                     
        vert2[2] = v1[8];                                                                                     
        vert2[3] = v1[0];                                                                                     
        vert2[4] = v1[1];                                                                                     
        vert2[5] = v1[2];                                                                                     
        vert2[6] = v2[6];                                                                                     
        vert2[7] = v2[7];                                                                                     
        vert2[8] = v2[8];                                                                                                      

    }
    else if(edges1[2]==edges2[1])                                                                       
    {  
        vert1[0] = v1[6];                                                                                     
        vert1[1] = v1[7];                                                                                     
        vert1[2] = v1[8];                                                                                     
        vert1[3] = v1[0];                                                                                     
        vert1[4] = v1[1];                                                                                     
        vert1[5] = v1[2];                                                                                     
        vert1[6] = v1[3];                                                                                     
        vert1[7] = v1[4];                                                                                     
        vert1[8] = v1[5];                                                                                     
        vert2[0] = v1[6];                                                                                     
        vert2[1] = v1[7];                                                                                     
        vert2[2] = v1[8];                                                                                     
        vert2[3] = v1[0];                                                                                     
        vert2[4] = v1[1];                                                                                     
        vert2[5] = v1[2];                                                                                     
        vert2[6] = v2[0];                                                                                     
        vert2[7] = v2[1];                                                                                     
        vert2[8] = v2[2];

    }
    else if(edges1[2]==edges2[2])                                                                       
    {                  
        vert1[0] = v1[6];                                                                                     
        vert1[1] = v1[7];                                                                                     
        vert1[2] = v1[8];                                                                                     
        vert1[3] = v1[0];                                                                                     
        vert1[4] = v1[1];                                                                                     
        vert1[5] = v1[2];                                                                                     
        vert1[6] = v1[3];                                                                                     
        vert1[7] = v1[4];                                                                                     
        vert1[8] = v1[5];                                                                                     
        vert2[0] = v1[6];                                                                                     
        vert2[1] = v1[7];                                                                                     
        vert2[2] = v1[8];                                                                                     
        vert2[3] = v1[0];                                                                                     
        vert2[4] = v1[1];                                                                                     
        vert2[5] = v1[2];                                                                                     
        vert2[6] = v2[3];                                                                                     
        vert2[7] = v2[4];                                                                                     
        vert2[8] = v2[5];

    }
    
    vert3[0] = vert1[0]-vert1[3];
    vert3[1] = vert1[1]-vert1[4];
    vert3[2] = vert1[2]-vert1[5];
    vert4[0] = vert1[3]-vert2[6];                                                                             
    vert4[1] = vert1[4]-vert2[7];                                                                             
    vert4[2] = vert1[5]-vert2[8];
    vert5[0] = vert1[6]-vert2[6];                                                                             
    vert5[1] = vert1[7]-vert2[7];                                                                             
    vert5[2] = vert1[8]-vert2[8];

    for(i=0; i< sz2[0]; i++)                                                                                  
    {                                                                                                         
        for(j=0; j< sz2[1]; j++)                                                                              
        {                                                                                                     
            for(k=0; k< sz[0]; k++)                                                                           
            {                                                                                                 
                for(l=0; l< sz[1]; l++)                                                                       
                {
                    point1[0] = -eta1[l]*vert3[0]-eta1[l]*eta2[k]*vert4[0]+(1-eta1[l]+eta1[l]*eta2[k])*vert5[0];
                    point1[1] = -eta1[l]*vert3[1]-eta1[l]*eta2[k]*vert4[1]+(1-eta1[l]+eta1[l]*eta2[k])*vert5[1];
                    point1[2] = -eta1[l]*vert3[2]-eta1[l]*eta2[k]*vert4[2]+(1-eta1[l]+eta1[l]*eta2[k])*vert5[2];
                    n1 = sqrt( (point1[0]*point1[0]+point1[1]*point1[1]+point1[2]*point1[2])*chi[i]*chi[i] );
                    n2 = sqrt( point1[0]*point1[0]+point1[1]*point1[1]+point1[2]*point1[2] );                 
                    res_j+= w1[l]*w2[k]*w3[j]*w4[i]*chi[i]*eta1[l]*(1-chi[i])*cexp((1*I)*kappa*n1)*phii[0]*phij[0]/(n2*4.0*PI);
                 }  
            }
        }
    }
    res[0] += creal(res_j);                                                                                   
    res[1] += cimag(res_j);                                                                                   
    res_j = 0.0+0.0*I; 

    for(i=0; i< sz2[0]; i++)                                                                                  
    {                                                                                                         
        for(j=0; j< sz2[1]; j++)                                                                              
        {                                                                                                     
            for(k=0; k< sz[0]; k++)                                                                           
            {                                                                                                 
                for(l=0; l< sz[1]; l++)                                                                       
                {                                                                                             
                    point1[0] = eta1[l]*vert3[0]+eta1[l]*eta2[k]*vert4[0]+(1-eta1[l])*vert5[0];
                    point1[1] = eta1[l]*vert3[1]+eta1[l]*eta2[k]*vert4[1]+(1-eta1[l])*vert5[1];
                    point1[2] = eta1[l]*vert3[2]+eta1[l]*eta2[k]*vert4[2]+(1-eta1[l])*vert5[2];
                    n1 = sqrt( (point1[0]*point1[0]+point1[1]*point1[1]+point1[2]*point1[2])*chi[i]*chi[i] ); 
                    n2 = sqrt( point1[0]*point1[0]+point1[1]*point1[1]+point1[2]*point1[2] );                 
                    res_j+= w1[l]*w2[k]*w3[j]*w4[i]*chi[i]*eta1[l]*(1-chi[i])*cexp((1*I)*kappa*n1)*phii[0]*phij[0]/(n2*4.0*PI);
                 }                                                                                            
            }                                                                                                 
        }                                                                                                     
    } 
    res[0] += creal(res_j);                                                                                   
    res[1] += cimag(res_j);                                                                                   
    res_j = 0.0+0.0*I; 
    for(i=0; i< sz2[0]; i++)                                                                                  
    {                                                                                                         
        for(j=0; j< sz2[1]; j++)                                                                              
        {                                                                                                     
            for(k=0; k< sz[0]; k++)                                                                           
            {                                                                                                 
                for(l=0; l< sz[1]; l++)                                                                       
                {                                                                                             
                    point1[0] = eta1[l]*eta2[k]*vert3[0]+eta1[l]*(1-eta2[k])*vert4[0]+(1-eta1[l])*vert5[0];               
                    point1[1] = eta1[l]*eta2[k]*vert3[1]+eta1[l]*(1-eta2[k])*vert4[1]+(1-eta1[l])*vert5[1];               
                    point1[2] = eta1[l]*eta2[k]*vert3[2]+eta1[l]*(1-eta2[k])*vert4[2]+(1-eta1[l])*vert5[2];               
                    n1 = sqrt( (point1[0]*point1[0]+point1[1]*point1[1]+point1[2]*point1[2])*chi[i]*chi[i] ); 
                    n2 = sqrt( point1[0]*point1[0]+point1[1]*point1[1]+point1[2]*point1[2] );                 
                    res_j+= w1[l]*w2[k]*w3[j]*w4[i]*chi[i]*eta1[l]*(1-chi[i])*cexp((1*I)*kappa*n1)*phii[0]*phij[0]/(n2*4.0*PI);
                 }                                                                                            
            }                                                                                                 
        }                                                                                                     
    } 
    res[0] += creal(res_j);                                                                                   
    res[1] += cimag(res_j);                                                                                   
    res_j = 0.0+0.0*I; 
    for(i=0; i< sz2[0]; i++)                                                                                  
    {                                                                                                         
        for(j=0; j< sz2[1]; j++)                                                                              
        {                                                                                                     
            for(k=0; k< sz[0]; k++)                                                                           
            {                                                                                                 
                for(l=0; l< sz[1]; l++)                                                                       
                {                                                                                             
                    point1[0] = eta1[l]*eta2[k]*vert3[0]-eta1[l]*(1-eta2[k])*vert4[0]+(1-eta1[l])*vert5[0];   
                    point1[1] = eta1[l]*eta2[k]*vert3[1]-eta1[l]*(1-eta2[k])*vert4[1]+(1-eta1[l])*vert5[1];   
                    point1[2] = eta1[l]*eta2[k]*vert3[2]-eta1[l]*(1-eta2[k])*vert4[2]+(1-eta1[l])*vert5[2];   
                    n1 = sqrt( (point1[0]*point1[0]+point1[1]*point1[1]+point1[2]*point1[2])*chi[i]*chi[i] ); 
                    n2 = sqrt( point1[0]*point1[0]+point1[1]*point1[1]+point1[2]*point1[2] );                 
                    res_j+= w1[l]*w2[k]*w3[j]*w4[i]*chi[i]*eta1[l]*(1-chi[i])*cexp((1*I)*kappa*n1)*phii[0]*phij[0]/(n2*4.0*PI);
                 }                                                                                            
            }                                                                                                 
        }                                                                                                     
    }
    res[0] += creal(res_j);                                                                                   
    res[1] += cimag(res_j);                                                                                   
    res_j = 0.0+0.0*I; 
    for(i=0; i< sz2[0]; i++)                                                                                  
    {                                                                                                         
        for(j=0; j< sz2[1]; j++)                                                                              
        {                                                                                                     
            for(k=0; k< sz[0]; k++)                                                                           
            {                                                                                                 
                for(l=0; l< sz[1]; l++)                                                                       
                {                                                                                             
                    point1[0] = -eta1[l]*eta2[k]*vert3[0]-eta1[l]*vert4[0]+vert5[0];   
                    point1[1] = -eta1[l]*eta2[k]*vert3[1]-eta1[l]*vert4[1]+vert5[1];   
                    point1[2] = -eta1[l]*eta2[k]*vert3[2]-eta1[l]*vert4[2]+vert5[2];   
                    n1 = sqrt( (point1[0]*point1[0]+point1[1]*point1[1]+point1[2]*point1[2])*chi[i]*chi[i] ); 
                    n2 = sqrt( point1[0]*point1[0]+point1[1]*point1[1]+point1[2]*point1[2] );                 
                    res_j+= w1[l]*w2[k]*w3[j]*w4[i]*chi[i]*eta1[l]*(1-chi[i])*cexp((1*I)*kappa*n1)*phii[0]*phij[0]/(n2*4.0*PI);
                 }                                                                                            
            }                                                                                                 
        }                                                                                                     
    } 
    res[0] += creal(res_j);                                                                                   
    res[1] += cimag(res_j);                                                                                   
    res_j = 0.0+0.0*I; 
    for(i=0; i< sz2[0]; i++)                                                                                  
    {                                                                                                         
        for(j=0; j< sz2[1]; j++)                                                                              
        {                                                                                                     
            for(k=0; k< sz[0]; k++)                                                                           
            {                                                                                                 
                for(l=0; l< sz[1]; l++)                                                                       
                {                                                                                             
                    point1[0] = eta1[l]*eta2[k]*vert3[0]+eta1[l]*vert4[0]+(1-eta1[l])*vert5[0];   
                    point1[1] = eta1[l]*eta2[k]*vert3[1]+eta1[l]*vert4[1]+(1-eta1[l])*vert5[1];   
                    point1[2] = eta1[l]*eta2[k]*vert3[2]+eta1[l]*vert4[2]+(1-eta1[l])*vert5[2];   
                    n1 = sqrt( (point1[0]*point1[0]+point1[1]*point1[1]+point1[2]*point1[2])*chi[i]*chi[i] ); 
                    n2 = sqrt( point1[0]*point1[0]+point1[1]*point1[1]+point1[2]*point1[2] );                 
                    res_j+= w1[l]*w2[k]*w3[j]*w4[i]*chi[i]*eta1[l]*(1-chi[i])*cexp((1*I)*kappa*n1)*phii[0]*phij[0]/(n2*4.0*PI);
                 }                                                                                            
            }                                                                                                 
        }                                                                                                     
    } 
    /*gsl_matrix_set(J1,0,0, vert1[3]-vert1[0]);                                                                         
    gsl_matrix_set(J1,0,1, vert1[6]-vert1[3]);                                                                         
    gsl_matrix_set(J1,1,0, vert1[4]-vert1[1]);                                                                         
    gsl_matrix_set(J1,1,1, vert1[7]-vert1[4]);                                                                         
    gsl_matrix_set(J1,2,0, vert1[5]-vert1[2]);                                                                         
    gsl_matrix_set(J1,2,1, vert1[8]-vert1[5]);                                                                               
    gsl_blas_dgemm(CblasTrans,CblasNoTrans, 1.0, J1, J1, 0.0, G1);                                               
    det1 = sqrt(fabs(gsl_matrix_get(G1,0,0)*gsl_matrix_get(G1,1,1)- gsl_matrix_get(G1,0,1)*gsl_matrix_get(G1,1,0))); 
    gsl_matrix_set(J2,0,0, vert2[3]-vert2[0]);                                                                      
    gsl_matrix_set(J2,0,1, vert2[6]-vert2[3]);                                                                      
    gsl_matrix_set(J2,1,0, vert2[4]-vert2[1]);                                                                      
    gsl_matrix_set(J2,1,1, vert2[7]-vert2[4]);                                                                      
    gsl_matrix_set(J2,2,0, vert2[5]-vert2[2]);                                                                      
    gsl_matrix_set(J2,2,1, vert2[8]-vert2[5]);
    gsl_blas_dgemm(CblasTrans,CblasNoTrans, 1.0, J2, J2, 0.0, G2);                                            
    det2 = sqrt(fabs(gsl_matrix_get(G2,0,0)*gsl_matrix_get(G2,1,1)- gsl_matrix_get(G2,0,1)*gsl_matrix_get(G2,1,0)));                                                 
    printf("%f %f \n", det1, det2);
    for(i=0; i< sz2[0]; i++)                                                                                  
    {                                                                                                         
        for(j=0; j< sz2[1]; j++)                                                                              
        {                                                                                                     
            for(k=0; k< sz[0]; k++)                                                                           
            {                                                                                                 
                for(l=0; l< sz[1]; l++)                                                                       
                {                                                                                                                                                    
                    z[0] = -eta1[l]*eta2[k];    
                    z[1] = eta1[l]*(1-eta2[k]);
                    z[2] = eta1[l]*eta3[j];                                                                  
                    z2[0] = -eta2[k];                                                                          
                    z2[1] = (1-eta2[k]);
                    z2[2] = eta3[j];       
                    x[0] = 1;
                    x[1] = z[2];                                                                              
                    y[0] = x[0]+z[0];                                                                         
                    y[1] = z[1]; 
                    point1[0] = chi[i]*((vert1[3]-vert1[0])*x[0]+(vert1[6]-vert1[3])*x[1]-(vert2[3]-vert2[0])*y[0]-(vert2[6]-vert2[3])*y[1]);
                    point1[1] = chi[i]*((vert1[4]-vert1[1])*x[0]+(vert1[7]-vert1[4])*x[1]-(vert2[4]-vert2[1])*y[0]-(vert2[7]-vert2[4])*y[1]);
                    point1[2] = chi[i]*((vert1[5]-vert1[2])*x[0]+(vert1[8]-vert1[5])*x[1]-(vert2[5]-vert2[2])*y[0]-(vert2[8]-vert2[5])*y[1]);
                    x[0] = chi[i]*x[0];                                                                       
                    x[1] = chi[i]*x[1];                                                                       
                    y[0] = chi[i]*y[0];                                                                       
                    y[1] = chi[i]*y[1];
                    //point1[0] = (vert1[3]-vert1[0])*z[0]+(vert1[6]-vert1[3])*z[2]+(vert2[6]-vert2[3])*z[1];
                    //point1[1] = (vert1[4]-vert1[1])*z[0]+(vert1[7]-vert1[4])*z[2]+(vert2[7]-vert2[4])*z[1]; 
                    //point1[2] = (vert1[5]-vert1[2])*z[0]+(vert1[8]-vert1[5])*z[2]+(vert2[8]-vert2[5])*z[1]; 
                    //point2[0] = (vert1[3]-vert1[0])*z2[0]+(vert1[6]-vert1[3])*z2[2]+(vert2[6]-vert2[3])*z2[1];         
                    //point2[1] = (vert1[4]-vert1[1])*z2[0]+(vert1[7]-vert1[4])*z2[2]+(vert2[7]-vert2[4])*z2[1];                                       
                    //point2[2] = (vert1[5]-vert1[2])*z2[0]+(vert1[8]-vert1[5])*z2[2]+(vert2[8]-vert2[5])*z2[1]; 
                    
                    n1 = sqrt( (point1[0]*point1[0]+point1[1]*point1[1]+point1[2]*point1[2])*chi[i]*chi[i] );                                         
                    //n2 = sqrt( point2[0]*point2[0]+point2[1]*point2[1]+point2[2]*point2[2] );                                                     
                    n2 = sqrt( point1[0]*point1[0]+point1[1]*point1[1]+point1[2]*point1[2] );
                    res_j+= w1[l]*w2[k]*w3[j]*w4[i]*chi[i]*chi[i]*eta1[l]*eta1[l]*cexp((1*I)*kappa*n1)*phii[0]*phij[0]/(n2*4.0*PI);
                                                                                                              
                }                                                                                             
            }                                                                                                 
        }                                                                                                     
    }
    res[0] += creal(res_j);                                                                                   
    res[1] += cimag(res_j);                                                                                   
    res_j = 0.0+0.0*I;                                                                                        
    for(i=0; i< sz2[0]; i++)                                                                                  
    {                                                                                                         
        for(j=0; j< sz2[1]; j++)                                                                              
        {                                                                                                     
            for(k=0; k< sz[0]; k++)                                                                           
            {                                                                                                 
                for(l=0; l< sz[1]; l++)                                                                       
                {                                                                                                                             
                    z[0] = -eta1[l]*eta2[k]*eta3[j];                                                                   
                    z[1] = eta2[k]*eta1[l]*(1-eta3[j]);
                    z[2] = eta1[l];                                                           
                    z2[0] = -eta2[k]*eta3[j];                                                                                
                    z2[1] = eta2[k]*(1-eta3[j]);
                    z2[2] = 1;   
                    x[0] = 1;
                    x[1] = z[2];                                                                              
                    y[0] = x[0]+z[0];                                                                         
                    y[1] = z[1]; 
                    point1[0] = chi[i]*((vert1[3]-vert1[0])*x[0]+(vert1[6]-vert1[3])*x[1]-(vert2[3]-vert2[0])*y[0]-(vert2[6]-vert2[3])*y[1]);
                    point1[1] = chi[i]*((vert1[4]-vert1[1])*x[0]+(vert1[7]-vert1[4])*x[1]-(vert2[4]-vert2[1])*y[0]-(vert2[7]-vert2[4])*y[1]);
                    point1[2] = chi[i]*((vert1[5]-vert1[2])*x[0]+(vert1[8]-vert1[5])*x[1]-(vert2[5]-vert2[2])*y[0]-(vert2[8]-vert2[5])*y[1]);
                    x[0] = chi[i]*x[0];                                                                       
                    x[1] = chi[i]*x[1];                                                                       
                    y[0] = chi[i]*y[0];                                                                       
                    y[1] = chi[i]*y[1];
                    //point1[0] = (vert1[3]-vert1[0])*z[0]+(vert1[6]-vert1[3])*z[2]+(vert2[6]-vert2[3])*z[1];   
                    //point1[1] = (vert1[4]-vert1[1])*z[0]+(vert1[7]-vert1[4])*z[2]+(vert2[7]-vert2[4])*z[1];   
                    //point1[2] = (vert1[5]-vert1[2])*z[0]+(vert1[8]-vert1[5])*z[2]+(vert2[8]-vert2[5])*z[1];   
                    //point2[0] = (vert1[3]-vert1[0])*z2[0]+(vert1[6]-vert1[3])*z2[2]+(vert2[6]-vert2[3])*z2[1];
                    //point2[1] = (vert1[4]-vert1[1])*z2[0]+(vert1[7]-vert1[4])*z2[2]+(vert2[7]-vert2[4])*z2[1];
                    //point2[2] = (vert1[5]-vert1[2])*z2[0]+(vert1[8]-vert1[5])*z2[2]+(vert2[8]-vert2[5])*z2[1];
                    
                    n1 = sqrt( (point1[0]*point1[0]+point1[1]*point1[1]+point1[2]*point1[2])*chi[i]*chi[i] );
                    n2 = sqrt( point1[0]*point1[0]+point1[1]*point1[1]+point1[2]*point1[2] ); 
                    //n2 = sqrt( point2[0]*point2[0]+point2[1]*point2[1]+point2[2]*point2[2] );                                
                    res_j+= w1[l]*w2[k]*w3[j]*w4[i]*chi[i]*chi[i]*eta1[l]*eta1[l]*eta2[k]*cexp((1*I)*kappa*n1)*phii[0]*phij[0]/(n2*4.0*PI);
                                                                                                              
                }                                                                                             
            }                                                                                                 
        }                                                                                                     
    }                                                                                                         
    res[0] += creal(res_j);                                                                                   
    res[1] += cimag(res_j);                                                                                   
    res_j = 0.0+0.0*I;                                                                                        
    for(i=0; i< sz2[0]; i++)                                                                                  
    {                                                                                                         
        for(j=0; j< sz2[1]; j++)                                                                              
        {                                                                                                     
            for(k=0; k< sz[0]; k++)                                                                           
            {                                                                                                 
                for(l=0; l< sz[1]; l++)                                                                       
                {                                                                                                                                             
                    z[0] = eta1[l]*eta2[k];                                                          
                    z[1] = eta2[k]*eta1[l]*eta3[j];
                    z[2] = eta1[l]*(1-eta2[k]);                                           
                    z2[0] = eta2[k];                                                                         
                    z2[1] = eta2[k]*eta3[j];
                    z2[2] = (1-eta2[k]);
                    x[0] = 1-eta1[l]*eta2[k];  
                    x[1] = z[2];                                                                              
                    y[0] = x[0]+z[0];                                                                         
                    y[1] = z[1];                                                                        
                    point1[0] = chi[i]*((vert1[3]-vert1[0])*x[0]+(vert1[6]-vert1[3])*x[1]-(vert2[3]-vert2[0])*y[0]-(vert2[6]-vert2[3])*y[1]);
                    point1[1] = chi[i]*((vert1[4]-vert1[1])*x[0]+(vert1[7]-vert1[4])*x[1]-(vert2[4]-vert2[1])*y[0]-(vert2[7]-vert2[4])*y[1]);
                    point1[2] = chi[i]*((vert1[5]-vert1[2])*x[0]+(vert1[8]-vert1[5])*x[1]-(vert2[5]-vert2[2])*y[0]-(vert2[8]-vert2[5])*y[1]);
                    x[0] = chi[i]*x[0];                                                                       
                    x[1] = chi[i]*x[1];                                                                       
                    y[0] = chi[i]*y[0];                                                                       
                    y[1] = chi[i]*y[1];
                    //point1[0] = (vert1[3]-vert1[0])*z[0]+(vert1[6]-vert1[3])*z[2]+(vert2[6]-vert2[3])*z[1];   
                    //point1[1] = (vert1[4]-vert1[1])*z[0]+(vert1[7]-vert1[4])*z[2]+(vert2[7]-vert2[4])*z[1];   
                    //point1[2] = (vert1[5]-vert1[2])*z[0]+(vert1[8]-vert1[5])*z[2]+(vert2[8]-vert2[5])*z[1];   
                    //point2[0] = (vert1[3]-vert1[0])*z2[0]+(vert1[6]-vert1[3])*z2[2]+(vert2[6]-vert2[3])*z2[1];
                    //point2[1] = (vert1[4]-vert1[1])*z2[0]+(vert1[7]-vert1[4])*z2[2]+(vert2[7]-vert2[4])*z2[1];
                    //point2[2] = (vert1[5]-vert1[2])*z2[0]+(vert1[8]-vert1[5])*z2[2]+(vert2[8]-vert2[5])*z2[1];
                    
                    n1 = sqrt( (point1[0]*point1[0]+point1[1]*point1[1]+point1[2]*point1[2])*chi[i]*chi[i] ); 
                    n2 = sqrt( point1[0]*point1[0]+point1[1]*point1[1]+point1[2]*point1[2] );
                    //n2 = sqrt( point2[0]*point2[0]+point2[1]*point2[1]+point2[2]*point2[2] );
                    res_j+= w1[l]*w2[k]*w3[j]*w4[i]*chi[i]*chi[i]*eta1[l]*eta1[l]*eta2[k]*cexp((1*I)*kappa*n1)*phii[0]*phij[0]/(n2*4.0*PI);
                                                                                                              
                }                                                                                             
            }                                                                                                 
        }                                                                                                     
    } 
    res[0] += creal(res_j);                                                                                   
    res[1] += cimag(res_j);                                                         
    res_j = 0.0+0.0*I;                          
    for(i=0; i< sz2[0]; i++)                                                                                  
    {                                                                                                         
        for(j=0; j< sz2[1]; j++)                                                                              
        {                                                                                                     
            for(k=0; k< sz[0]; k++)                                                                           
            {                                                                                                 
                for(l=0; l< sz[1]; l++)                                                                       
                {                                                                                             
                    z[0] = eta1[l]*eta2[k]*eta3[j];                                                                   
                    z[1] = eta1[l];                                                           
                    z[2] = eta1[l]*eta2[k]*(1-eta3[j]);                                                               
                    z2[0] = eta2[k]*eta3[j];                                                                           
                    z2[1] = 1;                                                                  
                    z2[2] = eta2[k]*(1-eta3[j]);
                    x[0] = 1-eta1[l]*eta2[k]*eta3[j]; 
                    x[1] = z[2];                                                                              
                    y[0] = x[0]+z[0];                                                                         
                    y[1] = z[1];                                                                    
                    point1[0] = chi[i]*((vert1[3]-vert1[0])*x[0]+(vert1[6]-vert1[3])*x[1]-(vert2[3]-vert2[0])*y[0]-(vert2[6]-vert2[3])*y[1]);
                    point1[1] = chi[i]*((vert1[4]-vert1[1])*x[0]+(vert1[7]-vert1[4])*x[1]-(vert2[4]-vert2[1])*y[0]-(vert2[7]-vert2[4])*y[1]);
                    point1[2] = chi[i]*((vert1[5]-vert1[2])*x[0]+(vert1[8]-vert1[5])*x[1]-(vert2[5]-vert2[2])*y[0]-(vert2[8]-vert2[5])*y[1]);
                    x[0] = chi[i]*x[0];                                                                       
                    x[1] = chi[i]*x[1];                                                                       
                    y[0] = chi[i]*y[0];                                                                       
                    y[1] = chi[i]*y[1];
                    //point1[0] = (vert1[3]-vert1[0])*z[0]+(vert1[6]-vert1[3])*z[2]+(vert2[6]-vert2[3])*z[1];   
                    //point1[1] = (vert1[4]-vert1[1])*z[0]+(vert1[7]-vert1[4])*z[2]+(vert2[7]-vert2[4])*z[1];   
                    //point1[2] = (vert1[5]-vert1[2])*z[0]+(vert1[8]-vert1[5])*z[2]+(vert2[8]-vert2[5])*z[1];   
                    //point2[0] = (vert1[3]-vert1[0])*z2[0]+(vert1[6]-vert1[3])*z2[2]+(vert2[6]-vert2[3])*z2[1];
                    //point2[1] = (vert1[4]-vert1[1])*z2[0]+(vert1[7]-vert1[4])*z2[2]+(vert2[7]-vert2[4])*z2[1];
                    //point2[2] = (vert1[5]-vert1[2])*z2[0]+(vert1[8]-vert1[5])*z2[2]+(vert2[8]-vert2[5])*z2[1];
                    
                    n1 = sqrt( (point1[0]*point1[0]+point1[1]*point1[1]+point1[2]*point1[2])*chi[i]*chi[i] ); 
                    n2 = sqrt( point1[0]*point1[0]+point1[1]*point1[1]+point1[2]*point1[2] );
                    //n2 = sqrt( point2[0]*point2[0]+point2[1]*point2[1]+point2[2]*point2[2] );
                    res_j+= w1[l]*w2[k]*w3[j]*w4[i]*chi[i]*chi[i]*eta1[l]*eta1[l]*eta2[k]*cexp((1*I)*kappa*n1)*phii[0]*phij[0]/(n2*4.0*PI);
                                                                                                              
                }                                                                                             
            }                                                                                                 
        }                                                                                                     
    }                                                                                                         
    res[0] += creal(res_j);                                                                                   
    res[1] += cimag(res_j);                                                                                   
    res_j = 0.0+0.0*I;
    for(i=0; i< sz2[0]; i++)                                                                                  
    {                                                                                                         
        for(j=0; j< sz2[1]; j++)                                                                              
        {                                                                                                     
            for(k=0; k< sz[0]; k++)                                                                           
            {                                                                                                 
                for(l=0; l< sz[1]; l++)                                                                       
                {                                                                                             
                    z[0] = eta1[l]*eta2[k]*eta3[j];                                                           
                    z[1] = eta1[l]*eta2[k];                                                                           
                    z[2] = eta1[l]*(1-eta3[j]*eta2[k]);                                                       
                    z2[0] = eta2[k]*eta3[j];                                                                  
                    z2[1] = eta2[k];                                                                                
                    z2[2] = 1-eta3[j]*eta2[k];
                    x[0] = 1-eta1[l]*eta2[k]*eta3[j];
                    x[1] = z[2];                                                                              
                    y[0] = x[0]+z[0];                                                                         
                    y[1] = z[1]; 
                    point1[0] = chi[i]*((vert1[3]-vert1[0])*x[0]+(vert1[6]-vert1[3])*x[1]-(vert2[3]-vert2[0])*y[0]-(vert2[6]-vert2[3])*y[1]);
                    point1[1] = chi[i]*((vert1[4]-vert1[1])*x[0]+(vert1[7]-vert1[4])*x[1]-(vert2[4]-vert2[1])*y[0]-(vert2[7]-vert2[4])*y[1]);
                    point1[2] = chi[i]*((vert1[5]-vert1[2])*x[0]+(vert1[8]-vert1[5])*x[1]-(vert2[5]-vert2[2])*y[0]-(vert2[8]-vert2[5])*y[1]);
                    x[0] = chi[i]*x[0];                                                                       
                    x[1] = chi[i]*x[1];                                                                       
                    y[0] = chi[i]*y[0];                                                                       
                    y[1] = chi[i]*y[1];
                    //point1[0] = (vert1[3]-vert1[0])*z[0]+(vert1[6]-vert1[3])*z[2]+(vert2[6]-vert2[3])*z[1];   
                    //point1[1] = (vert1[4]-vert1[1])*z[0]+(vert1[7]-vert1[4])*z[2]+(vert2[7]-vert2[4])*z[1];   
                    //point1[2] = (vert1[5]-vert1[2])*z[0]+(vert1[8]-vert1[5])*z[2]+(vert2[8]-vert2[5])*z[1];   
                    //point2[0] = (vert1[3]-vert1[0])*z2[0]+(vert1[6]-vert1[3])*z2[2]+(vert2[6]-vert2[3])*z2[1];
                    //point2[1] = (vert1[4]-vert1[1])*z2[0]+(vert1[7]-vert1[4])*z2[2]+(vert2[7]-vert2[4])*z2[1];
                    //point2[2] = (vert1[5]-vert1[2])*z2[0]+(vert1[8]-vert1[5])*z2[2]+(vert2[8]-vert2[5])*z2[1];
                    
                    n1 = sqrt( (point1[0]*point1[0]+point1[1]*point1[1]+point1[2]*point1[2])*chi[i]*chi[i] ); 
                    n2 = sqrt( point1[0]*point1[0]+point1[1]*point1[1]+point1[2]*point1[2] );
                    //n2 = sqrt( point2[0]*point2[0]+point2[1]*point2[1]+point2[2]*point2[2] );
                    res_j+= w1[l]*w2[k]*w3[j]*w4[i]*chi[i]*chi[i]*eta1[l]*eta1[l]*eta2[k]*cexp((1*I)*kappa*n1)*phii[0]*phij[0]/(n2*4.0*PI);
                                                                                                              
                }                                                                                             
            }                                                                                                 
        }                                                                                                     
    }               */                                                                                          
    res[0] += creal(res_j);                                                                                   
    res[1] += cimag(res_j);
    res[0] = res[0]*4*A_t1*A_t2;
    res[1] = res[1]*4*A_t1*A_t2;
                                                                                                              
}

void singular_maxwell_weakly_cv(double kappa, int *order, double *v1, double *v2,  double *p01, double *p02, double *res, gsl_vector *n, double len1, double len2, double A_t1, double A_t2, int *dofs1, int *dofs2, int opt1, int opt2)
{                                                                                                             
    double eta1[4], eta2[4], eta3[4], chi[4], w1[4], w2[4], w3[4], w4[4];                                     
    double vert1[9] = {v1[0], v1[1], v1[2], v1[3], v1[4], v1[5], v1[6], v1[7], v1[8]};                        
    double vert2[9] = {v2[0], v2[1], v2[2], v2[3], v2[4], v2[5], v2[6], v2[7], v2[8]};                        
    int sz[2], sz2[2];                                                                                        
    double x[2], y[2], point1[3], point2[3], point3[3], point4[3];                                        
    int i,j,k,l;                                                                                              
    double n1, n2;
    double complex res_j = 0.0+0.0*I;                                                                         
    res[0] = 0.0;                                                                                             
    res[1] = 0.0;                                                                                             
    quad_rules_ab( 0.0, 1.0, eta1, eta2, w1, w2, sz, order[2]);                                               
    quad_rules_ab( 0.0, 1.0, eta3, chi, w3, w4, sz2, order[2]);                                               
    double phii[3], phij[3];                                                                          


                                                                                                              
    if( dofs1[0] == dofs2[0] )                                                                                 
    {                                                                                                     

    }
    else if( dofs1[0] == dofs2[1] )                                                                            
    {                                                                                                         
        /*vert1[0] = v1[0];                                                                                     
        vert1[1] = v1[1];                                                                                     
        vert1[2] = v1[2];                                                                                     
        vert1[3] = v1[3];                                                                                     
        vert1[4] = v1[4];                                                                                     
        vert1[5] = v1[5];                                                                                     
        vert1[6] = v1[6];                                                                                     
        vert1[7] = v1[7];                                                                                     
        vert1[8] = v1[8];*/                                                                                     
        vert2[0] = v2[3];                                                                                     
        vert2[1] = v2[4];                                                                                     
        vert2[2] = v2[5];                                                                                     
        vert2[3] = v2[6];                                                                                     
        vert2[4] = v2[7];                                                                                     
        vert2[5] = v2[8];                                                                                     
        vert2[6] = v2[0];                                                                                     
        vert2[7] = v2[1];                                                                                     
        vert2[8] = v2[2];

    }                                                                                                         
    else if( dofs1[0] == dofs2[2] )                                                                            
    {                                                                                                         
        /*vert1[0] = v1[0];                                                                                     
        vert1[1] = v1[1];                                                                                     
        vert1[2] = v1[2];                                                                                     
        vert1[3] = v1[3];                                                                                     
        vert1[4] = v1[4];                                                                                     
        vert1[5] = v1[5];                                                                                     
        vert1[6] = v1[6];                                                                                     
        vert1[7] = v1[7];                                                                                     
        vert1[8] = v1[8];*/                                                                                     
        vert2[0] = v2[6];                                                                                     
        vert2[1] = v2[7];                                                                                     
        vert2[2] = v2[8];                                                                                     
        vert2[3] = v2[0];                                                                                     
        vert2[4] = v2[1];                                                                                     
        vert2[5] = v2[2];                                                                                     
        vert2[6] = v2[3];                                                                                     
        vert2[7] = v2[4];                                                                                     
        vert2[8] = v2[5];

    }
    else if( dofs1[1] == dofs2[0] )                                                                            
    {                                                                                                         
        vert1[0] = v1[3];                                                                                     
        vert1[1] = v1[4];                                                                                     
        vert1[2] = v1[5];                                                                                     
        vert1[3] = v1[6];                                                                                     
        vert1[4] = v1[7];                                                                                     
        vert1[5] = v1[8];                                                                                     
        vert1[6] = v1[0];                                                                                     
        vert1[7] = v1[1];                                                                                     
        vert1[8] = v1[2];                                                                                     
        /*vert2[0] = v2[0];                                                                                     
        vert2[1] = v2[1];                                                                                     
        vert2[2] = v2[2];                                                                                     
        vert2[3] = v2[3];                                                                                     
        vert2[4] = v2[4];                                                                                     
        vert2[5] = v2[5];                                                                                     
        vert2[6] = v2[6];                                                                                     
        vert2[7] = v2[7];                                                                                     
        vert2[8] = v2[8];*/

    }                                                                                                         
    else if( dofs1[1] == dofs2[1] )                                                                            
    {                                                                                                         
        vert1[0] = v1[3];                                                                                     
        vert1[1] = v1[4];                                                                                     
        vert1[2] = v1[5];                                                                                     
        vert1[3] = v1[6];                                                                                     
        vert1[4] = v1[7];                                                                                     
        vert1[5] = v1[8];                                                                                     
        vert1[6] = v1[0];                                                                                     
        vert1[7] = v1[1];                                                                                     
        vert1[8] = v1[2];                                                                                     
        vert2[0] = v2[3];                                                                                     
        vert2[1] = v2[4];                                                                                     
        vert2[2] = v2[5];                                                                                     
        vert2[3] = v2[6];                                                                                     
        vert2[4] = v2[7];                                                                                     
        vert2[5] = v2[8];                                                                                     
        vert2[6] = v2[0];                                                                                     
        vert2[7] = v2[1];                                                                                     
        vert2[8] = v2[2]; 

    }
    else if( dofs1[1] == dofs2[2] )                                                                            
    {                                                                                                         
        vert1[0] = v1[3];                                                                                     
        vert1[1] = v1[4];                                                                                     
        vert1[2] = v1[5];                                                                                     
        vert1[3] = v1[6];                                                                                     
        vert1[4] = v1[7];                                                                                     
        vert1[5] = v1[8];                                                                                     
        vert1[6] = v1[0];                                                                                     
        vert1[7] = v1[1];                                                                                     
        vert1[8] = v1[2];                                                                                     
        vert2[0] = v2[6];                                                                                     
        vert2[1] = v2[7];                                                                                     
        vert2[2] = v2[8];                                                                                     
        vert2[3] = v2[0];                                                                                     
        vert2[4] = v2[1];                                                                                     
        vert2[5] = v2[2];                                                                                     
        vert2[6] = v2[3];                                                                                     
        vert2[7] = v2[4];                                                                                     
        vert2[8] = v2[5];                                                                                     

    }                                                                                                         
    else if( dofs1[2] == dofs2[0] )                                                                            
    {                                                                                                         
        vert1[0] = v1[6];                                                                                     
        vert1[1] = v1[7];                                                                                     
        vert1[2] = v1[8];                                                                                     
        vert1[3] = v1[0];                                                                                     
        vert1[4] = v1[1];                                                                                     
        vert1[5] = v1[2];                                                                                     
        vert1[6] = v1[3];                                                                                     
        vert1[7] = v1[4];                                                                                     
        vert1[8] = v1[5];                                                                                     
        /*vert2[0] = v2[0];                                                                                     
        vert2[1] = v2[1];                                                                                     
        vert2[2] = v2[2];                                                                                     
        vert2[3] = v2[3];                                                                                     
        vert2[4] = v2[4];                                                                                     
        vert2[5] = v2[5];                                                                                     
        vert2[6] = v2[6];                                                                                     
        vert2[7] = v2[7];                                                                                     
        vert2[8] = v2[8]; */                                                                                    

    }
    else if( dofs1[2] == dofs2[1] )                                                                            
    {                                                                                                         
        vert1[0] = v1[6];                                                                                     
        vert1[1] = v1[7];                                                                                     
        vert1[2] = v1[8];                                                                                     
        vert1[3] = v1[0];                                                                                     
        vert1[4] = v1[1];                                                                                     
        vert1[5] = v1[2];                                                                                     
        vert1[6] = v1[3];                                                                                     
        vert1[7] = v1[4];                                                                                     
        vert1[8] = v1[5];                                                                                     
        vert2[0] = v2[3];                                                                                     
        vert2[1] = v2[4];                                                                                     
        vert2[2] = v2[5];                                                                                     
        vert2[3] = v2[6];                                                                                     
        vert2[4] = v2[7];                                                                                     
        vert2[5] = v2[8];                                                                                     
        vert2[6] = v2[0];                                                                                     
        vert2[7] = v2[1];                                                                                     
        vert2[8] = v2[2];                                                                                     

    }                                                                                                         
    else if( dofs1[2] == dofs2[2] )                                                                            
    {                                                                                                         
        vert1[0] = v1[6];                                                                                     
        vert1[1] = v1[7];                                                                                     
        vert1[2] = v1[8];                                                                                     
        vert1[3] = v1[0];                                                                                     
        vert1[4] = v1[1];                                                                                     
        vert1[5] = v1[2];                                                                                     
        vert1[6] = v1[3];                                                                                     
        vert1[7] = v1[4];                                                                                     
        vert1[8] = v1[5];                                                                                     
        vert2[0] = v2[6];                                                                                     
        vert2[1] = v2[7];                                                                                     
        vert2[2] = v2[8];                                                                                     
        vert2[3] = v2[0];                                                                                     
        vert2[4] = v2[1];                                                                                     
        vert2[5] = v2[2];                                                                                     
        vert2[6] = v2[3];                                                                                     
        vert2[7] = v2[4];                                                                                     
        vert2[8] = v2[5];                                                                                     

    }


    for(i=0; i< sz2[0]; i++)                                                                                  
    {                                                                                                         
        for(j=0; j< sz2[1]; j++)                                                                              
        {                                                                                                     
            for(k=0; k< sz[0]; k++)                                                                           
            {                                                                                                 
                for(l=0; l< sz[1]; l++)                                                                       
                {                                                                                             
                    x[0] = 1;                                                                                 
                    x[1] = eta1[l];                                                                           
                    y[0] = eta2[k];                                                                           
                    y[1] = eta2[k]*eta3[j];
                    point1[0] = chi[i]*(vert1[0]*(x[0]-y[0])+vert1[3]*(y[0]-y[1])+y[1]*vert1[6]-(x[0]-x[1])*vert2[3]-x[1]*vert2[6]);
                    point1[1] = chi[i]*(vert1[1]*(x[0]-y[0])+vert1[4]*(y[0]-y[1])+y[1]*vert1[7]-(x[0]-x[1])*vert2[4]-x[1]*vert2[7]);
                    point1[2] = chi[i]*(vert1[2]*(x[0]-y[0])+vert1[5]*(y[0]-y[1])+y[1]*vert1[8]-(x[0]-x[1])*vert2[5]-x[1]*vert2[8]);
                    point2[0] = (vert1[0]*(x[0]-y[0])+vert1[3]*(y[0]-y[1])+y[1]*vert1[6]-(x[0]-x[1])*vert2[3]-x[1]*vert2[6]);
                    point2[1] = (vert1[1]*(x[0]-y[0])+vert1[4]*(y[0]-y[1])+y[1]*vert1[7]-(x[0]-x[1])*vert2[4]-x[1]*vert2[7]);
                    point2[2] = (vert1[2]*(x[0]-y[0])+vert1[5]*(y[0]-y[1])+y[1]*vert1[8]-(x[0]-x[1])*vert2[5]-x[1]*vert2[8]);
                    point3[0] = (1-x[0])*vert1[0]+(x[0]-x[1])*vert1[3]+x[1]*vert1[6];                         
                    point3[1] = (1-x[0])*vert1[1]+(x[0]-x[1])*vert1[4]+x[1]*vert1[7];                         
                    point3[2] = (1-x[0])*vert1[2]+(x[0]-x[1])*vert1[5]+x[1]*vert1[8];                         
                    point4[0] = (1-y[0])*vert2[0]+(y[0]-y[1])*vert2[3]+y[1]*vert2[6];                         
                    point4[1] = (1-y[0])*vert2[1]+(y[0]-y[1])*vert2[4]+y[1]*vert2[7];                         
                    point4[2] = (1-y[0])*vert2[2]+(y[0]-y[1])*vert2[5]+y[1]*vert2[8];                                                                    
                    /*point1[0] = chi[i]*((vert1[3]-vert1[0])*x[0]+(vert1[6]-vert1[3])*x[1]-(vert2[3]-vert2[0])*y[0]-(vert2[6]-vert2[3])*y[1]);
                    point1[1] = chi[i]*((vert1[4]-vert1[1])*x[0]+(vert1[7]-vert1[4])*x[1]-(vert2[4]-vert2[1])*y[0]-(vert2[7]-vert2[4])*y[1]);
                    point1[2] = chi[i]*((vert1[5]-vert1[2])*x[0]+(vert1[8]-vert1[5])*x[1]-(vert2[5]-vert2[2])*y[0]-(vert2[8]-vert2[5])*y[1]);
                    point2[0] = (vert1[3]-vert1[0])*x[0]+(vert1[6]-vert1[3])*x[1]-(vert2[3]-vert2[0])*y[0]-(vert2[6]-vert2[3])*y[1];
                    point2[1] = (vert1[4]-vert1[1])*x[0]+(vert1[7]-vert1[4])*x[1]-(vert2[4]-vert2[1])*y[0]-(vert2[7]-vert2[4])*y[1];
                    point2[2] = (vert1[5]-vert1[2])*x[0]+(vert1[8]-vert1[5])*x[1]-(vert2[5]-vert2[2])*y[0]-(vert2[8]-vert2[5])*y[1];
                    point3[0] = chi[i]*((vert1[3]-vert1[0])*x[0]+(vert1[6]-vert1[3])*x[1]) + vert1[0];
                    point3[1] = chi[i]*((vert1[4]-vert1[1])*x[0]+(vert1[7]-vert1[4])*x[1]) + vert1[1]; 
                    point3[2] = chi[i]*((vert1[5]-vert1[2])*x[0]+(vert1[8]-vert1[5])*x[1]) + vert1[2]; 
                    point4[0] = chi[i]*((vert2[3]-vert2[0])*y[0]-(vert2[6]-vert2[3])*y[1]) + vert2[0];
                    point4[1] = chi[i]*((vert2[4]-vert2[1])*y[0]-(vert2[7]-vert2[4])*y[1]) + vert2[1];
                    point4[2] = chi[i]*((vert2[5]-vert2[2])*y[0]-(vert2[8]-vert2[5])*y[1]) + vert2[2];
                    */
                    ev_basis_function_maxwell(point3, phii, p01, A_t1, len1, NULL, 1, opt1, "RWG");           
                    ev_basis_function_maxwell(point4, phij, p02, A_t2, len2, NULL, 1, opt2, "RWG");
                    n1 = sqrt(point1[0]*point1[0]+point1[1]*point1[1]+point1[2]*point1[2]);                   
                    n2 = sqrt(point2[0]*point2[0]+point2[1]*point2[1]+point2[2]*point2[2]);                   
                    res_j+= w1[l]*w2[k]*w3[j]*w4[i]*chi[i]*chi[i]*eta2[k]*cexp((1*I)*kappa*n1)*dot_prod(3,0,0,phii,phij)/(n2*4.0*PI);
                                                                                                              
                }                                                                                             
            }                                                                                                 
        }                                                                                                     
    }                                                                                                         
    res[0] += creal(res_j);                                                                                   
    res[1] += cimag(res_j);                                                                                   
    res_j = 0.0+0.0*I;
    for(i=0; i< sz2[0]; i++)                                                                                  
    {                                                                                                         
        for(j=0; j< sz2[1]; j++)                                                                              
        {                                                                                                     
            for(k=0; k< sz[0]; k++)                                                                           
            {                                                                                                 
                for(l=0; l< sz[1]; l++)                                                                       
                {                                                                                             
                    x[0] = eta2[k];                                                                           
                    x[1] = eta2[k]*eta3[j];                                                                   
                    y[0] = 1;                                                                                 
                    y[1] = eta1[l];                                                                           
                    point1[0] = chi[i]*(vert1[0]*(x[0]-y[0])+vert1[3]*(y[0]-y[1])+y[1]*vert1[6]-(x[0]-x[1])*vert2[3]-x[1]*vert2[6]);
                    point1[1] = chi[i]*(vert1[1]*(x[0]-y[0])+vert1[4]*(y[0]-y[1])+y[1]*vert1[7]-(x[0]-x[1])*vert2[4]-x[1]*vert2[7]);
                    point1[2] = chi[i]*(vert1[2]*(x[0]-y[0])+vert1[5]*(y[0]-y[1])+y[1]*vert1[8]-(x[0]-x[1])*vert2[5]-x[1]*vert2[8]);
                    point2[0] = (vert1[0]*(x[0]-y[0])+vert1[3]*(y[0]-y[1])+y[1]*vert1[6]-(x[0]-x[1])*vert2[3]-x[1]*vert2[6]);
                    point2[1] = (vert1[1]*(x[0]-y[0])+vert1[4]*(y[0]-y[1])+y[1]*vert1[7]-(x[0]-x[1])*vert2[4]-x[1]*vert2[7]);
                    point2[2] = (vert1[2]*(x[0]-y[0])+vert1[5]*(y[0]-y[1])+y[1]*vert1[8]-(x[0]-x[1])*vert2[5]-x[1]*vert2[8]);
                    point3[0] = (1-x[0])*vert1[0]+(x[0]-x[1])*vert1[3]+x[1]*vert1[6];
                    point3[1] = (1-x[0])*vert1[1]+(x[0]-x[1])*vert1[4]+x[1]*vert1[7]; 
                    point3[2] = (1-x[0])*vert1[2]+(x[0]-x[1])*vert1[5]+x[1]*vert1[8];
                    point4[0] = (1-y[0])*vert2[0]+(y[0]-y[1])*vert2[3]+y[1]*vert2[6];                         
                    point4[1] = (1-y[0])*vert2[1]+(y[0]-y[1])*vert2[4]+y[1]*vert2[7];                         
                    point4[2] = (1-y[0])*vert2[2]+(y[0]-y[1])*vert2[5]+y[1]*vert2[8];  
                    /*point1[0] = chi[i]*((vert1[3]-vert1[0])*x[0]+(vert1[6]-vert1[3])*x[1]-(vert2[3]-vert2[0])*y[0]-(vert2[6]-vert2[3])*y[1]);
                    point1[1] = chi[i]*((vert1[4]-vert1[1])*x[0]+(vert1[7]-vert1[4])*x[1]-(vert2[4]-vert2[1])*y[0]-(vert2[7]-vert2[4])*y[1]);
                    point1[2] = chi[i]*((vert1[5]-vert1[2])*x[0]+(vert1[8]-vert1[5])*x[1]-(vert2[5]-vert2[2])*y[0]-(vert2[8]-vert2[5])*y[1]);
                    point2[0] = (vert1[3]-vert1[0])*x[0]+(vert1[6]-vert1[3])*x[1]-(vert2[3]-vert2[0])*y[0]-(vert2[6]-vert2[3])*y[1];
                    point2[1] = (vert1[4]-vert1[1])*x[0]+(vert1[7]-vert1[4])*x[1]-(vert2[4]-vert2[1])*y[0]-(vert2[7]-vert2[4])*y[1];
                    point2[2] = (vert1[5]-vert1[2])*x[0]+(vert1[8]-vert1[5])*x[1]-(vert2[5]-vert2[2])*y[0]-(vert2[8]-vert2[5])*y[1];
                    point3[0] = chi[i]*((vert1[3]-vert1[0])*x[0]+(vert1[6]-vert1[3])*x[1]) + vert1[0];        
                    point3[1] = chi[i]*((vert1[4]-vert1[1])*x[0]+(vert1[7]-vert1[4])*x[1]) + vert1[1];        
                    point3[2] = chi[i]*((vert1[5]-vert1[2])*x[0]+(vert1[8]-vert1[5])*x[1]) + vert1[2];        
                    point4[0] = chi[i]*((vert2[3]-vert2[0])*y[0]-(vert2[6]-vert2[3])*y[1]) + vert2[0];        
                    point4[1] = chi[i]*((vert2[4]-vert2[1])*y[0]-(vert2[7]-vert2[4])*y[1]) + vert2[1];        
                    point4[2] = chi[i]*((vert2[5]-vert2[2])*y[0]-(vert2[8]-vert2[5])*y[1]) + vert2[2]; */       
                    ev_basis_function_maxwell(point3, phii, p01, A_t1, len1, NULL, 1, opt1, "RWG");           
                    ev_basis_function_maxwell(point4, phij, p02, A_t2, len2, NULL, 1, opt2, "RWG");
                    n1 = sqrt(point1[0]*point1[0]+point1[1]*point1[1]+point1[2]*point1[2]);                   
                    n2 = sqrt(point2[0]*point2[0]+point2[1]*point2[1]+point2[2]*point2[2]);                   
                    res_j+= w1[l]*w2[k]*w3[j]*w4[i]*chi[i]*chi[i]*eta2[k]*cexp((1*I)*kappa*n1)*dot_prod(3,0,0,phii,phij)/(n2*4.0*PI);
                                                                                                              
                }                                                                                             
            }                                                                                                 
        }                                                                                                     
    }                                                                                                         
    res[0] += creal(res_j);                                                                                   
    res[1] += cimag(res_j);                                                                                   
    res[0] = res[0]*4*A_t1*A_t2;                                                                              
    res[1] = res[1]*4*A_t1*A_t2;

}

void singular_maxwell_hypersingular_cv(double kappa, int *order, double *v1, double *v2,  double *p01, double *p02, double *res, gsl_vector *n, double len1, double len2, double A_t1, double A_t2, int *dofs1, int *dofs2, int opt1, int opt2)
{                                                                                                             
    double eta1[4], eta2[4], eta3[4], chi[4], w1[4], w2[4], w3[4], w4[4];                                 
    double vert1[9] = {v1[0], v1[1], v1[2], v1[3], v1[4], v1[5], v1[6], v1[7], v1[8]};                        
    double vert2[9] = {v2[0], v2[1], v2[2], v2[3], v2[4], v2[5], v2[6], v2[7], v2[8]};                                                            
    int sz[2], sz2[2];                                                                                        
    double  x[2], y[2], point1[3], point2[3];                                        
    int i,j,k,l;                                                                                              
    double n1, n2;                                                                                     
    double complex res_j = 0.0+0.0*I;                                                                         
    res[0] = 0.0;                                                                                             
    res[1] = 0.0;                                                                                             
    quad_rules_ab( 0.0, 1.0, eta1, eta2, w1, w2, sz, order[2]);                                               
    quad_rules_ab( 0.0, 1.0, eta3, chi, w3, w4, sz2, order[2]);                                               
    double phii[sz[0]], phij[sz[0]];                                                                          
    divergence(A_t1,len1, NULL,"RWG", opt1, phii, 1);                                                         
    divergence(A_t2,len2, NULL,"RWG", opt2, phij, 1); 

    if( dofs1[0] == dofs2[0] )
    {

    }
    else if( dofs1[0] == dofs2[1] )                                                                                 
    {                                                                                                         
        /*vert1[0] = v1[0];                                                                                     
        vert1[1] = v1[1];                                                                                     
        vert1[2] = v1[2];                                                                                     
        vert1[3] = v1[3];                                                                                     
        vert1[4] = v1[4];                                                                                     
        vert1[5] = v1[5];                                                                                     
        vert1[6] = v1[6];                                                                                     
        vert1[7] = v1[7];                                                                                     
        vert1[8] = v1[8];*/                                                                                     
        vert2[0] = v2[3];                                                                                     
        vert2[1] = v2[4];                                                                                     
        vert2[2] = v2[5];                                                                                     
        vert2[3] = v2[6];                                                                                     
        vert2[4] = v2[7];                                                                                     
        vert2[5] = v2[8];                                                                                     
        vert2[6] = v2[0];                                                                                     
        vert2[7] = v2[1];                                                                                     
        vert2[8] = v2[2];

    } 
    else if( dofs1[0] == dofs2[2] )                                                                                 
    {                                                                                                         
        /*vert1[0] = v1[0];                                                                                     
        vert1[1] = v1[1];                                                                                     
        vert1[2] = v1[2];                                                                                     
        vert1[3] = v1[3];                                                                                     
        vert1[4] = v1[4];                                                                                     
        vert1[5] = v1[5];                                                                                     
        vert1[6] = v1[6];                                                                                     
        vert1[7] = v1[7];                                                                                     
        vert1[8] = v1[8]; */                                                                                    
        vert2[0] = v2[6];                                                                                     
        vert2[1] = v2[7];                                                                                     
        vert2[2] = v2[8];                                                                                     
        vert2[3] = v2[0];                                                                                     
        vert2[4] = v2[1];                                                                                     
        vert2[5] = v2[2];                                                                                     
        vert2[6] = v2[3];                                                                                     
        vert2[7] = v2[4];                                                                                     
        vert2[8] = v2[5];

    } 
    else if( dofs1[1] == dofs2[0] )                                                                                 
    {                                                                                                         
        vert1[0] = v1[3];                                                                                     
        vert1[1] = v1[4];                                                                                     
        vert1[2] = v1[5];                                                                                     
        vert1[3] = v1[6];                                                                                     
        vert1[4] = v1[7];                                                                                     
        vert1[5] = v1[8];                                                                                     
        vert1[6] = v1[0];                                                                                     
        vert1[7] = v1[1];                                                                                     
        vert1[8] = v1[2];                                                                                     
        /*vert2[0] = v2[0];                                                                                     
        vert2[1] = v2[1];                                                                                     
        vert2[2] = v2[2];                                                                                     
        vert2[3] = v2[3];                                                                                     
        vert2[4] = v2[4];                                                                                     
        vert2[5] = v2[5];                                                                                     
        vert2[6] = v2[6];                                                                                     
        vert2[7] = v2[7];                                                                                     
        vert2[8] = v2[8]; 
        */
    }                                                                                                         
    else if( dofs1[1] == dofs2[1] )                                                                            
    {                                                                                                         
        vert1[0] = v1[3];                                                                                     
        vert1[1] = v1[4];                                                                                     
        vert1[2] = v1[5];                                                                                     
        vert1[3] = v1[6];                                                                                     
        vert1[4] = v1[7];                                                                                     
        vert1[5] = v1[8];                                                                                     
        vert1[6] = v1[0];                                                                                     
        vert1[7] = v1[1];                                                                                     
        vert1[8] = v1[2];                                                                                     
        vert2[0] = v2[3];                                                                                     
        vert2[1] = v2[4];                                                                                     
        vert2[2] = v2[5];                                                                                     
        vert2[3] = v2[6];                                                                                     
        vert2[4] = v2[7];                                                                                     
        vert2[5] = v2[8];                                                                                     
        vert2[6] = v2[0];                                                                                     
        vert2[7] = v2[1];                                                                                     
        vert2[8] = v2[2];

    }
    else if( dofs1[1] == dofs2[2] )                                                                            
    {                                                                                                         
        vert1[0] = v1[3];                                                                                     
        vert1[1] = v1[4];                                                                                     
        vert1[2] = v1[5];                                                                                     
        vert1[3] = v1[6];                                                                                     
        vert1[4] = v1[7];                                                                                     
        vert1[5] = v1[8];                                                                                     
        vert1[6] = v1[0];                                                                                     
        vert1[7] = v1[1];                                                                                     
        vert1[8] = v1[2];                                                                                     
        vert2[0] = v2[6];                                                                                     
        vert2[1] = v2[7];                                                                                     
        vert2[2] = v2[8];                                                                                     
        vert2[3] = v2[0];                                                                                     
        vert2[4] = v2[1];                                                                                     
        vert2[5] = v2[2];                                                                                     
        vert2[6] = v2[3];                                                                                     
        vert2[7] = v2[4];                                                                                     
        vert2[8] = v2[5]; 

    } 
    else if( dofs1[2] == dofs2[0] )                                                                            
    {                                                                                                         
        vert1[0] = v1[6];                                                                                     
        vert1[1] = v1[7];                                                                                     
        vert1[2] = v1[8];                                                                                     
        vert1[3] = v1[0];                                                                                     
        vert1[4] = v1[1];                                                                                     
        vert1[5] = v1[2];                                                                                     
        vert1[6] = v1[3];                                                                                     
        vert1[7] = v1[4];                                                                                     
        vert1[8] = v1[5];                                                                                     
        /*vert2[0] = v2[0];                                                                                     
        vert2[1] = v2[1];                                                                                     
        vert2[2] = v2[2];                                                                                     
        vert2[3] = v2[3];                                                                                     
        vert2[4] = v2[4];                                                                                     
        vert2[5] = v2[5];                                                                                     
        vert2[6] = v2[6];                                                                                     
        vert2[7] = v2[7];                                                                                     
        vert2[8] = v2[8];
        */
    }
    else if( dofs1[2] == dofs2[1] )                                                                            
    {                                                                                                         
        vert1[0] = v1[6];                                                                                     
        vert1[1] = v1[7];                                                                                     
        vert1[2] = v1[8];                                                                                     
        vert1[3] = v1[0];                                                                                     
        vert1[4] = v1[1];                                                                                     
        vert1[5] = v1[2];                                                                                     
        vert1[6] = v1[3];                                                                                     
        vert1[7] = v1[4];                                                                                     
        vert1[8] = v1[5];                                                                                     
        vert2[0] = v2[3];                                                                                     
        vert2[1] = v2[4];                                                                                     
        vert2[2] = v2[5];                                                                                     
        vert2[3] = v2[6];                                                                                     
        vert2[4] = v2[7];                                                                                     
        vert2[5] = v2[8];                                                                                     
        vert2[6] = v2[0];                                                                                     
        vert2[7] = v2[1];                                                                                     
        vert2[8] = v2[2]; 

    }
    else if( dofs1[2] == dofs2[2] )                                                                            
    {                                                                                                         
        vert1[0] = v1[6];                                                                                     
        vert1[1] = v1[7];                                                                                     
        vert1[2] = v1[8];                                                                                     
        vert1[3] = v1[0];                                                                                     
        vert1[4] = v1[1];                                                                                     
        vert1[5] = v1[2];                                                                                     
        vert1[6] = v1[3];                                                                                     
        vert1[7] = v1[4];                                                                                     
        vert1[8] = v1[5];                                                                                     
        vert2[0] = v2[6];                                                                                     
        vert2[1] = v2[7];                                                                                     
        vert2[2] = v2[8];                                                                                     
        vert2[3] = v2[0];                                                                                     
        vert2[4] = v2[1];                                                                                     
        vert2[5] = v2[2];                                                                                     
        vert2[6] = v2[3];                                                                                     
        vert2[7] = v2[4];                                                                                     
        vert2[8] = v2[5]; 

    } 

    for(i=0; i< sz2[0]; i++)                                                                                  
    {                                                                                                         
        for(j=0; j< sz2[1]; j++)                                                                              
        {                                                                                                     
            for(k=0; k< sz[0]; k++)                                                                           
            {                                                                                                 
                for(l=0; l< sz[1]; l++)                                                                       
                {                                                                                             
                    x[0] = 1;
                    x[1] = eta1[l];                                                                  
                    y[0] = eta2[k];
                    y[1] = eta2[k]*eta3[j];                                                                                                        
                    point1[0] = chi[i]*(vert1[0]*(x[0]-y[0])+vert1[3]*(y[0]-y[1])+y[1]*vert1[6]-(x[0]-x[1])*vert2[3]-x[1]*vert2[6]);
                    point1[1] = chi[i]*(vert1[1]*(x[0]-y[0])+vert1[4]*(y[0]-y[1])+y[1]*vert1[7]-(x[0]-x[1])*vert2[4]-x[1]*vert2[7]);
                    point1[2] = chi[i]*(vert1[2]*(x[0]-y[0])+vert1[5]*(y[0]-y[1])+y[1]*vert1[8]-(x[0]-x[1])*vert2[5]-x[1]*vert2[8]);
                    point2[0] = (vert1[0]*(x[0]-y[0])+vert1[3]*(y[0]-y[1])+y[1]*vert1[6]-(x[0]-x[1])*vert2[3]-x[1]*vert2[6]);
                    point2[1] = (vert1[1]*(x[0]-y[0])+vert1[4]*(y[0]-y[1])+y[1]*vert1[7]-(x[0]-x[1])*vert2[4]-x[1]*vert2[7]);
                    point2[2] = (vert1[2]*(x[0]-y[0])+vert1[5]*(y[0]-y[1])+y[1]*vert1[8]-(x[0]-x[1])*vert2[5]-x[1]*vert2[8]);
                    /*point1[0] = chi[i]*((vert1[3]-vert1[0])*x[0]+(vert1[6]-vert1[3])*x[1]-(vert2[3]-vert2[0])*y[0]-(vert2[6]-vert2[3])*y[1]);
                    point1[1] = chi[i]*((vert1[4]-vert1[1])*x[0]+(vert1[7]-vert1[4])*x[1]-(vert2[4]-vert2[1])*y[0]-(vert2[7]-vert2[4])*y[1]);
                    point1[2] = chi[i]*((vert1[5]-vert1[2])*x[0]+(vert1[8]-vert1[5])*x[1]-(vert2[5]-vert2[2])*y[0]-(vert2[8]-vert2[5])*y[1]);
                    point2[0] = (vert1[3]-vert1[0])*x[0]+(vert1[6]-vert1[3])*x[1]-(vert2[3]-vert2[0])*y[0]-(vert2[6]-vert2[3])*y[1];
                    point2[1] = (vert1[4]-vert1[1])*x[0]+(vert1[7]-vert1[4])*x[1]-(vert2[4]-vert2[1])*y[0]-(vert2[7]-vert2[4])*y[1];
                    point2[2] = (vert1[5]-vert1[2])*x[0]+(vert1[8]-vert1[5])*x[1]-(vert2[5]-vert2[2])*y[0]-(vert2[8]-vert2[5])*y[1];
                    */
                    n1 = sqrt(point1[0]*point1[0]+point1[1]*point1[1]+point1[2]*point1[2]); 
                    n2 = sqrt(point2[0]*point2[0]+point2[1]*point2[1]+point2[2]*point2[2]);                                                      
                    res_j+= w1[l]*w2[k]*w3[j]*w4[i]*chi[i]*chi[i]*eta2[k]*cexp((1*I)*kappa*n1)*phii[0]*phij[0]/(n2*4.0*PI);
                                                                                                              
                }                                                                                             
            }                                                                                                 
        }                                                                                                     
    }                                                                                                         
    res[0] += creal(res_j);                                                                                   
    res[1] += cimag(res_j);                                                                                   
    res_j = 0.0+0.0*I;  
    for(i=0; i< sz2[0]; i++)                                                                                  
    {                                                                                                         
        for(j=0; j< sz2[1]; j++)                                                                              
        {                                                                                                     
            for(k=0; k< sz[0]; k++)                                                                           
            {                                                                                                 
                for(l=0; l< sz[1]; l++)                                                                       
                {                                                                                             
                    x[0] = eta2[k];                                                                                  
                    x[1] = eta2[k]*eta3[j];                                                                           
                    y[0] = 1;                                                                           
                    y[1] = eta1[l]; 
                    point1[0] = chi[i]*(vert1[0]*(x[0]-y[0])+vert1[3]*(y[0]-y[1])+y[1]*vert1[6]-(x[0]-x[1])*vert2[3]-x[1]*vert2[6]);
                    point1[1] = chi[i]*(vert1[1]*(x[0]-y[0])+vert1[4]*(y[0]-y[1])+y[1]*vert1[7]-(x[0]-x[1])*vert2[4]-x[1]*vert2[7]);
                    point1[2] = chi[i]*(vert1[2]*(x[0]-y[0])+vert1[5]*(y[0]-y[1])+y[1]*vert1[8]-(x[0]-x[1])*vert2[5]-x[1]*vert2[8]);
                    point2[0] = (vert1[0]*(x[0]-y[0])+vert1[3]*(y[0]-y[1])+y[1]*vert1[6]-(x[0]-x[1])*vert2[3]-x[1]*vert2[6]);
                    point2[1] = (vert1[1]*(x[0]-y[0])+vert1[4]*(y[0]-y[1])+y[1]*vert1[7]-(x[0]-x[1])*vert2[4]-x[1]*vert2[7]);
                    point2[2] = (vert1[2]*(x[0]-y[0])+vert1[5]*(y[0]-y[1])+y[1]*vert1[8]-(x[0]-x[1])*vert2[5]-x[1]*vert2[8]);
                                                                  
                    /*point1[0] = chi[i]*((vert1[3]-vert1[0])*x[0]+(vert1[6]-vert1[3])*x[1]-(vert2[3]-vert2[0])*y[0]-(vert2[6]-vert2[3])*y[1]);
                    point1[1] = chi[i]*((vert1[4]-vert1[1])*x[0]+(vert1[7]-vert1[4])*x[1]-(vert2[4]-vert2[1])*y[0]-(vert2[7]-vert2[4])*y[1]);
                    point1[2] = chi[i]*((vert1[5]-vert1[2])*x[0]+(vert1[8]-vert1[5])*x[1]-(vert2[5]-vert2[2])*y[0]-(vert2[8]-vert2[5])*y[1]);
                    point2[0] = (vert1[3]-vert1[0])*x[0]+(vert1[6]-vert1[3])*x[1]-(vert2[3]-vert2[0])*y[0]-(vert2[6]-vert2[3])*y[1];
                    point2[1] = (vert1[4]-vert1[1])*x[0]+(vert1[7]-vert1[4])*x[1]-(vert2[4]-vert2[1])*y[0]-(vert2[7]-vert2[4])*y[1];
                    point2[2] = (vert1[5]-vert1[2])*x[0]+(vert1[8]-vert1[5])*x[1]-(vert2[5]-vert2[2])*y[0]-(vert2[8]-vert2[5])*y[1];
                    */
                    n1 = sqrt(point1[0]*point1[0]+point1[1]*point1[1]+point1[2]*point1[2]);                   
                    n2 = sqrt(point2[0]*point2[0]+point2[1]*point2[1]+point2[2]*point2[2]);                   
                    res_j+= w1[l]*w2[k]*w3[j]*w4[i]*chi[i]*chi[i]*eta2[k]*cexp((1*I)*kappa*n1)*phii[0]*phij[0]/(n2*4.0*PI);
                                                                                                              
                }                                                                                             
            }                                                                                                 
        }                                                                                                     
    }                                                                                                         
    res[0] += creal(res_j);                                                                                   
    res[1] += cimag(res_j); 

    res[0] = res[0]*4*A_t1*A_t2;                                                                              
    res[1] = res[1]*4*A_t1*A_t2;


} 
void singular_maxwell_hypersingular(double *P0, double *V, double *weights, double A_t, double len, double *res, double *IE, int *sz,  int order, int opt, double kappa, char *tp, double phi)
{                                                                                                             
    double y[3], M[6], A[9], a;                                                                               
    double cres[2];                                                                                           
    int k;                                                                                                    
                                                                                                              
    cres[0] = 0.0;                                                                                            
    cres[1] = 0.0;                                                                                            
    res[0] = 0.0;                                                                                             
    res[1] = 0.0;                                                                                             
    for( k = 0; k < sz[0]; k++)                                                                               
    {                                                                                                         
            y[0] = P0[3*k];                                                                                   
            y[1] = P0[3*k+1];                                                                                 
            y[2] = P0[3*k+2];                                                                                 
            A[0] = P0[3*k];                                                                                   
            A[1] = P0[3*k+1];                                                                                 
            A[2] = P0[3*k+2];                                                                                 
            M[0] = V[0];                                                                                      
            M[1] = V[1];                                                                                      
            M[2] = V[2];                                                                                      
            M[3] = V[3];                                                                                      
            M[4] = V[4];                                                                                      
            M[5] = V[5];                                                                                      
            A[3] = V[0];                                                                                      
            A[4] = V[1];                                                                                      
            A[5] = V[2];                                                                                      
            A[6] = V[3];                                                                                      
            A[7] = V[4];                                                                                      
            A[8] = V[5];                                                                                      
            a = area(A);                                                                                
            evaluate_sing_kern_maxwell_hypersingular(y, M, V, A_t, len, order, opt, cres, kappa, tp);                                       
            res[0]+= a*cres[0];                                                                               
            res[1]+= a*cres[1];                                                                               
            M[0] = V[3];                                                                                      
            M[1] = V[4];                                                                                      
            M[2] = V[5];                                                                                      
            M[3] = V[6];                                                                                      
            M[4] = V[7];                                                                                      
            M[5] = V[8];                                                                                      
            A[3] = V[3];                                                                                      
            A[4] = V[4];                                                                                      
            A[5] = V[5];                                                                                      
            A[6] = V[6];                                                                                      
            A[7] = V[7];                                                                                      
            A[8] = V[8];                                                                                      
            a = area(A);
            evaluate_sing_kern_maxwell_hypersingular(y, M, V, A_t, len, order, opt, cres, kappa, tp);                                         
            res[0]+= a*cres[0];                                                                               
            res[1]+= a*cres[1];
            M[0] = V[6];                                                                                      
            M[1] = V[7];                                                                                      
            M[2] = V[8];                                                                                      
            M[3] = V[0];                                                                                      
            M[4] = V[1];                                                                                      
            M[5] = V[2];                                                                                      
            A[3] = V[6];                                                                                      
            A[4] = V[7];                                                                                      
            A[5] = V[8];                                                                                      
            A[6] = V[0];                                                                                      
            A[7] = V[1];                                                                                      
            A[8] = V[2];                                                                                      
            a = area(A);  
            evaluate_sing_kern_maxwell_hypersingular(y, M, V, A_t, len, order, opt, cres, kappa, tp);                                         
            res[0]+= cres[0]*a;                                                                               
            res[1]+= cres[1]*a;                                                                               
            res[0] = res[0]*phi;                                                                              
            res[1] = res[1]*phi;                                                                              
                                                                                                                                                                                           
    }                                                                                                         
} 

void integration_elements_PRIMAL( int sz_eta, int sz_chi, double *eta, double *chi, double *corners, double *IE, double *P0 )
{
    get_ie_t( sz_eta, corners, IE );
    points_t( corners, sz_eta, eta, chi, P0);
    
}


void integration_elements_DUAL( int sz_eta, int sz_chi, double *eta, double *chi, double *corners, double *IE, double *P0 )
{
    get_ie_q( sz_eta, sz_chi, corners, eta, chi, IE );
    points_q( corners, sz_eta, sz_chi, eta, chi, P0);
}

void integration_elements_DUAL_PRIMAL( int sz_eta, int sz_chi, double *corners_t, double *corners_q, double *eta, double *chi,double *eta_p, double *chi_p, double *IE, double *P0 )
{

    double points_aux[3*sz_eta*sz_chi],T[9];
    int i;
    
    get_ie_q( sz_eta, sz_chi, corners_q, eta, chi, IE );                                                  
    points_q( corners_q, sz_eta, sz_chi, eta, chi, P0);                                                   
    transpose( 3, 3, corners_t, T );
    
    if(T[0] ==0.0 && T[1] == 0.0 && T[2] ==0.0)
    {
        double T2[4], points_aux2[2*sz_eta*sz_chi];
        
        T2[0] = T[4]-T[3];
        T2[1] = T[5]-T[3];
        T2[2] = T[7]-T[6];
        T2[3] = T[8]-T[6];
        inverse(2, T2); 
        transpose( sz_eta*sz_chi, 3, P0, points_aux ); 
        for(i = 0; i<sz_eta*sz_chi; i++)
        {
            points_aux2[2*i] = points_aux[3*i+1]-T[3];
            points_aux2[2*i+1] = points_aux[3*i+2]-T[6];
        }
        matrix_prod( 2, 2, 2, sz_eta*sz_chi, T2, points_aux2 );
        for( i = 0; i < sz_eta*sz_chi; i++ )                                                                      
        {                                                                                                         
            chi_p[i] = points_aux2[2*i];                                                                        
            eta_p[i] = points_aux2[2*i+1];                                                                        
        }
    }
    else if(T[3] ==0.0 && T[4] == 0.0 && T[5] ==0.0)                                                               
    {                                                                                                         
        double T2[4], points_aux2[2*sz_eta*sz_chi];                                                           
                                                                                                              
        T2[0] = T[1]-T[0];                                                                                    
        T2[1] = T[2]-T[0];                                                                                    
        T2[2] = T[7]-T[6];                                                                                    
        T2[3] = T[8]-T[6];                                                                                    
        inverse(2, T2);                                                                                       
        transpose( sz_eta*sz_chi, 3, P0, points_aux );                                                        
        for(i = 0; i<sz_eta*sz_chi; i++)                                                                      
        {                                                                                                     
            points_aux2[2*i] = points_aux[3*i]-T[0];                                                        
            points_aux2[2*i+1] = points_aux[3*i+2]-T[6];                                                      
        }                                                                                                     
        matrix_prod( 2, 2, 2, sz_eta*sz_chi, T2, points_aux2 );                                               
        for( i = 0; i < sz_eta*sz_chi; i++ )                                                                  
        {                                                                                                     
            chi_p[i] = points_aux2[2*i];                                                                      
            eta_p[i] = points_aux2[2*i+1];                                                                    
        }                                                                                                     
    }
    else if(T[6] ==0.0 && T[7] == 0.0 && T[8] ==0.0)                                                          
    {                                                                                                         
        double T2[4], points_aux2[2*sz_eta*sz_chi];                                                           
                                                                                                              
        T2[0] = T[1]-T[0];                                                                                    
        T2[1] = T[2]-T[0];                                                                                    
        T2[2] = T[4]-T[3];                                                                                    
        T2[3] = T[5]-T[3];                                                                                    
        inverse(2, T2);                                                                                       
        transpose( sz_eta*sz_chi, 3, P0, points_aux );                                                        
        for(i = 0; i<sz_eta*sz_chi; i++)                                                                      
        {                                                                                                     
            points_aux2[2*i] = points_aux[3*i]-T[0];                                                          
            points_aux2[2*i+1] = points_aux[3*i+1]-T[3];                                                      
        }                                                                                                     
        matrix_prod( 2, 2, 2, sz_eta*sz_chi, T2, points_aux2 );                                               
        for( i = 0; i < sz_eta*sz_chi; i++ )                                                                  
        {                                                                                                     
            chi_p[i] = points_aux2[2*i];                                                                      
            eta_p[i] = points_aux2[2*i+1];                                                                    
        }                                                                                                     
    }
    else
    {
        double  points_aux2[3*sz_eta*sz_chi];
        inverse(3, T);
        transpose( sz_eta*sz_chi, 3, P0, points_aux );
        matrix_prod( 3, 3, 3, sz_eta*sz_chi, T, points_aux );
        transpose( 3, sz_eta*sz_chi, points_aux, points_aux2 ); 
        for( i = 0; i < sz_eta*sz_chi; i++ )
        {
            chi_p[i] = points_aux2[3*i+1];
            eta_p[i] = points_aux2[3*i+2];
        }
        for( i = 0; i < sz_eta*sz_chi; i++ )                                                                      
        {                                                                                                         
            chi_p[i] = points_aux2[3*i+1];                                                                        
            eta_p[i] = points_aux2[3*i+2];                                                                        
        }
    }

    
}

void change_of_basis(double *corners_t, double *eta_p, double *chi_p, int id, double *P0 )
{
    double points_aux[3],T[9];
    transpose( 3, 3, corners_t, T );
    if(T[0] ==0.0 && T[1] == 0.0 && T[2] ==0.0)
    {
        double T2[4], points_aux2[2];
        
        T2[0] = T[4]-T[3];
        T2[1] = T[5]-T[3];
        T2[2] = T[7]-T[6];
        T2[3] = T[8]-T[6];
        inverse(2, T2);
        transpose( 1, 3, P0, points_aux );
        points_aux2[0] = points_aux[1]-T[3];
        points_aux2[1] = points_aux[2]-T[6];
        matrix_prod( 2, 2, 2, 1, T2, points_aux2 );
        chi_p[id] = points_aux2[0];
        eta_p[id] = points_aux2[1];
    }
    else if(T[3] ==0.0 && T[4] == 0.0 && T[5] ==0.0)
    {
        double T2[4], points_aux2[2];
        
        T2[0] = T[1]-T[0];
        T2[1] = T[2]-T[0];
        T2[2] = T[7]-T[6];
        T2[3] = T[8]-T[6];
        inverse(2, T2);
        transpose( 1, 3, P0, points_aux );
        points_aux2[0] = points_aux[0]-T[0];
        points_aux2[1] = points_aux[2]-T[6];
        matrix_prod( 2, 2, 2, 1, T2, points_aux2 );
        chi_p[id] = points_aux2[0];
        eta_p[id] = points_aux2[1];
    }
    else if(T[6] ==0.0 && T[7] == 0.0 && T[8] ==0.0)
    {
        double T2[4], points_aux2[2];
        
        T2[0] = T[1]-T[0];
        T2[1] = T[2]-T[0];
        T2[2] = T[4]-T[3];
        T2[3] = T[5]-T[3];
        inverse(2, T2);
        transpose( 1, 3, P0, points_aux );
        points_aux2[0] = points_aux[0]-T[0];
        points_aux2[1] = points_aux[1]-T[3];
        matrix_prod( 2, 2, 2, 1, T2, points_aux2 );
        chi_p[id] = points_aux2[0];
        eta_p[id] = points_aux2[1];
    }
    else
    {
        double  points_aux2[3];
        inverse(3, T);
        transpose( 1, 3, P0, points_aux );
        matrix_prod( 3, 3, 3, 1, T, points_aux );
        transpose( 3, 1, points_aux, points_aux2 );
        chi_p[id] = points_aux2[1];
        eta_p[id] = points_aux2[2];
    }

}

