#include "ext.h"
pMeshes new_space( double *PyVp, double *PyVb, int *PyDp, int *PyDb, int *PyDq, double *PyTp, int *PyEdgep, int *PyEdgeb,  int size_of_Vp, int size_of_Vb, int size_of_Dp, int size_of_Db,  int size_of_Dq, int size_of_Tp, int size_of_Edgep, int size_of_Edgeb, char *trialbf, char *testbf )
{
    
    int                    i, j,k;
    int                options[2];
    pMeshes                  Mesh;
    int                 ndofb = 0;
    
   
    double **Vp = (double **)malloc( sizeof(double*) * size_of_Vp );
    double **Vb = (double **)malloc( sizeof(double*) * size_of_Vb );
    int **Dp = (int **)malloc( sizeof(int*) * size_of_Dp );
    int **Db = (int **)malloc( sizeof(int*) * size_of_Db );
    int **Dq = (int **)malloc( sizeof(int*) * size_of_Dq );
    double **Tp = (double **)malloc( sizeof(double*) * size_of_Tp );
    int **Ep = (int **)malloc( sizeof(int*) * size_of_Edgep );
    int **Eb = (int **)malloc( sizeof(int*) * size_of_Edgeb );


    for( i=0; i < size_of_Vp; i++ )
    {
        Vp[i] = (double *)malloc( sizeof(double) * 3 );
    }
    
    for( i=0; i < size_of_Vb; i++ )
    {
        Vb[i] = (double *)malloc( sizeof(double) * 3 );
    }
    
    for( i=0; i < size_of_Dp; i++ )
    {
        Dp[i] = (int *)malloc( sizeof(int) * 3 );
        Tp[i] = (double *)malloc( sizeof(double) * 9 );
    }
    
    for( i=0; i < size_of_Db; i++ )
    {
        Db[i] = (int *)malloc( sizeof(int) * 3 );
    }
    
    for( i=0; i < size_of_Dq; i++ )                                                                           
    {                                                                                                         
        Dq[i] = (int *)malloc( sizeof(int) * 4 );                                                             
    }
    
    for( i=0; i < size_of_Edgep; i++ )                                                                           
    {                                                                                                         
        Ep[i] = (int *)malloc( sizeof(int) * 2 );                                                             
    } 
    for( i=0; i < size_of_Edgeb; i++ )                                                                        
    {                                                                                                         
        Eb[i] = (int *)malloc( sizeof(int) * 2 );                                                             
    }

     for( i=0; i < size_of_Vp; i++ )
    {
        for( j=0; j < 3; j++ )
        {
            Vp[i][j] = PyVp[i*3+j];
        }
    }
    
    for( i=0; i < size_of_Vb; i++ )
    {
        for( j=0; j < 3; j++ )
        {
            Vb[i][j] = PyVb[i*3+j];
        }
    }
    
    for( i=0; i < size_of_Dp; i++ )
    {
        for( j=0; j < 3; j++ )
        {
            Dp[i][j] = PyDp[i*3+j];
            for( k = 0; k < 3; k++ )
            {
                Tp[i][3*j+k] = PyTp[i*9+3*j+k];
            }
        }
    }
    
    for( i=0; i < size_of_Db; i++ )
    {
        for( j=0; j < 3; j++ )
        {
            Db[i][j] = PyDb[i*3+j];
            if( Db[i][j]> ndofb )
            {
                ndofb = Db[i][j];
            }
        }
    }
    
    for( i=0; i < size_of_Dq; i++ )                                                                           
    {                                                                                                         
        for( j=0; j < 4; j++ )                                                                                
        {
            Dq[i][j] = PyDq[i*4+j];
        }
    }

    for( i=0; i < size_of_Edgep; i++ )                                                                           
    {                                                                                                         
        for( j=0; j < 2; j++ )                                                                                
        {                                                                                                     
            Ep[i][j] = PyEdgep[i*2+j];                                                                           
        }                                                                                                     
    } 
    for( i=0; i < size_of_Edgeb; i++ )                                                                        
    {                                                                                                         
        for( j=0; j < 2; j++ )                                                                                
        {                                                                                                     
            Eb[i][j] = PyEdgeb[i*2+j];                                                                        
        }                                                                                                     
    }

    if( strcmp(testbf, "P0" ) == 0)
    {
        options[1] = 0;
    }
    if (strcmp(trialbf, "P0") == 0 )
    {
        options[0] = 0;
    }
    if( strcmp(testbf, "P1" ) == 0 )
    {
        options[1] = 1;
    }
    if( strcmp(trialbf, "P1") == 0 )
    {
        options[0] = 1;
    }
    if( strcmp(testbf, "P0d" ) == 0)
    {
        options[1] = 2;
    }
    if( strcmp(trialbf, "P0d") == 0 )
    {
        options[0] = 2;
    }
    if( strcmp(testbf, "P1d" ) == 0 )
    {
        options[1] = 3;
    }
    if( strcmp(trialbf, "P1d") == 0 )
    {
        options[0] = 3;
    }
    if( strcmp(testbf, "P0b" ) == 0 )                                                                          
    {                                                                                                         
        options[1] = 4;                                                                                       
    }                                                                                                         
    if( strcmp(trialbf, "P0b") == 0 )                                                                          
    {                                                                                                         
        options[0] = 4;                                                                                       
    }
    if( strcmp(testbf, "P1b" ) == 0 )                                                                         
    {                                                                                                         
        options[1] = 5;                                                                                       
    }                                                                                                         
    if( strcmp(trialbf, "P1b") == 0 )                                                                         
    {                                                                                                         
        options[0] = 5;                                                                                       
    }
    if( strcmp(testbf, "RWG" ) == 0 )                                                                         
    {                                                                                                         
        options[1] = 6;                                                                                       
    }                                                                                                         
    if( strcmp(trialbf, "RWG") == 0 )                                                                         
    {                                                                                                         
        options[0] = 6;                                                                                       
    }
    if( strcmp(testbf, "BC" ) == 0 )                                                                         
    {                                                                                                         
        options[1] = 7;                                                                                       
    }                                                                                                         
    if( strcmp(trialbf, "BC") == 0 )                                                                         
    {                                                                                                         
        options[0] = 7;                                                                                       
    }
    Mesh = new_Mesh( Vp, Vb, Dp, Db, Dq, Tp, Ep, Eb, size_of_Dp, size_of_Db, size_of_Vp, size_of_Vb, size_of_Edgep, size_of_Edgeb, options);
    
    Mesh->trialbf = (char*)malloc((strlen(trialbf)+1) * sizeof(char));
    strcpy(Mesh->trialbf, trialbf); 
    Mesh->testbf = (char*)malloc((strlen(testbf)+1) * sizeof(char));                                        
    strcpy(Mesh->testbf, testbf);
    if( strcmp(trialbf, "P0b") == 0 )                                                                         
    {                                                                                                         
        strcpy(Mesh->trialbf, "P0");                                                                                       
    }
    if( strcmp(trialbf, "P0b") == 0 )                                                                         
    {                                                                                                         
        strcpy(Mesh->testbf, "P0");                                                                            
    }
    Mesh->Vp = Vp;
    Mesh->Vb = Vb;
    Mesh->ntp = size_of_Tp;
    Mesh->ndp = size_of_Dp;
    Mesh->ndb = size_of_Db;
    Mesh->nvp = size_of_Vp;
    Mesh->nvb = size_of_Vb;
    Mesh->nep = size_of_Edgep;
    Mesh->neb = size_of_Edgeb; 
    Mesh->ndofb = ndofb;
   
    for( i = 0; i < size_of_Dp; i++ )
    {
        free( Dp[i] );
        free( Tp[i] );
    }
    
    for( i = 0; i < size_of_Db; i++ )
    {
        free( Db[i] );
    }

    for( i = 0; i < size_of_Edgep; i++ )                                                                         
    {                                                                                                         
        free( Ep[i] );                                                                                        
    }
    for( i = 0; i < size_of_Edgeb; i++ )                                                                         
    {                                                                                                         
        free( Eb[i] );                                                                                        
    } 
    return Mesh;
}

pMesh new_mesh( double *PyVp, double *PyVb, int *PyDp, int *PyDb, int *PyDq, double *PyTp, double *PyTb, int *PyEdgep, int *PyEdgeb, int *PyEdgep2, int *PyEdgeb2, int size_of_Vp, int size_of_Vb, int size_of_Dp, int size_of_Db,  int size_of_Dq, int size_of_Tp, int size_of_Tb, int size_of_Edgep, int size_of_Edgeb)
{                                                                                                             
                                                                                                              
    int        i, j,k;                                                                                                                                                      
    pMesh           M;                                                                            
    int     ndofb = 0;                                                                            
                                                                                                              
    M = (pMesh) malloc(sizeof(Mesh));                                                                                                                
    double **Vp = (double **)malloc( sizeof(double*) * size_of_Vp );                                          
    double **Vb = (double **)malloc( sizeof(double*) * size_of_Vb );                                          
    int **Dp = (int **)malloc( sizeof(int*) * size_of_Dp );                                                   
    int **Db = (int **)malloc( sizeof(int*) * size_of_Db );                                                   
    int **Dq = (int **)malloc( sizeof(int*) * size_of_Dq );                                                   
    double **Tp = (double **)malloc( sizeof(double*) * size_of_Tp );                                          
    double **Tb = (double **)malloc( sizeof(double*) * size_of_Tb ); 
    int **Ep = (int **)malloc( sizeof(int*) * size_of_Edgep );                                                
    int **Eb = (int **)malloc( sizeof(int*) * size_of_Edgeb );                                                
    int **Ep2 = (int **)malloc( sizeof(int*) * size_of_Edgep );                                               
    int **Eb2 = (int **)malloc( sizeof(int*) * size_of_Edgeb );  
    for( i=0; i < size_of_Vp; i++ )                                                                           
    {                                                                                                         
        Vp[i] = (double *)malloc( sizeof(double) * 3 );                                                       
    }                                                                                                         
                                                                                                              
    for( i=0; i < size_of_Vb; i++ )                                                                           
    {                                                                                                         
        Vb[i] = (double *)malloc( sizeof(double) * 3 );                                                       
    }                                                                                                         
                                                                                                              
    for( i=0; i < size_of_Dp; i++ )                                                                           
    {                                                                                                         
        Dp[i] = (int *)malloc( sizeof(int) * 3 );                                                             
        Tp[i] = (double *)malloc( sizeof(double) * 9 );                                                       
    }                                                                                                         
                                                                                                              
    for( i=0; i < size_of_Db; i++ )                                                                           
    {                                                                                                         
        Db[i] = (int *)malloc( sizeof(int) * 3 );
        Tb[i] = (double *)malloc( sizeof(double) * 9 );
    }                                                                                                         
                                                                                                              
    for( i=0; i < size_of_Dq; i++ )                                                                           
    {                                                                                                         
        Dq[i] = (int *)malloc( sizeof(int) * 4 );                                                             
    }                                                                                                         
                                                                                                              
    for( i=0; i < size_of_Edgep; i++ )                                                                        
    {                                                                                                         
        Ep[i] = (int *)malloc( sizeof(int) * 2 );                                                             
    }                                                                                                         
    for( i=0; i < size_of_Edgeb; i++ )                                                                        
    {                                                                                                         
        Eb[i] = (int *)malloc( sizeof(int) * 2 );                                                             
    } 
    for( i=0; i < size_of_Edgep; i++ )                                                                        
    {                                                                                                         
        Ep2[i] = (int *)malloc( sizeof(int) * 2 );                                                            
    }                                                                                                         
    for( i=0; i < size_of_Edgeb; i++ )                                                                        
    {                                                                                                         
        Eb2[i] = (int *)malloc( sizeof(int) * 2 );                                                            
    } 
    for( i=0; i < size_of_Vp; i++ )                                                                          
    {                                                                                                         
        for( j=0; j < 3; j++ )                                                                                
        {                                                                                                     
            Vp[i][j] = PyVp[i*3+j];                                                                           
        }                                                                                                     
    }                                                                                                         
                                                                                                              
    for( i=0; i < size_of_Vb; i++ )                                                                           
    {                                                                                                         
        for( j=0; j < 3; j++ )                                                                                
        {                                                                                                     
            Vb[i][j] = PyVb[i*3+j];                                                                           
        }                                                                                                     
    }                                                                                                         
                                                                                                              
    for( i=0; i < size_of_Dp; i++ )                                                                           
    {                                                                                                         
        for( j=0; j < 3; j++ )                                                                                
        {                                                                                                     
            Dp[i][j] = PyDp[i*3+j];                                                                           
            for( k = 0; k < 3; k++ )                                                                          
            {                                                                                                 
                Tp[i][3*j+k] = PyTp[i*9+3*j+k];                                                               
            }                                                                                                 
        }                                                                                                     
    }                                                                                                         
                                                                                                              
    for( i=0; i < size_of_Db; i++ )                                                                           
    {                                                                                                         
        for( j=0; j < 3; j++ )                                                                                
        {                                                                                                     
            Db[i][j] = PyDb[i*3+j];                                                                           
            if( Db[i][j]> ndofb )                                                                             
            {                                                                                                 
                ndofb = Db[i][j];                                                                             
            }
            for( k = 0; k < 3; k++ )                                                                          
            {                                                                                                 
                Tb[i][3*j+k] = PyTb[i*9+3*j+k];                                                               
            }                                                                                                  
        }                                                                                                     
    }                                                                                                         
                                                                                                              
    for( i=0; i < size_of_Dq; i++ )                                                                           
    {                                                                                                         
        for( j=0; j < 4; j++ )                                                                                
        {                                                                                                     
            Dq[i][j] = PyDq[i*4+j];                                                                           
        }                                                                                                     
    }                                                                                                         
                                                                                                              
    for( i=0; i < size_of_Edgep; i++ )                                                                        
    {                                                                                                         
        for( j=0; j < 2; j++ )                                                                                
        {                                                                                                     
            Ep[i][j] = PyEdgep[i*2+j];                                                                        
        }                                                                                                     
    }
    for( i=0; i < size_of_Edgeb; i++ )                                                                        
    {                                                                                                         
        for( j=0; j < 2; j++ )                                                                                
        {                                                                                                     
            Eb[i][j] = PyEdgeb[i*2+j];                                                                        
        }                                                                                                     
    }
    for( i=0; i < size_of_Edgep; i++ )                                                                        
    {                                                                                                         
        for( j=0; j < 2; j++ )                                                                                
        {                                                                                                     
            Ep2[i][j] = PyEdgep2[i*2+j];                                                                      
        }                                                                                                     
    }                                                                                                         
    for( i=0; i < size_of_Edgeb; i++ )                                                                        
    {                                                                                                         
        for( j=0; j < 2; j++ )                                                                                
        {                                                                                                     
            Eb2[i][j] = PyEdgeb2[i*2+j];                                                                      
        }                                                                                                     
    } 
    M->Vp = Vp;
    M->Vb = Vb;
    M->Tp = Tp;
    M->Tb = Tb;
    M->Dp = Dp;
    M->Db = Db;
    M->Dq = Dq;
    M->Ep = Ep;
    M->Eb = Eb;
    M->Ep2 = Ep2;                                                                                               
    M->Eb2 = Eb2;
    M->ntp = size_of_Tp;                                                                                   
    M->ndp = size_of_Dp;                                                                                   
    M->ndb = size_of_Db;                                                                                   
    M->nvp = size_of_Vp;                                                                                   
    M->nvb = size_of_Vb;                                                                                   
    M->nep = size_of_Edgep;                                                                                
    M->neb = size_of_Edgeb;                                                                                
    M->ndofb = ndofb;
    return M;
}

pSpace new_space2(pMesh M, char *type)
{
    int opt;

    if( strcmp(type, "P0" ) == 0)                                                                           
    {                                                                                                         
        opt = 0;                                                                                       
    }                                                                                                         
    if( strcmp(type, "P1" ) == 0 )                                                                          
    {                                                                                                         
        opt = 1;                                                                                       
    }                                                                                                         
    if( strcmp(type, "P0d" ) == 0)                                                                          
    {                                                                                                         
        opt = 2;                                                                                       
    }                                                                                                         
    if( strcmp(type, "P1d" ) == 0 )                                                                         
    {                                                                                                         
        opt = 3;                                                                                       
    }                                                                                                         
    if( strcmp(type, "P0b" ) == 0 )                                                                         
    {                                                                                                         
        opt = 4;                                                                                       
    }                                                                                                         
    if( strcmp(type, "P1b" ) == 0 )                                                                         
    {                                                                                                         
        opt = 5;                                                                                       
    }
    if( strcmp(type, "RWG" ) == 0 )                                                                           
    {                                                                                                         
        opt = 6;                                                                                              
    }
    if( strcmp(type, "BC" ) == 0 )                                                                           
    {                                                                                                         
        opt = 7;                                                                                              
    }
    pSpace S;
    S = new_Space(M->Vp, M->Vb, M->Dp, M->Db, M->Dq, M->Tp, M->Tb, M->Ep, M->Eb, M->Ep2, M->Eb2, M->ndp, M->ndb, M->nvp, M->nvb, M->nep, M->neb, opt);
    S->M = M;
    return S;

}
psupermatrix hierarchical_mat_slp (pMeshes Mesh, double eta, int ls, double eps, int *order, double kappa, double *time)
{
    struct timeval st, e; 
    psupermatrix s;
    gettimeofday(&st, NULL);
    s = create_hmat(Mesh, eta, ls, eps, order, kappa);
    gettimeofday(&e, NULL);
    time[0] = (e.tv_sec - st.tv_sec) * 1000.0;      // sec to ms                                            
    time[0] += (e.tv_usec - st.tv_usec) / 1000.0;
    
    return s; 
}

psupermatrix hierarchical_mat_hlp (pMeshes Mesh, double eta, int ls, double eps, int *order, double kappa, double *time)
{
    struct timeval st, e;                                                                                     
    psupermatrix s;
    gettimeofday(&st, NULL);                                                                              
    s = create_hmat_hlp(Mesh, eta, ls, eps, order, kappa);                                       
    gettimeofday(&e, NULL);                                                                               
    time[0] = (e.tv_sec - st.tv_sec) * 1000.0;      // sec to ms                                            
    time[0] += (e.tv_usec - st.tv_usec) / 1000.0;
    
    return s;
}

double *get_mat(int N)
{
    double *M;
    M = (double*)malloc(N * N * sizeof(double));
    return M;
}

double get_element( double *M, int N, int row, int col)
{
    return M[row*N+col];
}

double get_gsl_real( gsl_matrix_complex *M, int row, int col )
{
    return GSL_REAL(gsl_matrix_complex_get(M, row, col)); 
}

double get_gsl_imag( gsl_matrix_complex *M, int row, int col )                                                
{                                                                                                             
   return GSL_IMAG(gsl_matrix_complex_get(M, row, col));                                                           
}

void del_gsl_mat(gsl_matrix_complex *M)
{
    gsl_matrix_complex_free(M);
}

void del_hierarchical_mat(psupermatrix s)
{
    free(s->index);
    del_supermatrix(s);
}

void del_potential(double *MR, double *MI)
{
    free(MR);
    free(MI);
}

void del_id(double *M)
{
    free(M);
}

void supermat2supermataux(psupermatrix s, psupermataux sa)
{
    int block_cols = s->block_cols;
    int block_rows = s->block_rows;
    int i,j;

    if(block_cols*block_rows>1)
    {
        for(j=0; j<block_cols; j++)
        {
            for(i=0; i<block_rows; i++)
            {
                supermat2supermataux( s->s[i+j*block_rows], sa);   
            }
        }
    }
    else
    {
        sa -> s = (psupermatrix*) realloc(sa -> s,(sa->size+1)*sizeof(psupermatrix));
        sa -> s[sa->size] = s; 
        sa->size = sa->size+1;
    }

}

void supermat2supermataux2(psupermatrix s, psupermataux sa, double complex *in, double complex *out)                                                    
{                                                                                                             
    int block_cols = s->block_cols;                                                                           
    int block_rows = s->block_rows;                                                                           
    int i,j;                                                                                                  
                                                                                                              
    if(block_cols*block_rows>1)                                                                               
    {                                                                                                         
        for(j=0; j<block_cols; j++)                                                                           
        {                                                                                                     
            for(i=0; i<block_rows; i++)                                                                       
            {                                                                                                 
                supermat2supermataux2( s->s[i+j*block_rows], sa, in, out);                                              
            }                                                                                                 
        }                                                                                                     
    }                                                                                                         
    else                                                                                                      
    {                                                                                                         
        sa -> s = (psupermatrix*) realloc(sa -> s,(sa->size+1)*sizeof(psupermatrix));                         
        sa -> s[sa->size] = s;                                                                                
        sa->s[sa->size]->vecin = (double complex **) malloc(sizeof(double complex*) * s->cols->M->PD1->nt);                                          
        sa->s[sa->size]->vecout = (double complex **) malloc(sizeof(double complex*) * s->rows->M->PD1->nt);
        /*for(i = 0; i<s->cols->M->PD1->nt; i++)                                                               
        { 
            sa->s[sa->size]->vecin[i] = &in[s->cols->M->PD1->triangles[i]];
            //printf("%d\n", &in[s->cols->M->PD1->triangles[i]]);
        }
        for(i = 0; i<s->rows->M->PD1->nt; i++)                                                               
        { 
            sa->s[sa->size]->vecout[i] = &out[s->rows->M->PD1->triangles[i]];                                                                                
        } */
        sa->size = sa->size+1;
    }                                                                                                         
                                                                                                              
} 

void del_supermataux(psupermataux s)
{
    int i;
    for(i=0;i<s->size;i++)
    {
        s->s[i] = NULL;
        free(s->s[i]);
    }
    free(s->index);
    free(s);
}

/*double *gmres_prec_hmat( double *MR, double *MI, double *MMR, double *MMI, psupermatrix P, double *br, double *bi, double *xr, double *xi,  double tol, int max_iter, int size_of_b, double *tm )
{
    
    int                    i, j ;
    int                 iter = 0;
    int                   status;
    double *r =calloc(size_of_b ,sizeof(double));
    gsl_complex             aux;
    struct timeval st, e; 
    
    const gsl_linalg_itersolve_type_c gmres_type =
    {
        "gmres",
        &gmres_alloc_c,
        &gmres_iterate_c,
        &gmres_normr_c,
        &gmres_free_c
    };
    const gsl_linalg_itersolve_type_c * gsl_linalg_itersolve_gmres =
    &gmres_type;
    
    psupermataux sa = (psupermataux*) malloc(sizeof(psupermataux));
    sa -> s = (psupermatrix*) malloc(sizeof(psupermatrix));
    sa->size = 0;
    supermat2supermataux(P, sa);
    sa->index = (int *) malloc((size_t) sizeof(int) * P->nv);
    sa->index = P->index;

    gsl_matrix_complex *M = gsl_matrix_complex_alloc(size_of_b, size_of_b);
    gsl_matrix_complex *A = gsl_matrix_complex_alloc(size_of_b, size_of_b);
    gsl_matrix_complex *MM = gsl_matrix_complex_alloc(size_of_b, size_of_b);
    gsl_vector_complex *b = gsl_vector_complex_alloc (size_of_b);
    gsl_vector_complex *baux = gsl_vector_complex_alloc(size_of_b);
    gsl_vector_complex *x = gsl_vector_complex_alloc (size_of_b);
    
    for( i=0; i < size_of_b; i++ )
    {
        for( j=0; j < size_of_b; j++ )
        {
            GSL_SET_COMPLEX(&aux, MR[i*size_of_b+j], MI[i*size_of_b+j]);
            gsl_matrix_complex_set(M, i, j, aux);
            GSL_SET_COMPLEX(&aux, MMR[i*size_of_b+j], MMI[i*size_of_b+j]);
            gsl_matrix_complex_set(MM, i, j, aux);
        }
    }

    for( i=0; i < size_of_b; i++ )
    {
        GSL_SET_COMPLEX(&aux, br[i], bi[i]);
        gsl_vector_complex_set(b, i ,aux);
        GSL_SET_COMPLEX(&aux, xr[i], xi[i]);                                                                  
        gsl_vector_complex_set(x, i ,aux);
    }
    
    const gsl_linalg_itersolve_type_c *T = gsl_linalg_itersolve_gmres;
    gsl_vector_complex *u = gsl_vector_complex_calloc(size_of_b);   
    gsl_linalg_itersolve_c *w;
    w = calloc(1, sizeof(gsl_linalg_itersolve_c));
    gsl_linalg_itersolve_alloc_c(T, w, size_of_b, size_of_b);

    do
    {
        gettimeofday(&st, NULL);
        gsl_blas_zgemm (CblasNoTrans, CblasNoTrans, gsl_complex_rect(1,0),MM, M, gsl_complex_rect(0,0),A);
        mat_mat(sa, A, tm);
        gsl_blas_zgemm (CblasNoTrans, CblasNoTrans, gsl_complex_rect(1,0),MM, A, gsl_complex_rect(0,0),M); 
        gsl_blas_zgemv (CblasNoTrans, gsl_complex_rect(1,0), MM, b, gsl_complex_rect(0,0), baux);
        mat_vec(sa, baux, b, tm);
        gsl_blas_zgemv (CblasNoTrans, gsl_complex_rect(1,0), MM, b, gsl_complex_rect(0,0), baux);
        gettimeofday(&e, NULL);                                                                               
        tm[2] = (e.tv_sec - st.tv_sec) * 1000.0;      // sec to ms                                            
        tm[2] += (e.tv_usec - st.tv_usec) / 1000.0;
        gettimeofday(&st, NULL); 
        status = gmres_iterate_c(M, NULL, 0, baux, tol, x, w->state,r); 
        gettimeofday(&e, NULL);
        tm[3] = (e.tv_sec - st.tv_sec) * 1000.0;      // sec to ms
        tm[3] += (e.tv_usec - st.tv_usec) / 1000.0;

    }
    while (status == 0 && ++iter<max_iter);
    
    for( i = 0; i<size_of_b; i++)
    {
        aux = gsl_vector_complex_get(x,i);
        xr[i] =GSL_REAL(aux);
        xi[i] =GSL_IMAG(aux);
    }
    
    gsl_matrix_complex_free(M);
    gsl_matrix_complex_free(A);
    gsl_matrix_complex_free(MM);
    gsl_vector_complex_free(b);
    gsl_vector_complex_free(baux);
    gsl_vector_complex_free(x);
    
    return r;
}

double *gmres_prec( double *MR, double *MI, double *MMR, double *MMI, double *PR, double *Pi, double *br, double *bi, double *xr, double *xi,  double tol, int max_iter, int size_of_b, double *tm )
{
    
    int                    i, j ;
    int                 iter = 0;
    int                   status;
    double *r =calloc(size_of_b ,sizeof(double));
    gsl_complex             aux;
    struct timeval st, e; 

    const gsl_linalg_itersolve_type_c gmres_type =
    {
        "gmres",
        &gmres_alloc_c,
        &gmres_iterate_c,
        &gmres_normr_c,
        &gmres_free_c
    };
    const gsl_linalg_itersolve_type_c * gsl_linalg_itersolve_gmres =
    &gmres_type;
    
    
    gsl_matrix_complex *M = gsl_matrix_complex_alloc(size_of_b, size_of_b);
    gsl_matrix_complex *MM = gsl_matrix_complex_alloc(size_of_b, size_of_b);
    gsl_matrix_complex *P = gsl_matrix_complex_alloc(size_of_b, size_of_b);
    gsl_vector_complex *b = gsl_vector_complex_alloc (size_of_b);
    gsl_vector_complex *x = gsl_vector_complex_alloc (size_of_b);    
    for( i=0; i < size_of_b; i++ )
    {
        for( j=0; j < size_of_b; j++ )
        {
            GSL_SET_COMPLEX(&aux, MR[i*size_of_b+j], MI[i*size_of_b+j]);
            gsl_matrix_complex_set(M, i, j, aux);
            GSL_SET_COMPLEX(&aux, MMR[i*size_of_b+j], MMI[i*size_of_b+j]);
            gsl_matrix_complex_set(MM, i, j, aux);
            GSL_SET_COMPLEX(&aux, PR[i*size_of_b+j], Pi[i*size_of_b+j]);
            gsl_matrix_complex_set(P, i, j, aux);
        }
    }
    
    
    
    for( i=0; i < size_of_b; i++ )
    {
        GSL_SET_COMPLEX(&aux, br[i], bi[i]);
        gsl_vector_complex_set(b, i ,aux);
        GSL_SET_COMPLEX(&aux, xr[i], xi[i]);                                                                  
        gsl_vector_complex_set(x, i ,aux);
    }
    
    const gsl_linalg_itersolve_type_c *T = gsl_linalg_itersolve_gmres;
    gsl_vector_complex *u = gsl_vector_complex_calloc(size_of_b);       
    gsl_linalg_itersolve_c *w;
    w = calloc(1, sizeof(gsl_linalg_itersolve_c));
    gsl_linalg_itersolve_alloc_c(T, w, size_of_b, size_of_b);
    do
    {
        gettimeofday(&st, NULL);
        status = gmres_iterate_c(M, P, 1, b, tol, x, w->state,r);
        gettimeofday(&e, NULL);                                                                                                                                    
        tm[0] = (e.tv_sec - st.tv_sec) * 1000.0;      // sec to ms                                            
        tm[0] += (e.tv_usec - st.tv_usec) / 1000.0;
        
    }
    while (status == 0 && ++iter<max_iter);
    
    for( i = 0; i<size_of_b; i++)
    {
        aux = gsl_vector_complex_get(x,i);
        xr[i] =GSL_REAL(aux);
        xi[i] =GSL_IMAG(aux);
    }
    
    gsl_matrix_complex_free(M);
    gsl_matrix_complex_free(P);
    gsl_matrix_complex_free(MM);
    gsl_vector_complex_free(b);
    gsl_vector_complex_free(x);
    return r;
}

double *gmres( double *MR, double *MI, double *br, double *bi, double *xr, double *xi,  double tol, int max_iter, int size_of_b, double *tm )
{
    
    int                    i, j ;
    int                 iter = 0;
    int                   status;
    double *r =calloc(size_of_b ,sizeof(double));
    gsl_complex             aux;
    struct timeval st, e;
    
    
    const gsl_linalg_itersolve_type_c gmres_type =
    {
        "gmres",
        &gmres_alloc_c,
        &gmres_iterate_c,
        &gmres_normr_c,
        &gmres_free_c
    };
    const gsl_linalg_itersolve_type_c * gsl_linalg_itersolve_gmres = &gmres_type;
    
    
    gsl_matrix_complex *M = gsl_matrix_complex_alloc(size_of_b, size_of_b);
    gsl_vector_complex *b = gsl_vector_complex_alloc (size_of_b);
    
    for( i=0; i < size_of_b; i++ )
    {
        for( j=0; j < size_of_b; j++ )
        {
            GSL_SET_COMPLEX(&aux, MR[i*size_of_b+j], MI[i*size_of_b+j]);
            gsl_matrix_complex_set(M, i, j, aux);
        }
    }
    
    for( i=0; i < size_of_b; i++ )
    {
        GSL_SET_COMPLEX(&aux, br[i], bi[i]);
        gsl_vector_complex_set(b, i ,aux);
    }
    
    const gsl_linalg_itersolve_type_c *T = gsl_linalg_itersolve_gmres;
    gsl_vector_complex *u = gsl_vector_complex_calloc(size_of_b);   
    gsl_linalg_itersolve_c *w;
    w = calloc(1, sizeof(gsl_linalg_itersolve_c));
    gsl_linalg_itersolve_alloc_c(T, w, size_of_b, size_of_b);
    do
    {
        gettimeofday(&st, NULL);
        status = gmres_iterate_c(M, NULL, 0, b, tol, u, w->state,r);
        gettimeofday(&e, NULL);
        tm[0] = (e.tv_sec - st.tv_sec) * 1000.0;      // sec to ms                                            
        tm[0] += (e.tv_usec - st.tv_usec) / 1000.0; 
    }
    while (status == 0 && ++iter<max_iter);
    
    for( i = 0; i<size_of_b; i++)
    {
        aux = gsl_vector_complex_get(u,i);
        xr[i] =GSL_REAL(aux);
        xi[i] =GSL_IMAG(aux);
    }
    
    gsl_matrix_complex_free(M);
    gsl_vector_complex_free(b);
    
    return r;
}
*/

void helmholtz_assemble_matrix_sl (double *MR, double *MI, pMeshes Mesh, double kappa, int *order, int opt, double *time)
{
    memset(MR, 0, sizeof(*MR));
    memset(MI, 0, sizeof(*MI));
    struct timeval st, e;
    if (opt == 0)
    {
        gettimeofday(&st, NULL); 
        assemble_SLP( Mesh, kappa, order, MR, MI, Mesh->trialbf, Mesh->testbf);
        gettimeofday(&e, NULL);                                                                                                                    
        time[0] = (e.tv_sec - st.tv_sec) * 1000.0;      // sec to ms                                              
        time[0] += (e.tv_usec - st.tv_usec) / 1000.0;
    }
    else if (opt == 1)
    {
        gettimeofday(&st, NULL); 
        assemble_SLD( Mesh, kappa, order, MR, MI, Mesh->trialbf, Mesh->testbf);
        gettimeofday(&e, NULL);                                           
        time[0] = (e.tv_sec - st.tv_sec) * 1000.0;      // sec to ms                                              
        time[0] += (e.tv_usec - st.tv_usec) / 1000.0;
    }
    

}

void helmholtz_assemble_matrix_hl (double *MR, double *MI, pMeshes Mesh, double kappa, int *order, double *time, int opt)
{                                                                           
    struct timeval st, e;
    memset(MR, 0, sizeof(*MR));
    memset(MI, 0, sizeof(*MI));
    if(opt==0)
    {
        gettimeofday(&st, NULL);
        assemble_HP( Mesh,kappa, order, MR, MI, Mesh->trialbf, Mesh->testbf);
        gettimeofday(&e, NULL);
        time[0] = (e.tv_sec - st.tv_sec) * 1000.0;      // sec to ms                                              
        time[0] += (e.tv_usec - st.tv_usec) / 1000.0;
    }
    else
    {
        gettimeofday(&st, NULL);                                                                              
        assemble_HP2( Mesh,kappa, order, MR, MI, Mesh->trialbf, Mesh->testbf);                                 
        gettimeofday(&e, NULL);                                                                               
        time[0] = (e.tv_sec - st.tv_sec) * 1000.0;      // sec to ms                                              
        time[0] += (e.tv_usec - st.tv_usec) / 1000.0; 
    }
    
    
}

void helmholtz_assemble_matrix_id (double *MR, pMeshes Mesh,  int order, int opt, double *time)
{
    memset(MR, 0, sizeof(*MR));                                                                                     
    struct timeval st, e;
    if (opt == 0)
    {
        gettimeofday(&st, NULL);
        assemble_IDP( Mesh, order, MR, Mesh->trialbf, Mesh->testbf);
        gettimeofday(&e, NULL);
        time[0] = (e.tv_sec - st.tv_sec) * 1000.0;      // sec to ms                                              
        time[0] += (e.tv_usec - st.tv_usec) / 1000.0;
    }
    else if (opt == 1)
    {
        gettimeofday(&st, NULL);
        assemble_IDD( Mesh, order, MR, Mesh->trialbf, Mesh->testbf);
        gettimeofday(&e, NULL);
        time[0] = (e.tv_sec - st.tv_sec) * 1000.0;      // sec to ms                                              
        time[0] += (e.tv_usec - st.tv_usec) / 1000.0; 
    }
    else if (opt == 2)
    {
        gettimeofday(&st, NULL);
        assemble_IDDB( Mesh, order, MR, Mesh->trialbf, Mesh->testbf );
        gettimeofday(&e, NULL);
        time[0] = (e.tv_sec - st.tv_sec) * 1000.0;      // sec to ms                                              
        time[0] += (e.tv_usec - st.tv_usec) / 1000.0;
    }
    else if (opt == 3)
    {
        gettimeofday(&st, NULL);
        assemble_IDM( Mesh, order, MR, Mesh->trialbf, Mesh->testbf);
        gettimeofday(&e, NULL);
        time[0] = (e.tv_sec - st.tv_sec) * 1000.0;      // sec to ms                                              
        time[0] += (e.tv_usec - st.tv_usec) / 1000.0;
    }
    else if (opt == 4)
    {
        gettimeofday(&st, NULL);
        assemble_IDDP( Mesh, order, MR, Mesh->trialbf, Mesh->testbf);
        gettimeofday(&e, NULL);
        time[0] = (e.tv_sec - st.tv_sec) * 1000.0;      // sec to ms                                              
        time[0] += (e.tv_usec - st.tv_usec) / 1000.0;
    }
    else if (opt == 5)                                                                                        
    {                                                                                                         
        gettimeofday(&st, NULL);                                                                              
        projection_matrix( Mesh, MR, 0 );                                                                   
        gettimeofday(&e, NULL);                                                                               
        time[0] = (e.tv_sec - st.tv_sec) * 1000.0;      // sec to ms                                              
        time[0] += (e.tv_usec - st.tv_usec) / 1000.0;                                                         
    }
    else if (opt == 6)
    {
        gettimeofday(&st, NULL);
        average_matrix_full( Mesh, MR, 1 );
        gettimeofday(&e, NULL);
        time[0] = (e.tv_sec - st.tv_sec) * 1000.0;      // sec to ms                                              
        time[0] += (e.tv_usec - st.tv_usec) / 1000.0;
    }


}

void helmholtz_project_mat(double *MRb, double *MIb, double *MR, double *MI,pMeshes Mesh, double *time)
{
    struct timeval st, e;                                                                                     
    memset(MR, 0, sizeof(*MR));                                                                               
    memset(MI, 0, sizeof(*MI));
    gettimeofday(&st, NULL);
    project_mat(Mesh,MRb, MIb, MR, MI,Mesh->trialbf, Mesh->testbf);
    gettimeofday(&e, NULL);                                                                               
    time[0] = (e.tv_sec - st.tv_sec) * 1000.0;      // sec to ms                                              
    time[0] += (e.tv_usec - st.tv_usec) / 1000.0;
}

void maxwell_weakly(double *MR, double *MI, pMeshes Mesh, double kappa, int *order, double *time, int opt)
{                                                                                                             
    struct timeval st, e;                                                                                     
    memset(MR, 0, sizeof(*MR));                                                                               
    memset(MI, 0, sizeof(*MI));                                                                                                                                                                      
    gettimeofday(&st, NULL);     
    assemble_maxwell_weakly( Mesh,kappa, order, MR, MI, Mesh->trialbf, Mesh->testbf);                                 
    gettimeofday(&e, NULL);                                                                               
    time[0] = (e.tv_sec - st.tv_sec) * 1000.0;      // sec to ms                                              
    time[0] += (e.tv_usec - st.tv_usec) / 1000.0; 

}

void maxwell_hypersingular(double *MR, double *MI, pSpace space_trial, pSpace space_test, double kappa, int *order, double *time, int opt)    
{                                                                                                             
    struct timeval st, e;                                                                                     
    memset(MR, 0, sizeof(*MR));                                                                               
    memset(MI, 0, sizeof(*MI));                                                                               
    gettimeofday(&st, NULL);                                                                                  
    assemble_maxwell_hypersingular( space_trial, space_test, kappa, order, MR, MI);                         
    gettimeofday(&e, NULL);                                                                                   
    time[0] = (e.tv_sec - st.tv_sec) * 1000.0;      // sec to ms                                              
    time[0] += (e.tv_usec - st.tv_usec) / 1000.0;                                                             
                                                                                                              
}

void maxwell_T(double *MR1, double *MI1, double *MR2, double *MI2, pSpace space_trial, pSpace space_test, double kappa, int *order, double *time, int opt)    
{                                                                                                             
    struct timeval st, e;                                                                                     
    memset(MR1, 0, sizeof(*MR1));                                                                               
    memset(MI1, 0, sizeof(*MI1));
    memset(MR2, 0, sizeof(*MR2));                                                                             
    memset(MI2, 0, sizeof(*MI2));
    gettimeofday(&st, NULL);                                                                                  
    assemble_maxwell_T( space_trial, space_test ,kappa, order, MR1, MI1, MR2, MI2 );                         
    gettimeofday(&e, NULL);                                                                                   
    time[0] = (e.tv_sec - st.tv_sec) * 1000.0;      // sec to ms                                              
    time[0] += (e.tv_usec - st.tv_usec) / 1000.0;                                                             
                                                                                                              
} 
void maxwell_id(double *MR, pMeshes Mesh, double kappa, int *order, double *time, int opt)           
{                                                                                                             
    struct timeval st, e;                                                                                     
    memset(MR, 0, sizeof(*MR));                                                                                                                                                           
    gettimeofday(&st, NULL);                                                                                  
    assemble_maxwell_identity( Mesh,kappa, order, MR, Mesh->trialbf, Mesh->testbf);                                
    gettimeofday(&e, NULL);                                                                                   
    time[0] = (e.tv_sec - st.tv_sec) * 1000.0;      // sec to ms                                              
    time[0] += (e.tv_usec - st.tv_usec) / 1000.0;                                                             
                                                                                                              
}

void maxwell_lb(double *MR, pMeshes Mesh, double kappa, int *order, double *time, int opt)                    
{                                                                                                             
    struct timeval st, e;                                                                                     
    memset(MR, 0, sizeof(*MR));                                                                               
    gettimeofday(&st, NULL);                                                                                  
    assemble_maxwell_laplace_beltrami( Mesh,kappa, order, MR, Mesh->trialbf, Mesh->testbf);                           
    gettimeofday(&e, NULL);                                                                                   
    time[0] = (e.tv_sec - st.tv_sec) * 1000.0;      // sec to ms                                              
    time[0] += (e.tv_usec - st.tv_usec) / 1000.0;                                                             
                                                                                                              
}

void mv(psupermatrix P, double *br, double *bi, double *t, int size_of_b)                             
{   
    struct timeval st, e;
    gsl_vector_complex *b = gsl_vector_complex_alloc (size_of_b);
    gsl_vector_complex *baux = gsl_vector_complex_calloc (size_of_b);
    gsl_complex aux;
    int i;                                                                              
    for( i=0; i < size_of_b; i++ )                                                                            
    {                                                                                                         
        GSL_SET_COMPLEX(&aux, br[i], bi[i]);                                                                  
        gsl_vector_complex_set(b, i ,aux);
    }
    psupermataux sa = (psupermataux) malloc(sizeof(supermataux));                                           
    sa -> s = (psupermatrix*) malloc(sizeof(psupermatrix));                                                   
    sa->size = 0;                                                                                             
    supermat2supermataux(P, sa);                                                                              
    sa->index = (int *) malloc((size_t) sizeof(int) * P->nv);                                                 
    sa->index = P->index;
    gettimeofday(&st, NULL);
    mat_vec(sa, b, baux, t);
    gettimeofday(&e, NULL);
    t[2] = (e.tv_sec - st.tv_sec) * 1000.0;      // sec to ms                                              
    t[2] += (e.tv_usec - st.tv_usec) / 1000.0;
    for( i=0; i < size_of_b; i++ )                                                                            
    {                                                                                                                                                          
        aux = gsl_vector_complex_get(baux, i); 
        br[i] = GSL_REAL(aux); 
        bi[i] = GSL_IMAG(aux);                                                                  
    } 
    gsl_vector_complex_free(b);
    gsl_vector_complex_free(baux);
}

gsl_matrix_complex *mtmt(psupermatrix P, double *MR, double *MI, double *t, int size_of_b)
{
    gsl_matrix_complex *M = gsl_matrix_complex_alloc(size_of_b, size_of_b);                                   
    gsl_complex aux;                                  
    struct timeval st, e;                                                                                                           
    int i,j;
    for( i=0; i < size_of_b; i++ )                                                                            
    {                                                                                                         
        for( j=0; j < size_of_b; j++ )                                                                        
        {                                                                                                     
            GSL_SET_COMPLEX(&aux, MR[i*size_of_b+j], MI[i*size_of_b+j]);                                      
            gsl_matrix_complex_set(M, i, j, aux);                                                             
        }                                                                                                     
    } 
    psupermataux sa = (psupermataux) malloc(sizeof(supermataux));                                           
    sa -> s = (psupermatrix*) malloc(sizeof(psupermatrix));                                                   
    sa->size = 0;                                                                                             
    supermat2supermataux(P, sa);                                                                              
    sa->index = (int *) malloc((size_t) sizeof(int) * P->nv);                                                 
    sa->index = P->index;
    gettimeofday(&st, NULL);
    mat_mat(sa, M, t);
    gettimeofday(&e, NULL);                                                                                   
    t[2] = (e.tv_sec - st.tv_sec) * 1000.0;      // sec to ms                                              
    t[2] += (e.tv_usec - st.tv_usec) / 1000.0; 
    for( i=0; i < size_of_b; i++ )                                                                            
    {                                                                                                         
        for( j=0; j < size_of_b; j++ )                                                                        
        {                                                                                                                      
            aux = gsl_matrix_complex_get(M, i, j);
            MR[i*size_of_b+j] = GSL_REAL(aux);
            MI[i*size_of_b+j] = GSL_IMAG(aux);  
        }                                                                                                     
    }
    del_supermataux(sa);
    //gsl_matrix_complex_free(M);
    return M;
}


void mv2(psupermatrix P, double *br, double *bi, double *t, int size_of_b)                                     
{                                                                                                             
    struct timeval st, e;                                                                                     
    gsl_vector_complex *b = gsl_vector_complex_alloc (size_of_b);                                             
    gsl_vector_complex *baux = gsl_vector_complex_calloc (size_of_b); 
    //double complex *in = (double *) malloc(sizeof(double) * size_of_b); 
    //double complex *out = (double *) calloc(size_of_b,sizeof(double));
    gsl_complex aux;                                                                                          
    int i;                                                                                                    
    for( i=0; i < size_of_b; i++ )                                                                            
    {                                                                                                         
        GSL_SET_COMPLEX(&aux, br[i], bi[i]);                                                                  
        gsl_vector_complex_set(b, i ,aux);
        //in[i] = br[i]+I*bi[i];
    }                                                                                                         
    psupermataux sa = (psupermataux) malloc(sizeof(supermataux));                                           
    sa -> s = (psupermatrix*) malloc(sizeof(psupermatrix));                                                   
    sa->size = 0;                                                                                             
    supermat2supermataux(P, sa);                                                                              
    sa->index = (int *) malloc((size_t) sizeof(int) * P->nv);                                                 
    sa->index = P->index;                                                                                     
    gettimeofday(&st, NULL);                                                                                  
    mat_vec2(sa, b, baux, t);                                                                                  
    gettimeofday(&e, NULL);                                                                                   
    t[2] = (e.tv_sec - st.tv_sec) * 1000.0;      // sec to ms                                                 
    t[2] += (e.tv_usec - st.tv_usec) / 1000.0;                                                                
    for( i=0; i < size_of_b; i++ )                                                                            
    {                                                                                                         
        aux = gsl_vector_complex_get(baux, i);                                                                
        br[i] = GSL_REAL(aux);                                                                                
        bi[i] = GSL_IMAG(aux);                                                                                
    }                                                                                                         
    gsl_vector_complex_free(b);                                                                               
    gsl_vector_complex_free(baux);                                                                            
}
