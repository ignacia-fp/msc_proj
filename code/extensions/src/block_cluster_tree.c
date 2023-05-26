
//  Block_cluster_tree.c
//  Block cluster tree
//
//  Created by María Ignacia Fierro on 11-12-17.
//  Copyright © 2017 María Ignacia Fierro. All rights reserved.
//

#include "block_cluster_tree.h"



double dnrm2_(int *n,double *x, int* incx);


/*****************************************************/
double distance_cluster(pcluster row, pcluster col)
{
    int ione = 1;
    int i;
    double *bmin_row, *bmax_row, *bmin_col, *bmax_col;
    int sz =2;
    int d = row -> d;
    double res = 0;
    //double *z = (double*) malloc(2*sizeof(double));
    bmin_row = row ->bmin;
    bmax_row = row ->bmax;
    bmin_col = col ->bmin;
    bmax_col = col ->bmax;
    double n;
    
    for( i = 0; i<d; i++ )
    {
        n =MIN(fabs(bmin_row[i]-bmin_col[i]),fabs(bmin_row[i]-bmax_col[i]));
        n =MIN(n,fabs(bmax_row[i]-bmin_col[i]));
        n =MIN(n,fabs(bmax_row[i]-bmax_col[i]));
        res+=n*n;
        //z[0] = bmin_row[i]-bmin_col[i];
        //z[1] = bmax_row[i]-bmax_col[i];
        //res+=dnrm2_(&sz,z,&ione);

    }
    res = sqrt(res);

    return res;
}

double diameter_cluster(pcluster c)
{
    double *bmin, *bmax;
    int i, d;
    bmin = c -> bmin;
    bmax = c -> bmax;
    d = c -> d;
    double diam = 0.0;
    
    for (i = 0; i < d; i++)
    {
        if(abs(bmax[i]-bmin[i])>diam)
        {
            diam = abs(bmax[i]-bmin[i]);
        }
    }
    
    return diam;
}

/***************************STRUCTURE INICIALIZATION**************************/
pcluster new_cluster(int start, int size, int d, int sons)
{
    int i;
    pcluster c = (pcluster) malloc(sizeof(cluster));
    c->d = d;
    c->son = (pcluster*) malloc(sons*sizeof(pcluster));
    c->start = start;
    c->size = size;
    c->sons = sons;
    c->bmin = (double*) malloc(d*sizeof(double));
    c->bmax = (double*) malloc(d*sizeof(double));

    for (i =0; i<sons; i++)
    {
        c->son[i] = (pcluster) malloc(sizeof(cluster));
    }
    for(i=0; i<d; i++)
    {
        c->bmin[i] = 0.0;
        c->bmax[i] = 0.0;
    }
    
    return c;
}

pclustertree new_clustertree(int n, int nidx)
{
    int i;
    pclustertree c = (pclustertree) malloc(sizeof(clustertree));
    c -> ndof = n;
    c -> nidx = nidx;
    c -> root = (pcluster) malloc(sizeof(cluster));
    return c;
}

pclusterfactory new_clusterfactory( int ndofp, int nidxp, int ndofb, int nidxb, int d)
{
    int i,j;
    pclusterfactory c = (pclusterfactory) malloc(sizeof(clusterfactory));
    c -> ndofp = ndofp;
    c -> nidxp = nidxp;
    c -> ndofb = ndofb;                                                                                       
    c -> nidxb = nidxb;
    c -> d = d;
    c->bmin = (double*) malloc(d*sizeof(double));
    c->bmax = (double*) malloc(d*sizeof(double));
    c->Vp = (double**) malloc(nidxp*sizeof(double*));
    c->Vb = (double**) malloc(nidxb*sizeof(double*));     
    for (i = 0; i<nidxp; i++)
    {
        c->Vp[i] = (double*) malloc(d*sizeof(double));
    }
    for (i = 0; i<nidxb; i++)                                                                                 
    {                                                                                                         
        c->Vb[i] = (double*) malloc(d*sizeof(double));                                                        
    }
    for(i=0; i<d; i++)
    {
        c->bmin[i] = 0.0;
        c->bmax[i] = 0.0;
    }
    for (i = 0; i<nidxp; i++)
    {
        for (j = 0; j<d; j++)
        {
            c->Vp[i][j] = 0.0;
        }
    }
    for (i = 0; i<nidxb; i++)                                                                                 
    {                                                                                                         
        for (j = 0; j<d; j++)                                                                                 
        {                                                                                                     
            c->Vb[i][j] = 0.0;                                                                                
        }                                                                                                     
    }  
    return c;
}

pblockcluster new_blockcluster(pcluster row, pcluster col, int block_rows, int block_cols, unsigned type)
{
    int i;
    pblockcluster b = (pblockcluster) malloc(sizeof(blockcluster));
    b->son = (pblockcluster*) malloc(block_cols*block_rows*sizeof(pblockcluster));
    b->row = (pcluster) malloc(sizeof(cluster));
    b->col = (pcluster) malloc(sizeof(cluster));
    b->row = row;
    b->col = col;
    b->block_rows = block_rows;
    b->block_cols = block_cols;
    b->type = type;
    for (i =0; i<block_cols*block_rows; i++)
    {
        b->son[i] = (pblockcluster) malloc(sizeof(blockcluster));
    }

    
    return b;
}

pfullmatrix new_fullmatrix_P0D(pcluster rows, pcluster cols, int *index, double kappa, int *order, double eps)
{
    pfullmatrix f = (pfullmatrix) malloc(sizeof(fullmatrix));
    
   
    f->rows = rows->M->PD0->ndof;                                                                             
    f->cols = cols->M->PD0->ndof;
    if(rows->M->PD0->type2 == 1)
    {
        f->e = assemble_SLD_f(rows, cols, index, kappa, order, rows->M->PD0->type, cols->M->PD0->type);
    }
    else if(rows->M->PD0->type2 == 2 || rows->M->PD0->type2 == 3)
    {
        f->e = assemble_SLD_f(rows, cols, index, kappa, order, cols->M->PD0->type, rows->M->PD0->type);
    }

    return f;
}

pfullmatrix new_fullmatrix_P1D(gsl_matrix_complex* M, pcluster rows, pcluster cols, int *index, double kappa, int *order, double eps)     
{                                                                                                             
    pfullmatrix f = (pfullmatrix) malloc(sizeof(fullmatrix));                                                 
                                                                                                                                                                      
    f->rows = rows->M->PD1->ndof;                                                                         
    f->cols = cols->M->PD1->ndof;                                                                         
    if(rows->M->PD1->type2 == 1)                                                                          
    {                                                                                                     
        f->e = assemble_SLD_f(rows, cols, index, kappa, order, rows->M->PD1->type, cols->M->PD1->type);   
    }                                                                                                     
    else if(rows->M->PD1->type2 == 2 || rows->M->PD1->type2 == 3)                                         
    {
        f->e = assemble_SLD_P1_f(M, rows->M, cols->M,order, kappa, 0); 

    }                                                                                                                                                                                                              
                                                                                                              
    return f;                                                                                                 
}
prkmatrix new_rkmatrix_P0D(pcluster rows, pcluster cols, int *index, double kappa, int *order, double eps)
{
    prkmatrix r; 
    r = ACA(NULL, rows, cols, index, kappa, order, rows->M->PD0->type, cols->M->PD0->type, rows->M->PD0->ndof, cols->M->PD0->ndof, eps);
}

prkmatrix new_rkmatrix_P1D(gsl_matrix_complex* M, pcluster rows, pcluster cols, int *index, double kappa, int *order, double eps)         
{
    prkmatrix r;
    int i,j;
    r = ACA(M, rows, cols, index, kappa, order, rows->M->PD1->type, cols->M->PD1->type, rows->M->PD1->nt, cols->M->PD1->nt, eps);

   /* r = ACA(rows, cols, index, kappa, order, rows->M->PD1->type, cols->M->PD1->type, rows->M->P1P->ndof, cols->M->P1P->ndof, eps);
    gsl_matrix_complex *P1 = average_matrix(rows->M,0);                                                  
    gsl_matrix_complex *P2 = average_matrix(cols->M,0);
    gsl_complex alpha = gsl_complex_rect(1,0);
    gsl_complex beta = gsl_complex_rect(0,0); 
    gsl_matrix_complex *res = gsl_matrix_complex_alloc(P1->size1, r->a->size1);
    gsl_matrix_complex *res2 = gsl_matrix_complex_alloc(r->b->size1, P2->size1);
    gsl_matrix_complex *res3 = gsl_matrix_complex_alloc(r->a->size1, P1->size1);
    gsl_blas_zgemm(CblasNoTrans, CblasTrans, alpha, P1, r->a, beta, res);                                  
    gsl_matrix_complex_free(P1);
    gsl_matrix_complex_transpose_memcpy(res3, res);
    gsl_matrix_complex_free(res);
    r->a = res3;                                                                                                                                                     
    gsl_blas_zgemm(CblasNoTrans, CblasTrans, alpha, r->b, P2, beta, res2); 
    r->b = res2;
    gsl_matrix_complex_free(P2);*/
    return r;
}

pfullmatrix new_fullmatrix_P1D_hlp(gsl_matrix_complex* M, pcluster rows, pcluster cols, int *index, double kappa, int *order, double eps)
{                                                                                                             
    pfullmatrix f = (pfullmatrix) malloc(sizeof(fullmatrix));                                                 
                                                                                                              
    f->rows = rows->M->PD1->ndof;                                                                             
    f->cols = cols->M->PD1->ndof;                                                                             
    f->e = assemble_HLP_f(M, rows->M, cols->M,order, kappa, 0);
    return f;                                                                                                 
}
prkmatrix new_rkmatrix_P1D_hlp(gsl_matrix_complex* M,pcluster rows, pcluster cols, int *index, double kappa, int *order, double eps)
{                                                                                                             
    prkmatrix r;                                                                                              
    int i,j;                                                                                                  
    r = ACA_hlp(M, rows, cols, index, kappa, order, rows->M->PD1->type, cols->M->PD1->type, rows->M->PD1->nt, cols->M->PD1->nt, eps);
}

psupermatrix new_supermatrix(int block_rows, int block_cols, pcluster rows, pcluster cols,  psupermatrix *s, prkmatrix rk, pfullmatrix full)
{
    psupermatrix S = (psupermatrix) malloc(sizeof(supermatrix));
    S -> rows = (pcluster) malloc(sizeof(cluster));
    S -> cols = (pcluster) malloc(sizeof(cluster));
    S -> rk = (prkmatrix) malloc(sizeof(rkmatrix));
    S -> full = (pfullmatrix) malloc(sizeof(fullmatrix));
    S -> s = (psupermatrix*) malloc(block_cols*block_rows*sizeof(psupermatrix));
    S -> rows = rows;
    S -> cols = cols;
    S -> rk = rk;
    S -> full = full;
    S -> block_rows = block_rows;
    S -> block_cols = block_cols;
    return S;
}


/*****************************STRUCTURE DELETION**************************/


void del_pcluster(pcluster c)
{
    int i;
    free(c->bmin);
    free(c->bmax);
    del_mesh(c->M);
    for(i=0; i<c->sons; i++)
    {
        if(c->son[i]!=NULL)
        {
            del_pcluster(c->son[i]);
        }
    }
    free(c);
}

void del_clustertree(pclustertree c)
{
    del_pcluster(c->root);
    free(c);
}

void del_clusterfactory(pclusterfactory c)
{
    int i;
    free(c->bmin);
    free(c->bmax);
    for (i =0; i<c->d;i++)
    {
        free(c->Vp[i]);
    }
    free(c);
}

void del_rkmatrix(prkmatrix r)
{
    //gsl_matrix_complex_free(r->a);
    //gsl_matrix_complex_free(r->b);
    free(r);
}

void del_fullmatrix(pfullmatrix f)
{
    //gsl_matrix_complex_free(f->e);
    free(f);
}

void del_supermatrix(psupermatrix s)
{
    int i,j;
    int block_cols = s->block_cols;
    int block_rows = s->block_rows;
    if(block_cols*block_rows>0)
    {
        for(j=0; j<block_cols; j++)
        {
            for(i=0; i<block_rows; i++)
            {
                del_supermatrix(s->s[i+j*block_rows]);
            }
        }
        //del_pcluster(s->rows);
        //del_pcluster(s->cols);
        //free(s);
    }
    else
    {
        if(s->rk)
        {
            del_rkmatrix(s->rk);
        }
        else
        {
            del_fullmatrix(s->full);
        }
        gsl_vector_complex_free(s->in);
        gsl_vector_complex_free(s->out);
        free(s->hfi);
        free(s->hfo);
    }
}


/*****************************************UTILLITIES********************/

pcluster split_boundingbox_P0D(pclusterfactory factory, pMeshes Mesh, int d, int* indexp, int leafsize, int startp, int sizep)
{
    int i, h, j, lp, rp, jnext;
    double pom;
    pcluster root;
    double *bmin;                                                                                             
    double *bmax;                                                                                             
    double **Vp;
    bmin = factory -> bmin;                                                                                   
    bmax = factory -> bmax;
    Vp = factory -> Vp;

    for(j=0; j<d; j++)
    {
        bmin[j] = bmax[j] = Vp[indexp[startp]][j];
        for(i = startp; i< startp+sizep; i++)
        { 
            if(Vp[indexp[i]][j] <bmin[j])
            {
                bmin[j] = Vp[indexp[i]][j];
            } 
            else if(bmax[j]< Vp[indexp[i]][j])
            {
                bmax[j] = Vp[indexp[i]][j];
            }
        }
    } 
    
    if(sizep <= leafsize)
    {   
        root = new_cluster(startp,sizep,d,0);
        for(j=0; j<d; j++)
        {
            root->bmin[j] = bmin[j];
            root->bmax[j] = bmax[j];
        }
        root -> M = sub_Mesh_P0D(Mesh, Mesh->trialbf, indexp, startp, sizep); 
    }
    else
    { 

        root = new_cluster(startp, sizep, d, 2);
        for(j=0;j<d;j++)
        {
            root->bmin[j] = bmin[j];
            root->bmax[j] = bmax[j];
        }
        root -> M = sub_Mesh_P0D(Mesh, Mesh->trialbf, indexp, startp, sizep); 
        pom = bmax[0] - bmin[0];
        jnext = 0;
        for(j = 1; j<d; j++)
        {
            if((bmax[j]-bmin[j])>pom)
            {
                pom = bmax[j] - bmin[j];
                jnext = j;
            }
        }
        pom = bmin[jnext] + 0.5*pom;   
        lp = startp;
        rp = startp+sizep-1;

        while(lp < rp) {
      
            while(lp < (sizep+startp) && Vp[indexp[lp]][jnext] <= pom)
            lp++;
            while(rp >= startp && Vp[indexp[rp]][jnext] > pom)
            rp--;
            if(lp < rp)
            {
                h = indexp[lp];
                indexp[lp] = indexp[rp];
                indexp[rp] = h;
            }
        }

        
        root->son[0] = split_boundingbox_P0D(factory, Mesh, d, indexp, leafsize, startp, lp-startp);
        root->son[1] = split_boundingbox_P0D(factory, Mesh, d, indexp, leafsize, lp, startp+sizep-lp);
  }
    return root;
}

pcluster split_boundingbox_P1D(pclusterfactory factory, pMeshes Mesh, int d, int* indexp, int *indexb, int *notriangle, int leafsize, int startp, int sizep, int startb, int sizeb, int lnt)
{                                                                                                             
    int i, h, j, lp, rp, lb, rb, jnext;                                                                               
    double pom;                                                                                               
    pcluster root;                                                                                            
    double *bmin;                                                                                             
    double *bmax;                                                                                             
    double **Vp;                                                                                              
    double **Vb;
    bmin = factory -> bmin;                                                                                   
    bmax = factory -> bmax;                                                                                   
    Vp = factory -> Vp;                                                                                       
    Vb = factory -> Vb;                                                                                                          

/*    for(j=0; j<d; j++)                                                                                        
    {                                                                                                         
        bmin[j] = bmax[j] = Vp[indexp[startp]][j];                                                            
        for(i = startp; i< startp+sizep; i++)                                                                 
        {                                                                                                     
            if(Vp[indexp[i]][j] <bmin[j])                                                                     
            {                                                                                                 
                bmin[j] = Vp[indexp[i]][j];                                                                   
            }                                                                                                 
            else if(bmax[j]< Vp[indexp[i]][j])                                                                
            {                                                                                                 
                bmax[j] = Vp[indexp[i]][j];                                                                   
            }                                                                                                 
        }                                                                                                     
    }                                                                                                         
   */  
    for(j=0; j<d; j++)                                                                                        
    {                                                                                                         
        bmin[j] = bmax[j] = Vb[indexb[startb]][j];                                                            
        for(i = startb; i< startb+sizeb; i++)                                                                 
        {                                                                                                     
            if(Vb[indexb[i]][j] <bmin[j])                                                                     
            {                                                                                                 
                bmin[j] = Vb[indexb[i]][j];                                                                   
            }                                                                                                 
            else if(bmax[j]< Vb[indexb[i]][j])                                                                
            {                                                                                                 
                bmax[j] = Vb[indexb[i]][j];                                                                   
            }                                                                                                 
        }                                                                                                     
    }
                                                                                                         
    if(sizeb <= leafsize)                                                                                     
    {                                                                                                         
        root = new_cluster(startb,sizeb,d,0);                                                                 
        for(j=0; j<d; j++)                                                                                    
        {                                                                                                     
            root->bmin[j] = bmin[j];                                                                          
            root->bmax[j] = bmax[j];                                                                          
        }                                                                                                     
        root -> M = sub_Mesh_P1D2(Mesh, Mesh->trialbf, indexp, indexb, notriangle, startp, sizep, startb, sizeb, lnt);                    
    } 
    else                                                                                                      
    {                                                                                                         
                                                                                                              
        root = new_cluster(startb, sizeb, d, 2);                                                              
        for(j=0;j<d;j++)                                                                                      
        {                                                                                                     
            root->bmin[j] = bmin[j];                                                                          
            root->bmax[j] = bmax[j];                                                                          
        }                                                                                                     
        root -> M = sub_Mesh_P1D2(Mesh, Mesh->trialbf, indexp, indexb, notriangle, startp, sizep, startb, sizeb, lnt);                    
        pom = bmax[0] - bmin[0];                                                                              
        jnext = 0;                                                                                            
        for(j = 1; j<d; j++)                                                                                  
        {                                                                                                     
            if((bmax[j]-bmin[j])>pom)                                                                         
            {                                                                                                 
                pom = bmax[j] - bmin[j];                                                                      
                jnext = j;                                                                                    
            }                                                                                                 
        }                                                                                                     
        pom = bmin[jnext] + 0.5*pom;                                                                          
        lb = startb;                                                                                          
        rb = startb+sizeb-1;                                                                                  
                                                                                                              
        while(lb < rb) {                                                                                      
                                                                                                              
            while(lb < (sizeb+startb) && Vb[indexb[lb]][jnext] <= pom)                                        
            lb++;                                                                                             
            while(rb >= startb && Vb[indexb[rb]][jnext] > pom)                                                
            rb--;                                                                                             
            if(lb < rb)                                                                                       
            {                                                                                                 
                h = indexb[lb];                                                                               
                indexb[lb] = indexb[rb];                                                                      
                indexb[rb] = h;                                                                               
            }                                                                                                 
        }                                                                                                     

        lp = startp;                                                                                          
        rp = startp+sizep-1;                                                                                  
                                                                                                              
        while(lp < rp) {                                                                                      
                                                                                                              
            while(lp < (sizep+startp) && Vp[indexp[lp]][jnext] <= pom)                                        
            lp++;                                                                                             
            while(rp >= startp && Vp[indexp[rp]][jnext] > pom)                                                
            rp--;                                                                                             
            if(lp < rp)                                                                                       
            {                                                                                                 
                h = indexp[lp];                                                                               
                indexp[lp] = indexp[rp];                                                                      
                indexp[rp] = h;                                                                               
            }                                                                                                 
        }                                                                                                              
                                                                                                              
        root->son[0] = split_boundingbox_P1D(factory, Mesh, d, indexp, indexb, NULL, leafsize, startp, lp-startp, startb, lb-startb, 0);      
        root->son[1] = split_boundingbox_P1D(factory, Mesh, d, indexp,  indexb, NULL, leafsize, lp, startp+sizep-lp, lb, startb+sizeb-lb,0);                                                                                                     
  }                                                                                                          
    return root;                                                                                              
}

pclustertree create_subclustertree(pclusterfactory factory, pMeshes Mesh, int d, const int *indexp,int ndp, int leafsize)
{
    pclustertree ct;
    
    ct = new_clustertree(ndp, factory->nidxp);
    ct->root = split_boundingbox_P0D(factory, Mesh, d, indexp, leafsize, 0, ndp);

    return ct;
}

pclustertree create_subclustertree_PD1(pclusterfactory factory, pMeshes Mesh, int d, const int *indexp,int ndp, const int *indexb,int ndb,int leafsize)
{                                                                                                             
    pclustertree ct;                                                                                          
                                                                                                              
    ct = new_clustertree(ndp, factory->nidxp);                                                                
    ct->root = split_boundingbox_P1D(factory, Mesh, d, indexp, indexb,  NULL, leafsize, 0, ndp, 0, ndb, 0);               
    
    return ct;                                                                                                
} 
pblockcluster build_blockcluster(pcluster row, pcluster col, double eta)
{
    int i,j;
    pblockcluster bc;
    double dist, diam_row,diam_col;

    dist = distance_cluster(row, col);
    diam_row = diameter_cluster(row);
    diam_col = diameter_cluster(col);

    if(diam_row < eta*dist || diam_col < eta*dist)
    {
        bc = new_blockcluster(row, col, 0, 0, HLIB_BLOCK_MINADM | HLIB_BLOCK_WEAKADM);
    }
    else if(row->sons > 0 && col->sons > 0)
    {   
        bc = new_blockcluster(row, col, row->sons, col->sons, 0);
        for(j=0; j<col->sons; j++)
        {
            for(i=0; i<row->sons; i++)
            {
                bc->son[i+j*bc->block_rows] = build_blockcluster(row->son[i], col->son[j], eta);
            }
        }
    }
    else
    {
        bc = new_blockcluster(row, col, 0, 0, 0);
    }
    return bc;
}

psupermatrix build_supermatrix_from_blockcluster_P0D(pblockcluster bc, int *index, double kappa, int *order, double eps)
{
    psupermatrix s;
    pcluster cols, rows;
    int block_rows, block_cols;
    int i, j;
    int depth = 0;
    rows = bc -> row;
    cols = bc -> col;
    block_rows = bc -> block_rows;
    block_cols = bc -> block_cols;
    prkmatrix r;
    pfullmatrix f;


        if(block_cols*block_rows>0)
        {
            s = new_supermatrix(block_rows, block_cols, rows, cols, NULL, NULL, NULL);
            for(j=0; j<block_cols; j++)
            {
                for(i=0; i<block_rows; i++)
                {
                    s->s[i+j*block_rows] = build_supermatrix_from_blockcluster_P0D(bc->son[i+j*block_rows], index, kappa, order, eps);
                    depth = MAX(depth,s->s[i+j*block_rows]->depth); 
                }
            }
            s->depth = depth + 1;
        }
        else
        {
        
            if(bc->type & HLIB_BLOCK_WEAKADM)
            {
                r = new_rkmatrix_P0D(rows, cols, index, kappa, order, eps);
                s = new_supermatrix(1, 1, rows, cols, NULL, r, NULL);
            }
            else
            {
                f = new_fullmatrix_P0D(rows, cols, index, kappa, order, eps);
                s = new_supermatrix(1, 1, rows, cols, NULL, NULL, f);
            }
            s->depth = 0;
            s->in = gsl_vector_complex_alloc(cols->size);
            s->out = gsl_vector_complex_alloc(rows->size);
            s->hfi = (int *) malloc(sizeof(int) * cols->size);
            s->hfo = (int *) malloc(sizeof(int) * rows->size);
            for(i = 0; i<cols->size; i++)                                                                         
            {                                                                                                     
                s->hfi[i] = cols->M->PD0->sbf[i]->dv;
            }
            for(i = 0; i<rows->size; i++)                                                                         
            {                                                                                                     
                s->hfo[i] = rows->M->PD0->sbf[i]->dv;
            } 
        }
        return s;                                                                                                 
}
    

psupermatrix build_supermatrix_from_blockcluster_P1D(gsl_matrix_complex* M,pblockcluster bc, int *index, double kappa, int *order, double eps, int opt)
{                                                                                                             
    psupermatrix s;                                                                                           
    pcluster cols, rows;                                                                                      
    int block_rows, block_cols;                                                                               
    int i, j;                                                                                                 
    rows = bc -> row;                                                                                         
    cols = bc -> col;                                                                                         
    block_rows = bc -> block_rows;                                                                            
    block_cols = bc -> block_cols;                                                                            
    prkmatrix r;                                                                                              
    pfullmatrix f;
    int depth = 0; 
        if(block_cols*block_rows>0)                                                                           
        {                                                                                                     
            s = new_supermatrix(block_rows, block_cols, rows, cols, NULL, NULL, NULL);                        
            for(j=0; j<block_cols; j++)                                                                       
            {                                                                                                 
                for(i=0; i<block_rows; i++)                                                                   
                {                                                                                             
                    s->s[i+j*block_rows] = build_supermatrix_from_blockcluster_P1D(M, bc->son[i+j*block_rows], index, kappa, order, eps, opt);
                    depth = MAX(depth,s->s[i+j*block_rows]->depth);                                                                                                              
                }                                                                                             
            }             
            s->depth = depth + 1;
        }                                                                                                     
        else                                                                                                  
        {                                                                                                     
            
            if(opt==0)
            {                                                                                                  
            if(bc->type & HLIB_BLOCK_WEAKADM)                                                                 
            {                                                                                                 
                r = new_rkmatrix_P1D(M, rows, cols, index, kappa, order, eps);                                   
                s = new_supermatrix(1, 1, rows, cols, NULL, r, NULL);                                         
            }                                                                                                 
            else                                                                                              
            {                                                                                                 
                f = new_fullmatrix_P1D(M, rows, cols, index, kappa, order, eps);                                 
                s = new_supermatrix(1, 1, rows, cols, NULL, NULL, f);                                         
            }
            }
            else
            {
            if(bc->type & HLIB_BLOCK_WEAKADM)                                                                 
            {   
               // f = new_fullmatrix_P1D_hlp(M, rows, cols, index, kappa, order, eps); 
                r = new_rkmatrix_P1D_hlp(M, rows, cols, index, kappa, order, eps);                                
                s = new_supermatrix(1, 1, rows, cols, NULL, r, NULL); 
                //s = new_supermatrix(1, 1, rows, cols, NULL, NULL,f); 
            }                                                                                                 
            else                                                                                              
            {                                                                                                 
                f = new_fullmatrix_P1D_hlp(M, rows, cols, index, kappa, order, eps);                              
                s = new_supermatrix(1, 1, rows, cols, NULL, NULL, f);                                         
            }
            }
            s->depth = 0; 
            s->in = gsl_vector_complex_alloc(cols->M->PD1->nt);                                                     
            s->out = gsl_vector_complex_alloc(rows->M->PD1->nt);
            s->hfi = (int *) malloc(sizeof(int) * cols->M->PD1->nt);                                                
            s->hfo = (int *) malloc(sizeof(int) * rows->M->PD1->nt);
            for(i = 0; i<cols->M->PD1->nt; i++)                                                                     
            {                                                                                                 
                s->hfi[i] = cols->M->PD1->triangles[i];   
            }                                                                                                 
            for(i = 0; i<rows->M->PD1->nt; i++)                                                                     
            {                                                                                                 
                s->hfo[i] = rows->M->PD1->triangles[i]; 
            }
        }
    
    return s;
}



psupermatrix create_hmat(pMeshes Mesh, double eta, int ls, double eps, int *order, double kappa)
{
    pclusterfactory  factory;
    pblockcluster        bct;
    psupermatrix           s;
    pclustertree  rows, cols;
    int              *indexp;
    int                    i;

    factory = new_clusterfactory(Mesh->ndp, Mesh->nvp, 0, 0, 3);
    indexp = (int *) malloc((size_t) sizeof(int) * Mesh->nvp);


    for(i=0; i<Mesh->nvp; i++)
    {
        factory->Vp[i][0] = Mesh->Vp[i][0];
        factory->Vp[i][1] = Mesh->Vp[i][1];
        factory->Vp[i][2] = Mesh->Vp[i][2];
        indexp[i] = i;
    }
        
    rows  = create_subclustertree(factory, Mesh, 3, indexp, Mesh->nvp, ls);
        
    for(i=0; i<Mesh->nvp; i++)
    {
        indexp[i] = i;
    }
    cols  = create_subclustertree(factory, Mesh, 3, indexp, Mesh->nvp, ls);
    bct = build_blockcluster(rows->root, cols->root, eta); 
    if(strcmp(bct->row->M->trialbf, "P0d") ==0 && strcmp(bct->col->M->testbf, "P0d") ==0 )                                                                  
    {
        s =  build_supermatrix_from_blockcluster_P0D(bct, indexp, kappa, order, eps);
    }
    else
    {
        gsl_matrix_complex* M;
        M = gsl_matrix_complex_calloc(Mesh->P1P->ndof, Mesh->P1P->ndof);
        s =  build_supermatrix_from_blockcluster_P1D(M, bct, indexp, kappa, order, eps,0);
    }
    s->nv = Mesh->nvp;
    s->index = (int *) malloc((size_t) sizeof(int) * Mesh->nvp);                                            
    s->index = indexp;                                                                                      
    del_clusterfactory(factory);
    return s;
    
}

psupermatrix create_hmat_hlp(pMeshes Mesh, double eta, int ls, double eps, int *order, double kappa)               
{                                                                                                             
    pclusterfactory  factory;                                                                                 
    pblockcluster        bct;                                                                                 
    psupermatrix           s;                                                                                 
    pclustertree  rows, cols;                                                                                 
    int              *indexp;                                                                                 
    int              *indexb;  
    int                    i;                                                                                 
                                                                                                              
    factory = new_clusterfactory(Mesh->ndp, Mesh->nvp, Mesh->ndb, Mesh->nvb, 3);                                                    
    indexp = (int *) malloc((size_t) sizeof(int) * Mesh->nvp);                                                
    indexb = (int *) malloc((size_t) sizeof(int) * Mesh->nvb);                                                                                                                
                                                                                                              
    for(i=0; i<Mesh->nvp; i++)                                                                                
    {                                                                                                         
        factory->Vp[i][0] = Mesh->Vp[i][0];                                                                   
        factory->Vp[i][1] = Mesh->Vp[i][1];                                                                   
        factory->Vp[i][2] = Mesh->Vp[i][2];                                                                   
        indexp[i] = i;                                                                                        
    } 
    for(i=0; i<Mesh->nvb; i++)                                                                                
    {                                                                                                         
        factory->Vb[i][0] = Mesh->Vb[i][0];                                                                   
        factory->Vb[i][1] = Mesh->Vb[i][1];                                                                   
        factory->Vb[i][2] = Mesh->Vb[i][2];                                                                   
        indexb[i] = i;                                                                                        
    }                                                                                                         
                                                                                                              
    rows  = create_subclustertree_PD1(factory, Mesh, 3, indexp, Mesh->nvp, indexb, Mesh->nvb, ls);                                   
                                                                                                              
    for(i=0; i<Mesh->nvp; i++)                                                                                
    {                                                                                                         
        indexp[i] = i;                                                                                        
    }                               
    for(i=0; i<Mesh->nvb; i++)                                                                                
    {                                                                                                         
        indexb[i] = i;                                                                                        
    }                                                                           
    cols  = create_subclustertree_PD1(factory, Mesh, 3, indexp, Mesh->nvp, indexb, Mesh->nvb, ls);                                   
    bct = build_blockcluster(rows->root, cols->root, eta);                                                    
    gsl_matrix_complex* M;                                                                                
    M = gsl_matrix_complex_calloc(Mesh->P1P->ndof, Mesh->P1P->ndof);                                      
    s =  build_supermatrix_from_blockcluster_P1D(M, bct, indexp, kappa, order, eps,1);                                                                                                                             
    gsl_matrix_complex_free(M);  
    s->nv = Mesh->nvp;                                                                                        
    s->index = (int *) malloc((size_t) sizeof(int) * Mesh->nvp);                                              
    s->index = indexp;                                                                                        
//    del_clusterfactory(factory);                                                                              
    return s;                                                                                                 
                                                                                                              
}  
