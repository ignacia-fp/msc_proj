#include "basis_functions.h"
//

/********************+************STRUCTURE DELETION********************************/


void delete_Triangle(pTriangle T)
{
    free(T->v);
    free(T->dofs);
    free(T);
}

void delete_Trianglep(pTrianglep T)
{
    int i;
    for( i = 0; i < 6; i++)
    {
        delete_Triangle(T -> T[i]);
    }
    free(T->v);
    free(T->dofs);
    free(T);
}

void delete_Quadrilateral(pQuadrilateral Q)
{
    free(Q->v);
    free(Q->dofs);
    free(Q->dofsp);
    free(Q);
}

void delete_sP0P(psP0P s)
{
    delete_Triangle(s -> T);
    free(s);
}

void delete_sP1P(psP1P s)
{
    int i;
    
    for( i = 0; i < s->len; i++)
    {
        delete_Triangle(s -> T[i]);
    }
    free( s->num );
    free(s);
    
}

void delete_sPD0(psPD0 s)
{
    int i;
    
    for( i = 0; i < s->len; i++)
    {
        delete_Quadrilateral(s -> Q[i]);
    }
    free( s->num );
    free(s);
    
}

void delete_sPD1(psPD1 s)
{
    int i;
    
    for( i = 0; i < s->len; i++)
    {
        delete_Trianglep(s -> T[i]);
    }
    free( s->num );
    free(s);
    
}

void delete_P0Pmesh(pP0Pmesh M)
{
    int i;
    for( i = 0; i < M->nt; i++)
    {
        delete_sP0P(M->sbf[i]);
    }
    free(M -> type);
    free(M);
    
}

void delete_P1Pmesh(pP1Pmesh M)
{
    int i;
    for( i = 0; i < M->ndof; i++)
    {
        delete_sP1P(M->sbf[i]);
    }
    free(M -> type);
    free(M ->globaldofs);
    free(M);
    
}

void delete_PD0mesh(pPD0mesh M)
{
    int i;
    for( i = 0; i < M->ndof; i++)
    {
        delete_sPD0(M->sbf[i]);
    }
    free(M -> type);
    free(M);
    
}

void delete_PD1mesh(pPD1mesh M)
{
    int i;
    for( i = 0; i < M->ndof; i++)
    {
        delete_sPD1(M->sbf[i]);
    }
    free(M -> triangles);
    free(M -> dofsb);
    free(M ->type);
    free(M);
    
}

void del_mesh(pMeshes M)
{
    int i;
    if(M -> P0P!=NULL)
    {
        delete_P0Pmesh(M->P0P);
    }
    else if(M -> P1P!=NULL)
    {
        delete_P1Pmesh(M->P1P);
    }
    else if(M -> PD0!=NULL)
    {
        delete_PD0mesh(M->PD0);
    }
    else if(M -> PD1!=NULL)
    {
        delete_PD1mesh(M->PD1);
    }
    for(i=0;i<M -> nvp; i++)
    {
        free(M -> Vp[i]);
    }
    for(i=0;i<M -> nvb; i++)
    {
        free(M -> Vb[i]);
    }
    free(M->Vp);
    free(M->Vb);
    free(M->testbf);
    free(M->trialbf);
    free(M);
}

/********************+************BASIS FUNCTIONS SUPPORT********************************/
pTriangle new_Triangle(double *vertices, int *dofs)
{
    int i;
    pTriangle t = (pTriangle) malloc(sizeof(Triangle));
    for( i = 0; i<3; i++)                                                                                     
    { 
        t->dofs[i] = dofs[i];
    }
    t->a = area(vertices)/2;

    for( i = 0; i<9; i++)
    {
        t->v[i] = vertices[i];
    }
    t->n = normal(NULL, t->v); 
    t->c = curl_bf(t->v, 0, t->n); 
    return t;
}

pTrianglep new_Trianglep(double *vertices, int *dofs, double **verticesb, int **dofsb, int nt)                                                           
{                                                                                                             
    int i;                                                                                                    
    pTrianglep t = (pTrianglep) malloc(sizeof(Trianglep));                                                       
    double v [9];
    for( i = 0; i<3; i++)                                                                                     
    {                                                                                                         
        t->dofs[i] = dofs[i];                                                                                 
    }                                                                                                         
    t->a = area(vertices)/2;                                                                                  
                                                                                                              
    for( i = 0; i<9; i++)                                                                                     
    {                                                                                                         
        t->v[i] = vertices[i];                                                                                
    }    

    t->T = (pTriangle*) malloc(6*sizeof(pTriangle)); 
    t->num = (int*)malloc(6*sizeof(int));
    t->n = normal(NULL, t->v);
    t->c = curl_bf(t->v, 0, t->n);
    v[0] = verticesb[dofsb[5][0]][0];                                                                     
    v[1] = verticesb[dofsb[5][0]][1];                                                                     
    v[2] = verticesb[dofsb[5][0]][2];                                                                     
    v[3] = verticesb[dofsb[5][1]][0];                                                                     
    v[4] = verticesb[dofsb[5][1]][1];                                                                     
    v[5] = verticesb[dofsb[5][1]][2];                                                                     
    v[6] = verticesb[dofsb[5][2]][0];                                                                     
    v[7] = verticesb[dofsb[5][2]][1];                                                                     
    v[8] = verticesb[dofsb[5][2]][2]; 
    t -> T[0] = new_Triangle(v, dofsb[5]);
    t-> num[0] = nt*6+5;
    v[0] = verticesb[dofsb[0][0]][0];                                                                         
    v[1] = verticesb[dofsb[0][0]][1];                                                                         
    v[2] = verticesb[dofsb[0][0]][2];                                                                         
    v[3] = verticesb[dofsb[0][1]][0];                                                                         
    v[4] = verticesb[dofsb[0][1]][1];                                                                         
    v[5] = verticesb[dofsb[0][1]][2];                                                                         
    v[6] = verticesb[dofsb[0][2]][0];                                                                         
    v[7] = verticesb[dofsb[0][2]][1];                                                                         
    v[8] = verticesb[dofsb[0][2]][2];
    t->num[1] = nt*6;
    t -> T[1] = new_Triangle(v, dofsb[0]);  
    for( i = 2; i<6; i++)                                                                                     
    {
        v[0] = verticesb[dofsb[i-1][0]][0];                                                                      
        v[1] = verticesb[dofsb[i-1][0]][1];                                                                      
        v[2] = verticesb[dofsb[i-1][0]][2];                                                                      
        v[3] = verticesb[dofsb[i-1][1]][0];                                                                      
        v[4] = verticesb[dofsb[i-1][1]][1];                                                                      
        v[5] = verticesb[dofsb[i-1][1]][2];                                                                      
        v[6] = verticesb[dofsb[i-1][2]][0];                                                                      
        v[7] = verticesb[dofsb[i-1][2]][1];                                                                      
        v[8] = verticesb[dofsb[i-1][2]][2];
        t->num[i] = nt*6+(i-1);
        t -> T[i] = new_Triangle(v, dofsb[i-1]);  
    }
    return t;                                                                                                 
}

pQuadrilateral new_Quadrilateral(double *vertices, int *dofs, int *dofsp)
{                                                                                                                                                                                                          
    pQuadrilateral t = (pQuadrilateral) malloc(sizeof(Quadrilateral));
    int i;
    for( i = 0; i<3; i++)                                                                                     
    { 
        t->dofsp[i] = dofsp[i];                                                                                 
    }                                                                                                         
    for( i = 0; i<4; i++)                                                                                     
    {                                                                                                         
        t->dofs[i] = dofs[i];                                                                               
    }                                                                                                                                                                                          
    for( i = 0; i<12; i++)                                                                                     
    {                                                                                                         
        t->v[i] = vertices[i];                                                                                
    }                                                                                                         
    return t;                                                                                                 
} 


psP0P new_sP0P(double *vertices, int *dofs)
{
     
    psP0P s = (psP0P) malloc(sizeof(sP0P));
    s -> T = new_Triangle(vertices, dofs);
        
    return s;
}

psP1P new_sP1P(double **vertices, int **dofs, int len, int *num)
{
    int i;
        
    psP1P s = (psP1P) malloc(sizeof(sP1P));
    s -> T = (pTriangle*) malloc(len*sizeof(pTriangle));
    s -> num = (int*) malloc(len*sizeof(int));
    double v[9];
    s ->dv = dofs[0][0]; 
    for(i = 0; i< len; i++)
    {
        s -> num[i] = num[i];
        v[0] = vertices[dofs[i][0]][0];
        v[1] = vertices[dofs[i][0]][1];
        v[2] = vertices[dofs[i][0]][2];
        v[3] = vertices[dofs[i][1]][0];                                                                       
        v[4] = vertices[dofs[i][1]][1];                                                                       
        v[5] = vertices[dofs[i][1]][2];
        v[6] = vertices[dofs[i][2]][0];                                                                       
        v[7] = vertices[dofs[i][2]][1];                                                                       
        v[8] = vertices[dofs[i][2]][2];
        s -> T[i] = new_Triangle(v, dofs[i]); 
    }
    
    s -> len = len;
    
    return s;
}

psPD0 new_sPD0(double **vertices, int **dofs, int **dofsp, int len, int *num)
{
    int i;
        
    psPD0 s = (psPD0) malloc(sizeof(sPD0));
    s -> Q = (pQuadrilateral*) malloc(len*sizeof(pQuadrilateral));
    s -> num = (int*) malloc(2*len*sizeof(int));
    double v[12];

    for(i = 0; i< 2*len; i++)                                                                                   
    {                                                                                                         
        s -> num[i] = num[i];
    }

    for(i = 0; i< len; i++)
    {    
        v[0] = vertices[dofs[i][0]][0];
        v[1] = vertices[dofs[i][0]][1];
        v[2] = vertices[dofs[i][0]][2];
        v[3] = vertices[dofs[i][1]][0];                                                                       
        v[4] = vertices[dofs[i][1]][1];                                                                       
        v[5] = vertices[dofs[i][1]][2];
        v[6] = vertices[dofs[i][2]][0];                                                                       
        v[7] = vertices[dofs[i][2]][1];                                                                       
        v[8] = vertices[dofs[i][2]][2];
        v[9] = vertices[dofs[i][3]][0];                                                                       
        v[10] = vertices[dofs[i][3]][1];                                                                       
        v[11] = vertices[dofs[i][3]][2];
        s -> Q[i] = new_Quadrilateral(v, dofs[i], dofsp[i]);
    }
   
    s -> len = len;
        
    return s;
}

psPD1 new_sPD1(double **vertices, double **verticesb, int **dofsp, int ***dofsb, int len, int *num)                                 
{                                                                                                             
   int i;                                                                                                    
                                                                                                              
    psPD1 s = (psPD1) malloc(sizeof(sPD1));                                                                   
    s -> T = (pTrianglep*) malloc(len*sizeof(pTrianglep));                                                      
    s -> num = (int*) malloc(len*sizeof(int));                                                                
    double v[9];                                                                                                                                                                                         
    for(i = 0; i< len; i++)                                                                                   
    {      
        s -> num[i] = num[i];  
        v[0] = vertices[dofsp[i][0]][0];                                                                       
        v[1] = vertices[dofsp[i][0]][1];                                                                       
        v[2] = vertices[dofsp[i][0]][2];                                                                       
        v[3] = vertices[dofsp[i][1]][0];                                                                       
        v[4] = vertices[dofsp[i][1]][1];                                                                       
        v[5] = vertices[dofsp[i][1]][2];                                                                       
        v[6] = vertices[dofsp[i][2]][0];                                                                       
        v[7] = vertices[dofsp[i][2]][1];                                                                       
        v[8] = vertices[dofsp[i][2]][2];                                                                       
        s -> T[i] = new_Trianglep(v, dofsp[i], verticesb, dofsb[num[i]],0);                                                                 
    }                                                                                                         
                                                                                                              
    s -> len = len;                                                                                           
                                                                                                              
    return s; 
} 

psRWG new_sRWG(double *v1, double *v2, int *dofs1, int *dofs2, int *E, int i)
{
    psRWG s = (psRWG) malloc(sizeof(sRWG));                                                                   
    s -> T1 = new_Triangle(v1, dofs1);                                                                    
    s -> T2 = new_Triangle(v2, dofs2);
    double p1[3], p2[3];

    if((dofs1[0]==dofs2[0] && dofs1[1] == dofs2[1])||(dofs1[0]==dofs2[1] && dofs1[1] == dofs2[0]))
    {
        p1[0] = v1[0];
        p1[1] = v1[1];
        p1[2] = v1[2];
        p2[0] = v1[3];                                                                                        
        p2[1] = v1[4];                                                                                        
        p2[2] = v1[5];
        s->p1[0] = v1[6];
        s->p1[1] = v1[7]; 
        s->p1[2] = v1[8]; 
        s->p2[0] = v2[6];                                                                                     
        s->p2[1] = v2[7];                                                                                     
        s->p2[2] = v2[8];
        E[2*i] = dofs1[0];
        E[2*i+1] = dofs1[1];
    }
    else if((dofs1[1]==dofs2[0] && dofs1[2] == dofs2[1])||(dofs1[1]==dofs2[1] && dofs1[2] == dofs2[0]))                                                            
    {                                                                                                         
        p1[0] = v1[3];                                                                                        
        p1[1] = v1[4];                                                                                        
        p1[2] = v1[5];                                                                                        
        p2[0] = v1[6];                                                                                        
        p2[1] = v1[7];                                                                                        
        p2[2] = v1[8]; 
        s->p1[0] = v1[0];                                                                                     
        s->p1[1] = v1[1];                                                                                     
        s->p1[2] = v1[2];                                                                                     
        s->p2[0] = v2[6];                                                                                     
        s->p2[1] = v2[7];                                                                                     
        s->p2[2] = v2[8];
        E[2*i] = dofs1[1];                                                                                    
        E[2*i+1] = dofs1[2]; 
    } 
    else if((dofs1[0]==dofs2[0] && dofs1[2] == dofs2[1])||(dofs1[0]==dofs2[1] && dofs1[2] == dofs2[0]))                                                    
    {                                                                                                         
        p1[0] = v1[0];                                                                                        
        p1[1] = v1[1];                                                                                        
        p1[2] = v1[2];                                                                                        
        p2[0] = v1[6];                                                                                        
        p2[1] = v1[7];                                                                                        
        p2[2] = v1[8]; 
        s->p1[0] = v1[3];                                                                                     
        s->p1[1] = v1[4];                                                                                     
        s->p1[2] = v1[5];                                                                                     
        s->p2[0] = v2[6];                                                                                     
        s->p2[1] = v2[7];                                                                                     
        s->p2[2] = v2[8]; 
        E[2*i] = dofs1[0];                                                                                    
        E[2*i+1] = dofs1[2]; 
    } 
    else if((dofs1[0]==dofs2[1] && dofs1[1] == dofs2[2])||(dofs1[0]==dofs2[2] && dofs1[1] == dofs2[1]))                                                          
    {                                                                                                         
        p1[0] = v1[0];                                                                                        
        p1[1] = v1[1];                                                                                        
        p1[2] = v1[2];                                                                                        
        p2[0] = v1[3];                                                                                        
        p2[1] = v1[4];                                                                                        
        p2[2] = v1[5]; 
        s->p1[0] = v1[6];                                                                                     
        s->p1[1] = v1[7];                                                                                     
        s->p1[2] = v1[8];                                                                                     
        s->p2[0] = v2[0];                                                                                     
        s->p2[1] = v2[1];                                                                                     
        s->p2[2] = v2[2]; 
        E[2*i] = dofs1[0];                                                                                    
        E[2*i+1] = dofs1[1]; 
    }                                                                                                         
    else if((dofs1[1]==dofs2[1] && dofs1[2] == dofs2[2])||(dofs1[1]==dofs2[2] && dofs1[2] == dofs2[1]))                                                       
    {                                                                                                         
        p1[0] = v1[3];                                                                                        
        p1[1] = v1[4];                                                                                        
        p1[2] = v1[5];                                                                                        
        p2[0] = v1[6];                                                                                        
        p2[1] = v1[7];                                                                                        
        p2[2] = v1[8]; 
        s->p1[0] = v1[0];                                                                                     
        s->p1[1] = v1[1];                                                                                     
        s->p1[2] = v1[2];                                                                                     
        s->p2[0] = v2[0];                                                                                     
        s->p2[1] = v2[1];                                                                                     
        s->p2[2] = v2[2]; 
        E[2*i] = dofs1[1];                                                                                    
        E[2*i+1] = dofs1[2]; 
    }                                                                                                         
    else if((dofs1[0]==dofs2[1] && dofs1[2] == dofs2[2])||(dofs1[0]==dofs2[2] && dofs1[2] == dofs2[1]))                                                      
    {                                                                                                         
        p1[0] = v1[0];                                                                                        
        p1[1] = v1[1];                                                                                        
        p1[2] = v1[2];                                                                                        
        p2[0] = v1[6];                                                                                        
        p2[1] = v1[7];                                                                                        
        p2[2] = v1[8];   
        s->p1[0] = v1[3];                                                                                     
        s->p1[1] = v1[4];                                                                                     
        s->p1[2] = v1[5];                                                                                     
        s->p2[0] = v2[0];                                                                                     
        s->p2[1] = v2[1];                                                                                     
        s->p2[2] = v2[2]; 
        E[2*i] = dofs1[0];                                                                                    
        E[2*i+1] = dofs1[2]; 
    }
    else if((dofs1[0]==dofs2[0] && dofs1[1] == dofs2[2])||(dofs1[0]==dofs2[2] && dofs1[1] == dofs2[0]))       
    {                                                                                                         
        p1[0] = v1[0];                                                                                        
        p1[1] = v1[1];                                                                                        
        p1[2] = v1[2];                                                                                        
        p2[0] = v1[3];                                                                                        
        p2[1] = v1[4];                                                                                        
        p2[2] = v1[5]; 
        s->p1[0] = v1[6];                                                                                     
        s->p1[1] = v1[7];                                                                                     
        s->p1[2] = v1[8];                                                                                     
        s->p2[0] = v2[3];                                                                                     
        s->p2[1] = v2[4];                                                                                     
        s->p2[2] = v2[5]; 
        E[2*i] = dofs1[0];                                                                                    
        E[2*i+1] = dofs1[1]; 
    }                                                                                                         
    else if((dofs1[1]==dofs2[0] && dofs1[2] == dofs2[2])||(dofs1[1]==dofs2[2] && dofs1[2] == dofs2[0]))       
    {                                                                                                         
        p1[0] = v1[3];                                                                                        
        p1[1] = v1[4];                                                                                        
        p1[2] = v1[5];                                                                                        
        p2[0] = v1[6];                                                                                        
        p2[1] = v1[7];                                                                                        
        p2[2] = v1[8];   
        s->p1[0] = v1[0];                                                                                     
        s->p1[1] = v1[1];                                                                                     
        s->p1[2] = v1[2];                                                                                     
        s->p2[0] = v2[3];                                                                                     
        s->p2[1] = v2[4];                                                                                     
        s->p2[2] = v2[5]; 
        E[2*i] = dofs1[1];                                                                                    
        E[2*i+1] = dofs1[2]; 
    }                                                                                                         
    else if((dofs1[0]==dofs2[0] && dofs1[2] == dofs2[2])||(dofs1[0]==dofs2[2] && dofs1[2] == dofs2[0]))       
    {                                                                                                         
        p1[0] = v1[0];                                                                                        
        p1[1] = v1[1];                                                                                        
        p1[2] = v1[2];                                                                                        
        p2[0] = v1[6];                                                                                        
        p2[1] = v1[7];                                                                                        
        p2[2] = v1[8]; 
        s->p1[0] = v1[3];                                                                                     
        s->p1[1] = v1[4];                                                                                     
        s->p1[2] = v1[5];                                                                                     
        s->p2[0] = v2[3];                                                                                     
        s->p2[1] = v2[4];                                                                                     
        s->p2[2] = v2[5]; 
        E[2*i] = dofs1[0];                                                                                    
        E[2*i+1] = dofs1[2]; 
    }  
    s->len = len_edge(p1, p2);
    return s;
}

psBC new_sBC(psPD0 dof1, psPD0 dof2, int ndof1, int ndof2, pTriangle *T)
{
    int i,j;
    psBC s = (psBC) malloc(sizeof(sBC));
    s->len1 = dof1->len;
    s->sbf1 = (psRWG*) malloc(dof1->len*sizeof(psRWG));
    s->sbf2 = (psRWG*) malloc(dof1->len*sizeof(psRWG)); 
    s->sbf3 = (psRWG*) malloc(dof2->len*sizeof(psRWG));                                                       
    s->sbf4 = (psRWG*) malloc(dof2->len*sizeof(psRWG)); 
    s->sbf5 = (psRWG*) malloc(2*sizeof(psRWG));
    s->coef1 = gsl_vector_complex_alloc(dof1->len);
    s->coef2 = gsl_vector_complex_alloc(dof1->len);
    s->coef3 = gsl_vector_complex_alloc(dof2->len);
    s->coef4 = gsl_vector_complex_alloc(dof2->len);
    s->coef5 = gsl_vector_complex_alloc(2);
    pQuadrilateral Qaux = (pQuadrilateral) malloc(sizeof(Quadrilateral)); 
    double p1[3], p2[3];
    int num1, num2;

    for( i = 0; i<dof1->len; i++)                                                                           
    {                                                                                                         
        for( j = 0; j<dof2->len; j++)                                                                       
        { 
            if(dof1->Q[i]->dofs[1]==dof2->Q[j]->dofs[3] && dof1->Q[i]->dofs[2]==dof2->Q[j]->dofs[2] )
            {
                Qaux = dof1->Q[i];
                num1 = dof1->num[2*i];
                num2 = dof1->num[2*i+1];
                dof1->num[2*i] = dof1->num[0];
                dof1->num[2*i+1] = dof1->num[1];
                dof1->num[0] = num1;
                dof1->num[1] = num2;
                dof1->Q[i] = dof1->Q[0];
                dof1->Q[0] = Qaux;
                Qaux = dof2->Q[j];
                num1 = dof2->num[2*j];                                                                        
                num2 = dof2->num[2*j+1];                                                                      
                dof2->num[2*j] = dof2->num[0];                                                                
                dof2->num[2*j+1] = dof2->num[1]; 
                dof2->num[0] = num1;                                                                          
                dof2->num[1] = num2; 
                dof2->Q[j] = dof2->Q[0];                                                                      
                dof2->Q[0] = Qaux; 
                i = dof1->len;
                j = dof2->len; 
            }
        }

    }

    for( i = 0; i<dof1->len-1; i++)                                                                             
    {
        for( j = i+1; j<dof1->len; j++)                                                                             
        {
            if(dof1->Q[i]->dofs[3]== dof1->Q[j]->dofs[1])
            {
                Qaux = dof1->Q[i+1];
                num1 = dof1->num[2*i+2];                                                                        
                num2 = dof1->num[2*i+3];
                dof1->num[2*i+2] = dof1->num[2*j];
                dof1->num[2*i+3] = dof1->num[2*j+1];
                dof1->num[2*j] = num1;
                dof1->num[2*j+1] = num2; 
                dof1->Q[i+1] = dof1->Q[j];
                dof1->Q[j] = Qaux;
            }
        }
    }

    for( i = 0; i<dof2->len-1; i++)                                                                           
    {                                                                                                         
        for( j = i+1; j<dof2->len; j++)                                                                       
        {                                                                                                     
            if(dof2->Q[i]->dofs[1]== dof2->Q[j]->dofs[3])                                                     
            {                                                                                                 
                Qaux = dof2->Q[i+1];  
                num1 = dof2->num[2*i+2];                                                                      
                num2 = dof2->num[2*i+3];                                                                      
                dof2->num[2*i+2] = dof2->num[2*j];                                                            
                dof2->num[2*i+3] = dof2->num[2*j+1];                                                          
                dof2->num[2*j] = num1;                                                                        
                dof2->num[2*j+1] = num2;                                                                         
                dof2->Q[i+1] = dof2->Q[j];                                                                    
                dof2->Q[j] = Qaux;                                                                            
            }                                                                                                 
        }                                                                                                     
    }
    printf("%d\n",ndof1);
    for( i = 0; i<dof1->len; i++)
    {
        s->sbf1[i] = (psRWG) malloc(sizeof(sRWG));  
        s->sbf1[i]->T1 = T[dof1->num[2*i]];
        s->sbf1[i]->T2 = T[dof1->num[2*i+1]];
        s->sbf1[i]->num[0] = dof1->num[2*i];
        s->sbf1[i]->num[1] = dof1->num[2*i+1];
        s->sbf1[i]->p1[0] = dof1->Q[i]->v[3];
        s->sbf1[i]->p1[1] = dof1->Q[i]->v[4];
        s->sbf1[i]->p1[2] = dof1->Q[i]->v[5];
        s->sbf1[i]->p2[0] = dof1->Q[i]->v[9];                                                                          
        s->sbf1[i]->p2[1] = dof1->Q[i]->v[10];                                                                          
        s->sbf1[i]->p2[2] = dof1->Q[i]->v[11]; 
        p1[0] = dof1->Q[i]->v[0];
        p1[1] = dof1->Q[i]->v[1];
        p1[2] = dof1->Q[i]->v[2];
        p2[0] = dof1->Q[i]->v[6];
        p2[1] = dof1->Q[i]->v[7];
        p2[2] = dof1->Q[i]->v[8];
        s->sbf1[i]->len = len_edge(p1, p2);  
        gsl_vector_complex_set(s->coef1,i, gsl_complex_rect((dof1->len-(2*i+1))/(2*s->sbf1[i]->len*dof1->len),0));
        printf("1: %d %d %d %d\n", dof1->Q[i]->dofs[0], dof1->Q[i]->dofs[1], dof1->Q[i]->dofs[2], dof1->Q[i]->dofs[3]);
    }
    s->sbf2[0] = (psRWG) malloc(sizeof(sRWG));  
    s->sbf2[0]->T1 = T[dof1->num[2*dof1->len-1]];
    s->sbf2[0]->T2 = T[dof1->num[0]];
    s->sbf2[0]->num[0] = dof1->num[2*dof1->len-1];                                                             
    s->sbf2[0]->num[1] = dof1->num[0];
    s->sbf2[0]->p1[0] = dof1->Q[dof1->len-1]->v[6];                                                                 
    s->sbf2[0]->p1[1] = dof1->Q[dof1->len-1]->v[7];                                                                 
    s->sbf2[0]->p1[2] = dof1->Q[dof1->len-1]->v[8];                                                                 
    s->sbf2[0]->p2[0] = dof1->Q[0]->v[6];                                                                          
    s->sbf2[0]->p2[1] = dof1->Q[0]->v[7];                                                                          
    s->sbf2[0]->p2[2] = dof1->Q[0]->v[8];
    p1[0] = dof1->Q[0]->v[0];                                                                             
    p1[1] = dof1->Q[0]->v[1];                                                                             
    p1[2] = dof1->Q[0]->v[2];                                                                             
    p2[0] = dof1->Q[0]->v[9];                                                                             
    p2[1] = dof1->Q[0]->v[10];                                                                             
    p2[2] = dof1->Q[0]->v[11];     
    s->sbf2[0]->len = len_edge(p1, p2);
    gsl_vector_complex_set(s->coef2,0, gsl_complex_rect((dof1->len-2)/(2*s->sbf2[0]->len*dof1->len),0));
    for( i = 1; i<dof1->len; i++)                                                                             
    {                                                                                                         
        s->sbf2[i] = (psRWG) malloc(sizeof(sRWG));                                                            
        s->sbf2[i]->T1 = T[dof1->num[2*i]];                                                                   
        s->sbf2[i]->T2 = T[dof1->num[2*i+1]];                                                                 
        s->sbf2[i]->num[0] = dof1->num[2*i];                                                                  
        s->sbf2[i]->num[1] = dof1->num[2*i+1];                                                                
        s->sbf2[i]->p1[0] = dof1->Q[i-1]->v[6];                                                                 
        s->sbf2[i]->p1[1] = dof1->Q[i-1]->v[7];                                                                 
        s->sbf2[i]->p1[2] = dof1->Q[i-1]->v[8];                                                                 
        s->sbf2[i]->p2[0] = dof1->Q[i]->v[6];                                                                          
        s->sbf2[i]->p2[1] = dof1->Q[i]->v[7];                                                                          
        s->sbf2[i]->p2[2] = dof1->Q[i]->v[8];                                                                
        p1[0] = dof1->Q[i-1]->v[0];                                                                             
        p1[1] = dof1->Q[i-1]->v[1];                                                                             
        p1[2] = dof1->Q[i-1]->v[2];                                                                             
        p2[0] = dof1->Q[i]->v[9];                                                                             
        p2[1] = dof1->Q[i]->v[10];                                                                             
        p2[2] = dof1->Q[i]->v[11];                                                                             
        s->sbf2[i]->len = len_edge(p1, p2);                                                                   
        gsl_vector_complex_set(s->coef2,i, gsl_complex_rect((dof1->len-(2*i+2))/(2*s->sbf2[i]->len*dof1->len),0));
   }

    printf("%d\n",ndof2);                                                                                     
    for( i = 0; i<dof2->len; i++)                                                                             
    {                                                                                                         
        s->sbf3[i] = (psRWG) malloc(sizeof(sRWG));                                                            
        s->sbf3[i]->T1 = T[dof2->num[2*i]];                                                                   
        s->sbf3[i]->T2 = T[dof2->num[2*i+1]];                                                                 
        s->sbf3[i]->num[0] = dof2->num[2*i];                                                                  
        s->sbf3[i]->num[1] = dof2->num[2*i+1];                                                                
        s->sbf3[i]->p1[0] = dof2->Q[i]->v[3];                                                                 
        s->sbf3[i]->p1[1] = dof2->Q[i]->v[4];                                                                 
        s->sbf3[i]->p1[2] = dof2->Q[i]->v[5];                                                                 
        s->sbf3[i]->p2[0] = dof2->Q[i]->v[9];                                                                 
        s->sbf3[i]->p2[1] = dof2->Q[i]->v[10];                                                                
        s->sbf3[i]->p2[2] = dof2->Q[i]->v[11];                                                                
        p1[0] = dof2->Q[i]->v[0];                                                                             
        p1[1] = dof2->Q[i]->v[1];                                                                             
        p1[2] = dof2->Q[i]->v[2];                                                                             
        p2[0] = dof2->Q[i]->v[6];                                                                             
        p2[1] = dof2->Q[i]->v[7];                                                                             
        p2[2] = dof2->Q[i]->v[8];                                                                             
        s->sbf3[i]->len = len_edge(p1, p2);                                                            
        gsl_vector_complex_set(s->coef3,i, gsl_complex_rect((dof2->len-(2*i+1))/(2*s->sbf3[i]->len*dof2->len),0));
        printf("2: %d %d %d %d\n", dof2->Q[i]->dofs[0], dof2->Q[i]->dofs[1], dof2->Q[i]->dofs[2], dof2->Q[i]->dofs[3]);
    }                                                                                                         
    s->sbf4[0] = (psRWG) malloc(sizeof(sRWG));                                                                
    s->sbf4[0]->T1 = T[dof2->num[2*dof2->len-1]];                                                             
    s->sbf4[0]->T2 = T[dof2->num[0]];                                                                         
    s->sbf4[0]->num[0] = dof2->num[2*dof2->len-1];                                                            
    s->sbf4[0]->num[1] = dof2->num[0];                                                                        
    s->sbf4[0]->p1[0] = dof2->Q[dof2->len-1]->v[6];                                                           
    s->sbf4[0]->p1[1] = dof2->Q[dof2->len-1]->v[7];                                                           
    s->sbf4[0]->p1[2] = dof2->Q[dof2->len-1]->v[8];                                                           
    s->sbf4[0]->p2[0] = dof2->Q[0]->v[6];                                                                     
    s->sbf4[0]->p2[1] = dof2->Q[0]->v[7];                                                                     
    s->sbf4[0]->p2[2] = dof2->Q[0]->v[8];                                                                     
    p1[0] = dof2->Q[0]->v[0];                                                                                 
    p1[1] = dof2->Q[0]->v[1];                                                                                 
    p1[2] = dof2->Q[0]->v[2];                                                                                 
    p2[0] = dof2->Q[0]->v[9];                                                                                 
    p2[1] = dof2->Q[0]->v[10];                                                                                
    p2[2] = dof2->Q[0]->v[11];                                                                                
    s->sbf4[0]->len = len_edge(p1, p2);                                             
    gsl_vector_complex_set(s->coef4,0, gsl_complex_rect((dof2->len-2)/(2*s->sbf4[0]->len*dof2->len),0));
    for( i = 1; i<dof2->len; i++)                                                                             
    {                                                                                                         
        s->sbf4[i] = (psRWG) malloc(sizeof(sRWG));                                                            
        s->sbf4[i]->T1 = T[dof2->num[2*i]];                                                                   
        s->sbf4[i]->T2 = T[dof2->num[2*i+1]];                                                                 
        s->sbf4[i]->num[0] = dof2->num[2*i];                                                                  
        s->sbf4[i]->num[1] = dof2->num[2*i+1];                                                                
        s->sbf4[i]->p1[0] = dof2->Q[i-1]->v[6];                                                               
        s->sbf4[i]->p1[1] = dof2->Q[i-1]->v[7];                                                               
        s->sbf4[i]->p1[2] = dof2->Q[i-1]->v[8];                                                               
        s->sbf4[i]->p2[0] = dof2->Q[i]->v[6];                                                                 
        s->sbf4[i]->p2[1] = dof2->Q[i]->v[7];                                                                 
        s->sbf4[i]->p2[2] = dof2->Q[i]->v[8];                                                                 
        p1[0] = dof2->Q[i-1]->v[0];                                                                           
        p1[1] = dof2->Q[i-1]->v[1];                                                                           
        p1[2] = dof2->Q[i-1]->v[2];
        p2[0] = dof2->Q[i]->v[6];                                                                             
        p2[1] = dof2->Q[i]->v[7];                                                                             
        p2[2] = dof2->Q[i]->v[8];                                                                             
        s->sbf4[i]->len = len_edge(p1, p2);                                                  
        gsl_vector_complex_set(s->coef4,i, gsl_complex_rect((dof2->len-(2*i+2))/(2*s->sbf4[i]->len*dof2->len),0));
    }
    s->sbf5[0] = (psRWG) malloc(sizeof(sRWG)); 
    printf("11 %d %d %d %d \n",dof1->Q[0]->dofs[1] ,T[dof1->num[0]]->dofs[0],T[dof1->num[0]]->dofs[1],T[dof1->num[0]]->dofs[2]);
    printf("22 %d %d %d %d \n",dof1->Q[0]->dofs[1] ,T[dof1->num[1]]->dofs[0],T[dof1->num[1]]->dofs[1],T[dof1->num[1]]->dofs[2]);
    if(dof1->Q[0]->dofs[1] == T[dof1->num[0]]->dofs[0] || dof1->Q[0]->dofs[1] == T[dof1->num[0]]->dofs[1] || dof1->Q[0]->dofs[1] == T[dof1->num[0]]->dofs[2] )
    {
        s->sbf5[0]->T1 = T[dof1->num[0]];
        s->sbf5[0]->num[0] = dof1->num[0];
    }
    else
    {
        s->sbf5[0]->T1 = T[dof1->num[1]];   
        s->sbf5[0]->num[0] = dof1->num[1]; 
    }
    printf("33 %d %d %d %d\n", dof2->Q[0]->dofs[3], T[dof2->num[1]]->dofs[0] , T[dof2->num[1]]->dofs[1], T[dof2->num[1]]->dofs[2] );
    printf("44 %d %d %d %d\n", dof2->Q[0]->dofs[3], T[dof2->num[0]]->dofs[0] , T[dof2->num[0]]->dofs[1], T[dof2->num[0]]->dofs[2] );
    if(dof2->Q[0]->dofs[3] == T[dof2->num[1]]->dofs[0] || dof2->Q[0]->dofs[3] == T[dof2->num[1]]->dofs[1] || dof2->Q[0]->dofs[3] == T[dof2->num[1]]->dofs[2] )
    {                                                                                                         
        s->sbf5[0]->T2 = T[dof2->num[1]];
        s->sbf5[0]->num[1] = dof2->num[1]; 
    }                                                                                                         
    else                                                                                                      
    {                                                                                                         
        s->sbf5[0]->T2 = T[dof2->num[0]]; 
        s->sbf5[0]->num[1] = dof2->num[0]; 
    } 
    s->sbf5[0]->p1[0] = dof1->Q[0]->v[0];                                                                     
    s->sbf5[0]->p1[1] = dof1->Q[0]->v[1];                                                                     
    s->sbf5[0]->p1[2] = dof1->Q[0]->v[2];                                                                     
    s->sbf5[0]->p2[0] = dof2->Q[0]->v[0];                                                                     
    s->sbf5[0]->p2[1] = dof2->Q[0]->v[1];                                                                     
    s->sbf5[0]->p2[2] = dof2->Q[0]->v[2];  
    p1[0] = dof1->Q[0]->v[3];
    p1[1] = dof1->Q[0]->v[4];
    p1[2] = dof1->Q[0]->v[5];
    p2[0] = dof1->Q[0]->v[6];
    p2[1] = dof1->Q[0]->v[7];
    p2[2] = dof1->Q[0]->v[8];
    s->sbf5[0]->len = len_edge(p1, p2);
    gsl_vector_complex_set(s->coef5,0, gsl_complex_rect(1/(2*s->sbf5[0]->len),0)); 
    s->sbf5[1] = (psRWG) malloc(sizeof(sRWG)); 
    printf("%d %d %d %d\n", dof1->Q[dof1->len-1]->dofs[1], T[dof1->num[2*dof1->len-1]]->dofs[0], T[dof1->num[2*dof1->len-1]]->dofs[1], T[dof1->num[2*dof1->len-1]]->dofs[2]);
    printf("%d %d %d %d\n", dof1->Q[dof1->len-1]->dofs[1], T[dof1->num[2*dof1->len-2]]->dofs[0], T[dof1->num[2*dof1->len-2]]->dofs[1], T[dof1->num[2*dof1->len-2]]->dofs[2]);
    if(dof1->Q[dof1->len-1]->dofs[1] == T[dof1->num[2*dof1->len-1]]->dofs[0] || dof1->Q[dof1->len-1]->dofs[1] == T[dof1->num[2*dof1->len-1]]->dofs[1] || dof1->Q[dof1->len-1]->dofs[1] == T[dof1->num[2*dof1->len-1]]->dofs[2] )
    {                                                                                                         
        s->sbf5[1]->T1 = T[dof1->num[2*dof1->len-1]];  
        s->sbf5[1]->num[0] = dof1->num[2*dof1->len-1];
    }                                                                                                         
    else                                                                                                      
    {                                                                                                         
        s->sbf5[1]->T1 = T[dof1->num[2*dof1->len-2]]; 
        s->sbf5[1]->num[0] = dof1->num[2*dof1->len-2]; 
    }
    printf("%d %d %d %d\n",dof2->Q[dof2->len-1]->dofs[1] ,T[dof2->num[2*dof2->len-2]]->dofs[0],T[dof2->num[2*dof2->len-2]]->dofs[1], T[dof2->num[2*dof2->len-2]]->dofs[2]  );
    printf("%d %d %d %d\n",dof2->Q[dof2->len-1]->dofs[1] ,T[dof2->num[2*dof2->len-1]]->dofs[0],T[dof2->num[2*dof2->len-1]]->dofs[1], T[dof2->num[2*dof2->len-1]]->dofs[2]  );
    if(dof2->Q[dof2->len-1]->dofs[1] == T[dof2->num[2*dof2->len-2]]->dofs[0] || dof2->Q[dof2->len-1]->dofs[1] == T[dof2->num[2*dof2->len-2]]->dofs[1] || dof2->Q[dof2->len-1]->dofs[1] == T[dof2->num[2*dof2->len-2]]->dofs[2] )
    {                                                                                                         
        s->sbf5[1]->T2 = T[dof2->num[2*dof2->len-2]];
        s->sbf5[1]->num[1] = dof2->num[2*dof1->len-2]; 
    }                                                                                                         
    else                                                                                                      
    {                                                                                                         
        s->sbf5[1]->T2 = T[dof2->num[2*dof2->len-1]];   
        s->sbf5[1]->num[1] = dof2->num[2*dof1->len-1]; 
    } 
    s->sbf5[1]->p1[0] = dof1->Q[dof1->len-1]->v[0];                                                                     
    s->sbf5[1]->p1[1] = dof1->Q[dof1->len-1]->v[1];                                                                     
    s->sbf5[1]->p1[2] = dof1->Q[dof1->len-1]->v[2];                                                                     
    s->sbf5[1]->p2[0] = dof2->Q[dof2->len-1]->v[0];                                                                     
    s->sbf5[1]->p2[1] = dof2->Q[dof2->len-1]->v[1];                                                                     
    s->sbf5[1]->p2[2] = dof2->Q[dof2->len-1]->v[2];                                                                     
    p1[0] = dof1->Q[dof1->len-1]->v[6];                                                                                 
    p1[1] = dof1->Q[dof1->len-1]->v[7];                                                                                 
    p1[2] = dof1->Q[dof1->len-1]->v[8];                                                                                 
    p2[0] = dof1->Q[dof1->len-1]->v[9];                                                                                 
    p2[1] = dof1->Q[dof1->len-1]->v[10];                                                                                 
    p2[2] = dof1->Q[dof1->len-1]->v[11];                                                                                 
    s->sbf5[1]->len = len_edge(p1, p2);   
    gsl_vector_complex_set(s->coef5,1, gsl_complex_rect(1/(2*s->sbf5[1]->len),0));  
    return s;
}
/********************+************BASIS FUNCTIONS MESHES********************************/

void complete_sRWG( pTriangle T1, pTriangle T2, int sz, int *E )
{
    int i;

    for( i=0; i<sz; i++ )
    {
        
        if((T1->dofs[0] ==E[2*i] && T1->dofs[1] ==E[2*i+1])||(T1->dofs[1] ==E[2*i] && T1->dofs[0] ==E[2*i+1]))
        {
            T1->edges[0] = i;
        }
        if((T1->dofs[1] ==E[2*i] && T1->dofs[2] ==E[2*i+1])||(T1->dofs[2] ==E[2*i] && T1->dofs[1] ==E[2*i+1]))
        {                                                                                                     
            T1->edges[1] = i;                                                                                 
        }
        if((T1->dofs[2] ==E[2*i] && T1->dofs[0] ==E[2*i+1])||(T1->dofs[0] ==E[2*i] && T1->dofs[2] ==E[2*i+1]))
        {                                                                                                     
            T1->edges[2] = i;                                                                                 
        }
        if((T2->dofs[0] ==E[2*i] && T2->dofs[1] ==E[2*i+1])||(T2->dofs[1] ==E[2*i] && T2->dofs[0] ==E[2*i+1]))    
        {                                                                                                     
            T2->edges[0] = i;                                                                                  
        }                                                                                                     
        if((T2->dofs[1] ==E[2*i] && T2->dofs[2] ==E[2*i+1])||(T2->dofs[2] ==E[2*i] && T2->dofs[1] ==E[2*i+1]))    
        {                                                                                                     
            T2->edges[1] = i;                                                                                  
        }                                                                                                     
        if((T2->dofs[2] ==E[2*i] && T2->dofs[0] ==E[2*i+1])||(T2->dofs[0] ==E[2*i] && T2->dofs[2] ==E[2*i+1]))    
        {                                                                                                     
            T2->edges[2] = i;                                                                                  
        }
    }
}

void complete_Triangle( pTriangle T, int sz, int **E )                                              
{                                                                                                             
    int i;                                                                                                    
                                                                                                              
    for( i=0; i<sz; i++ )                                                                                     
    {                                                                                                         
                                                                                                              
        if((T->dofs[0] ==E[i][0] && T->dofs[1] ==E[i][1])||(T->dofs[1] ==E[i][0] && T->dofs[0] ==E[i][1]))
        {                                                                                                     
            T->edges[0] = i;                                                                                 
        }                                                                                                     
        if((T->dofs[1] ==E[i][0] && T->dofs[2] ==E[i][1])||(T->dofs[2] ==E[i][0] && T->dofs[1] ==E[i][1]))
        {                                                                                                     
            T->edges[1] = i;                                                                                 
        }                                                                                                     
        if((T->dofs[2] ==E[i][0] && T->dofs[0] ==E[i][1])||(T->dofs[0] ==E[i][0] && T->dofs[2] ==E[i][1]))
        {                                                                                                     
            T->edges[2] = i;                                                                                 
        }                                                                                                     
    }                                                                                                         
}

pP0Pmesh new_P0Pmesh(double **triangles, int **dofs, int nt)
{

    int i;
    pP0Pmesh M = (pP0Pmesh) malloc(sizeof(P0Pmesh));
    M -> type = "P0";
    M -> nt = nt;
    M -> sbf = (psP0P*) malloc(nt*sizeof(psP0P));

    for (i = 0; i < nt; i++)
    {   
        M -> sbf[i] = new_sP0P(triangles[i],dofs[i]);
    }
    
    return M;
}
               
pP1Pmesh new_P1Pmesh(double **vertices, int **dofs, int ndof, int nt)
{
    int i, j, ntd;
    
    int ***A = (int***) malloc(ndof*sizeof(int**));
    int **num = (int**) malloc(ndof*sizeof(int*));
    int *len = (int*) malloc(ndof*sizeof(int));
    
    pP1Pmesh M = (pP1Pmesh) malloc(sizeof(P1Pmesh));
    
    M -> type = "P1";
    M -> ndof = ndof;
    M -> sbf = (psP1P*) malloc(ndof*sizeof(psP1P));
    M ->globaldofs = (int*) malloc(ndof*sizeof(int));   
    for (i = 0; i< ndof; i++)                                                                                
    {                                                                                                         
        M ->globaldofs[i] = i;                                                                               
    } 
        
    for (i = 0; i < ndof; i++)
    {
        A[i] = (int**) malloc(sizeof(int*));
        num[i] = (int*) malloc(sizeof(int));
        ntd = 0;
        for  (j = 0; j < nt; j++)
        {
            if ( i == dofs[j][0] )
            { 
                A[i] = (int **)realloc(A[i],sizeof(int*)*(ntd+1));
                num[i] = (int *)realloc(num[i],sizeof(int)*(ntd+1));
                A[i][ntd] = (int*) malloc(3*sizeof(int));
                num[i][ntd] = j;
                A[i][ntd][0] = dofs[j][0];
                A[i][ntd][1] = dofs[j][1];
                A[i][ntd][2] = dofs[j][2];
                ntd++;
            }
            if ( i == dofs[j][1] )                                                                            
            {                                                                                                 
                A[i] = (int **)realloc(A[i],sizeof(int*)*(ntd+1));                                            
                num[i] = (int *)realloc(num[i],sizeof(int)*(ntd+1));                                          
                A[i][ntd] = (int*) malloc(3*sizeof(int));                                                     
                num[i][ntd] = j;                                                                                                                                                  
                A[i][ntd][0] = dofs[j][1];                                                                    
                A[i][ntd][1] = dofs[j][2];                                                                    
                A[i][ntd][2] = dofs[j][0];                                                                    
                ntd++;                                                                                        
            } 
            if ( i == dofs[j][2] )                                                                            
            {                                                                                                 
                A[i] = (int **)realloc(A[i],sizeof(int*)*(ntd+1));                                            
                num[i] = (int *)realloc(num[i],sizeof(int)*(ntd+1));                                          
                A[i][ntd] = (int*) malloc(3*sizeof(int));                                                     
                num[i][ntd] = j;                                                                                                                                                   
                A[i][ntd][0] = dofs[j][2];                                                                    
                A[i][ntd][1] = dofs[j][0];                                                                    
                A[i][ntd][2] = dofs[j][1];                                                                    
                ntd++;                                                                                        
            } 

        }
        len[i] = ntd;
    }
    
    for (i = 0; i < ndof; i++)
    {
        M -> sbf[i] = new_sP1P(vertices, A[i], len[i], num[i]);
    }
    
    for (i = 0; i < ndof; i++)
    {
        free(num[i]);
        
        for(j = 0; j<len[i]; j++)
        {
            free(A[i][j]);
            
        }
    }
    
    free(len);
    
    return M;
}


pPD0mesh new_PD0mesh(double **verticesp, double **verticesb, int **dofsp, int **dofsb, int **dofsq, int ndof,  int ntb, int ntp)
{
    int i, j, ntd, nb, idx;
    int **A = (int**) malloc(ndof*sizeof(int*));
    int ***B = (int***) malloc(ndof*sizeof(int**));
    int ***C = (int***) malloc(ndof*sizeof(int**)); 
    int *len = (int*) malloc(ndof*sizeof(int));
    int *len2 = (int*) malloc(ndof*sizeof(int));
    pPD0mesh M = (pPD0mesh) malloc(sizeof(PD0mesh));    
    
    M -> type = "P0d";
    M -> ndof = ndof;
    M -> sbf = (psPD0*) malloc(ndof*sizeof(psPD0));

    idx = 0;


    for (i = 0; i < ndof; i++)
    {
        A[i] = (int*) malloc(sizeof(int));
        C[i] = (int**) malloc(sizeof(int*));
        nb = 0;
        ntd = 0;
        
        do
        {
            A[i] = (int *)realloc(A[i],sizeof(int)*(nb+2));
            A[i][nb] = 2*dofsq[idx][3];
            A[i][nb+1] = 2*dofsq[idx][3]+1;
            nb = nb+2;
            C[i] = (int**)realloc(C[i],sizeof(int*)*(ntd+1));
            C[i][ntd] = (int *)malloc(sizeof(int)*4);
            C[i][ntd][0] = dofsq[idx][0];
            C[i][ntd][1] = dofsq[idx][1];
            C[i][ntd][2] = dofsq[idx][2];
            C[i][ntd][3] = dofsq[idx][3];
            ntd++;
            idx++;
            if(idx ==ntb/2)
            {
                break;
            }
        }while(dofsq[idx][0] == dofsq[idx-1][0]); 
        len[i] = nb; 

    }

    for (i = 0; i < ndof; i++)
    {
        B[i] = (int**) malloc((len[i]/2)*sizeof(int*));
        idx = 0;
      
        for (j = 0; j < (len[i]/2); j++)
        {
            B[i][j] = (int*) malloc(4*sizeof(int)); 
            B[i][j][0] = dofsb[A[i][idx]][0];
            B[i][j][1] = dofsb[A[i][idx+1]][1];
            B[i][j][2] = dofsb[A[i][idx]][1];
            B[i][j][3] = dofsb[A[i][idx]][2];
            
            idx = idx+2;
        }
        
        len2[i] =len[i]/2;
    }
        
    
                
    for (i = 0; i < ndof; i++)
    {
        M -> sbf[i] =new_sPD0(verticesb, B[i], C[i], len2[i], A[i]);
        M -> sbf[i] -> dv = i;
    }
    
    for (i = 0; i < ndof; i++)
    {
        free(A[i]);
        
        for(j = 0; j<len2[i]; j++)
        {
            free(B[i][j]);
            free(C[i][j]);
            
        }
    }
    
    free(len);
    free(len2);
   
    return M;
}

pPD1mesh new_PD1mesh(double **vertices, double **verticesb, int **dofs, int **dofsb, int ndof, int nt, int ndofb, int ntb)
{                                                                                                             

    int i, j, k, ntd;                                                                                            
    double v[9];
    int ***A = (int***) malloc(ndof*sizeof(int**));                                                           
    int **num = (int**) malloc(ndof*sizeof(int*));                                                            
    int *len = (int*) malloc(ndof*sizeof(int));                                                               
    int ***db = (int***) malloc(nt*sizeof(int**));

    pPD1mesh M = (pPD1mesh) malloc(sizeof(PD1mesh));                                                          
                                                                                                              
    M -> type = "P1d";                                                                                         
    M -> ndof = ndof;                                                                                         
    M -> sbf = (psPD1*) malloc(ndof*sizeof(psPD1));                                                           
    M -> T = (pTrianglep*) malloc(nt*sizeof(pTrianglep));
    M -> triangles = (int*) malloc(nt*sizeof(int));
    M->lensbf = (int**) malloc(nt*sizeof(int*));
    
    for (i = 0; i < nt; i++)                                                                                
    {
        db[i] = (int**) malloc(6*sizeof(int*));
        for  (j = 0; j < 6; j++)                                                                             
        {
            db[i][j] = (int*) malloc(3*sizeof(int));
            for( k = 0; k<3; k++)
            {
                db[i][j][k] = dofsb[6*i+j][k];
            }

        }
    }
                                                                                                              
    for (i = 0; i < ndof; i++)                                                                                
    {                                                                                                         
        A[i] = (int**) malloc(sizeof(int*));                                                                  
        num[i] = (int*) malloc(sizeof(int));                                                                  
        ntd = 0;                                                                                              
        for  (j = 0; j < nt; j++)                                                                             
        {                                                                                                     
            if ( i == dofs[j][0] )                                                                            
            {                                                                                                 
                A[i] = (int **)realloc(A[i],sizeof(int*)*(ntd+1));                                            
                num[i] = (int *)realloc(num[i],sizeof(int)*(ntd+1));                                          
                A[i][ntd] = (int*) malloc(3*sizeof(int));                                                     
                num[i][ntd] = j;                                                                              
                A[i][ntd][0] = dofs[j][0];                                                                    
                A[i][ntd][1] = dofs[j][1];                                                                    
                A[i][ntd][2] = dofs[j][2];                                                                    
                ntd++;                                                                                        
            }                                                                                                 
            if ( i == dofs[j][1] )                                                                            
            {                                                                                                 
                A[i] = (int **)realloc(A[i],sizeof(int*)*(ntd+1));                                            
                num[i] = (int *)realloc(num[i],sizeof(int)*(ntd+1));                                          
                A[i][ntd] = (int*) malloc(3*sizeof(int));                                                     
                num[i][ntd] = j;                                                                              
                A[i][ntd][0] = dofs[j][1];                                                                    
                A[i][ntd][1] = dofs[j][2];                                                                    
                A[i][ntd][2] = dofs[j][0];                                                                    
                ntd++;                                                                                        
            }                                                                                                 
            if ( i == dofs[j][2] )                                                                            
            {                                                                                                 
                A[i] = (int **)realloc(A[i],sizeof(int*)*(ntd+1));                                            
                num[i] = (int *)realloc(num[i],sizeof(int)*(ntd+1));                                          
                A[i][ntd] = (int*) malloc(3*sizeof(int));                                                     
                num[i][ntd] = j;                                                                              
                A[i][ntd][0] = dofs[j][2];                                                                    
                A[i][ntd][1] = dofs[j][0];                                                                    
                A[i][ntd][2] = dofs[j][1];                                                                    
                ntd++;                                                                                        
            }                                                                                                 
                                                                                                              
        }                                                                                                     
        len[i] = ntd; 
    }                                                                                                         
                                                                                                              
    for (i = 0; i < ndof; i++)                                                                                
    {                                                                                                         
        M -> sbf[i] = new_sPD1(vertices, verticesb, A[i], db, len[i], num[i]);                                               
    }                                                                                                        
                                                                                                              
    for (i = 0; i < ndof; i++)                                                                                
    {                                                                                                         
        free(num[i]);                                                                                         
                                                                                                              
        for(j = 0; j<len[i]; j++)                                                                             
        {                                                                                                     
            free(A[i][j]);                                                                                    
                                                                                                              
        }                                                                                                     
    }                                                                                                         
    
    for(i = 0; i< nt; i++)                                                                                   
    {                                                                                                                                                                                       
        v[0] = vertices[dofs[i][0]][0];                                                                      
        v[1] = vertices[dofs[i][0]][1];                                                                      
        v[2] = vertices[dofs[i][0]][2];                                                                      
        v[3] = vertices[dofs[i][1]][0];                                                                      
        v[4] = vertices[dofs[i][1]][1];                                                                      
        v[5] = vertices[dofs[i][1]][2];                                                                      
        v[6] = vertices[dofs[i][2]][0];                                                                      
        v[7] = vertices[dofs[i][2]][1];                                                                      
        v[8] = vertices[dofs[i][2]][2];                                                                      
        M -> T[i] = new_Trianglep(v, dofs[i], verticesb, db[i], i);   
        M -> triangles[i] = i;
    }
    for(i = 0; i<nt; i++)                                                                                     
    {
        M->lensbf[i] = (int*) malloc(3*sizeof(int));                                                                                                        
        for( j = 0; j< 3; j++)                                                                                
        {                                                                                                     
            M->lensbf[i][j] = M->sbf[M->T[i]->dofs[j]]->len;                                                  
        }                                                                                                     
    }
    M->nt = nt;                                                                                                         
    free(len);                                                                                                
                                                                                                              
    return M; 

} 

pRWGmesh new_RWGmesh(double **triangles, int **edges, int **dofs,int nt, int ne)
{
    int i,j;
    pRWGmesh M = (pRWGmesh) malloc(sizeof(RWGmesh));                                                          
    M -> type = "RWG";                                                                                         
    M -> nedge = ne;                                                                                             
    M -> sbf = (psRWG*) malloc(ne*sizeof(psRWG));                                                             
    int *E = (int*) malloc(2*ne*sizeof(int)); ;
    int aux[4], aux2[4];
    int sz_idx = 0;
    int sz_sing = 0;

    for (i = 0; i < ne; i++)                                                                                  
    {                                                                                                         
        M -> sbf[i] = new_sRWG(triangles[edges[i][0]],triangles[edges[i][1]],dofs[edges[i][0]], dofs[edges[i][1]], E, i);   
        M -> sbf[i] -> edge = i;
        M -> sbf[i] -> num[0] = edges[i][0];
        M -> sbf[i] -> num[1] = edges[i][1]; 
    }                                                                                                         
                                                                                                              
    for (i = 0; i < ne; i++)
    {
        complete_sRWG(M -> sbf[i]->T1, M -> sbf[i]->T2, 2*ne, E);
    }
    
    M->sing1 = (int*) malloc(sizeof(int));                                          
    M->sing2 = (int*) malloc(sizeof(int));                                          
    M->sing3 = (int*) malloc(sizeof(int));                                          
    M->sing4 = (int*) malloc(sizeof(int)); 
    M->idx_sing = (int*) malloc(sizeof(int));
    M->idx1 = (int*) malloc(sizeof(int));


    for (i = 0; i < ne; i++)                                                                    
    {                                                                                                     
        for(j=0; j<ne; j++)                                                                     
        {                                                                                                 
                aux[0] = check_if_equal(M->sbf[i]->T1->dofs,M->sbf[j]->T1->dofs);
                aux[1] = check_if_equal(M->sbf[i]->T1->dofs,M->sbf[j]->T2->dofs);
                aux[2] = check_if_equal(M->sbf[i]->T2->dofs,M->sbf[j]->T1->dofs);
                aux[3] = check_if_equal(M->sbf[i]->T2->dofs,M->sbf[j]->T2->dofs);
                if(aux[0]==3 || aux[1] ==3 ||aux[2]==3 ||aux[3]==3)
                {   
                    M->sing1 = (int *)realloc(M->sing1,sizeof(int)*(sz_sing+1));
                    M->sing2 = (int *)realloc(M->sing2,sizeof(int)*(sz_sing+1)); 
                    M->sing3 = (int *)realloc(M->sing3,sizeof(int)*(sz_sing+1)); 
                    M->sing4 = (int *)realloc(M->sing4,sizeof(int)*(sz_sing+1)); 
                    M->idx_sing = (int *)realloc(M->idx_sing,sizeof(int)*(sz_sing+1)); 
                    M->sing1[sz_sing] = aux[0];
                    M->sing2[sz_sing] = aux[1];
                    M->sing3[sz_sing] = aux[2];
                    M->sing4[sz_sing] = aux[3];
                    M->idx_sing[sz_sing] = i*ne+j; 
                    sz_sing++;
                    
                }
                else
                {
                    M->idx1 = (int *)realloc(M->idx1,sizeof(int)*(sz_idx+1));
                    M->idx1[sz_idx] = i*ne+j; 
                    sz_idx++;
                }
        }                                                                                                 
    }
    M->sz_sing = sz_sing;
    M->sz_idx = sz_idx;

    return M;
}

pBCmesh new_BCmesh(pTriangle *T, int **edgesp2, int ne, pPD0mesh PD0)
{
    pBCmesh M = (pBCmesh) malloc(sizeof(BCmesh));
    M -> type = "BC";                                                                                        
    M -> nedge = ne;                                                                                          
    M -> sbf = (psBC*) malloc(ne*sizeof(psBC));  
    int i;

    for (i = 0; i < ne; i++)                                                                                  
    {                 

        M -> sbf[i] = new_sBC(PD0->sbf[edgesp2[i][0]], PD0->sbf[edgesp2[i][1]], edgesp2[i][0], edgesp2[i][1], T);
        M -> sbf[i] -> edge = i;                                                                                     
    } 
}

void complete_P1D(pPD1mesh PD1, pP1Pmesh P1P)
{

    int i,j,k,l,id;
    double eta1[4], eta2[4], eta3[4], chi1[4], chi2[4], chi3[4];
    int sz1[2], sz2[2], sz3[2];
    
    quad_rules_DUAL(eta1, chi1, sz1, 1);
    quad_rules_DUAL(eta2, chi2, sz2, 2);
    quad_rules_DUAL(eta3, chi3, sz3, 3);
    gsl_matrix_view c;    
    gsl_vector *c2;
    PD1 -> sbfb = (psPD1b*) malloc(PD1->nt*sizeof(psPD1b)); 
    for( i =0; i< PD1->nt; i++)
    {
        PD1->sbfb[i] = (psPD1b) malloc(sizeof(sPD1b)); 
        PD1->sbfb[i]->idx[0] = P1P->globaldofs[PD1->T[i]->dofs[0]];                                                         
        PD1->sbfb[i]->idx[1] = P1P->globaldofs[PD1->T[i]->dofs[1]];                                                         
        PD1->sbfb[i]->idx[2] = P1P->globaldofs[PD1->T[i]->dofs[2]];                                                         
        PD1->sbfb[i]->idx[3] = P1P->globaldofs[PD1->T[i]->T[0]->dofs[2]];                                                   
        PD1->sbfb[i]->idx[4] = P1P->globaldofs[PD1->T[i]->T[0]->dofs[1]];                                                   
        PD1->sbfb[i]->idx[5] = P1P->globaldofs[PD1->T[i]->T[5]->dofs[2]];                                                   
        PD1->sbfb[i]->idx[6] = P1P->globaldofs[PD1->T[i]->T[3]->dofs[2]];                                                   
        PD1->sbfb[i]->idx2[0] = PD1->T[i]->dofs[0]; 
        PD1->sbfb[i]->idx2[1] = PD1->T[i]->dofs[1];
        PD1->sbfb[i]->idx2[2] = PD1->T[i]->dofs[2];
        PD1->sbfb[i]->idx2[3] = PD1->T[i]->T[0]->dofs[2];
        PD1->sbfb[i]->idx2[4] = PD1->T[i]->T[0]->dofs[1];
        PD1->sbfb[i]->idx2[5] = PD1->T[i]->T[5]->dofs[2];
        PD1->sbfb[i]->idx2[6] = PD1->T[i]->T[3]->dofs[2];

        PD1 -> sbfb[i] -> sbf  = (psPD1b2*) malloc(7*sizeof(psPD1b2));                                                         
                                                                                   
        for(j = 0; j<3; j++)                                                                                      
        {                                                                                                         
            PD1->sbfb[i]->sbf[j] = (psPD1b2) malloc(sizeof(sPD1b2)); 
            PD1->sbfb[i]->sbf[j]->len = P1P->sbf[PD1->sbfb[i]->idx[j]]->len/2;                                                           
            PD1->sbfb[i]->sbf[j]-> Q = (pQuadrilateral*) malloc(PD1->sbfb[i]->sbf[j]->len*sizeof(pQuadrilateral));               
            id = 0;
            for(k = 0; k< P1P->sbf[PD1->sbfb[i]->idx[j]]->len; k = k+2)                                                          
            {      
                PD1->sbfb[i]->sbf[j]->Q[id] = (pQuadrilateral) malloc(sizeof(Quadrilateral));   
                PD1->sbfb[i]->sbf[j]->Q[id]->v[0] = P1P->sbf[PD1->sbfb[i]->idx[j]]->T[k]->v[0];                                         
                PD1->sbfb[i]->sbf[j]->Q[id]->v[1] = P1P->sbf[PD1->sbfb[i]->idx[j]]->T[k]->v[1];                                         
                PD1->sbfb[i]->sbf[j]->Q[id]->v[2] = P1P->sbf[PD1->sbfb[i]->idx[j]]->T[k]->v[2];                                         
                PD1->sbfb[i]->sbf[j]->Q[id]->v[3] = P1P->sbf[PD1->sbfb[i]->idx[j]]->T[k+1]->v[3];                                       
                PD1->sbfb[i]->sbf[j]->Q[id]->v[4] = P1P->sbf[PD1->sbfb[i]->idx[j]]->T[k+1]->v[4];                                       
                PD1->sbfb[i]->sbf[j]->Q[id]->v[5] = P1P->sbf[PD1->sbfb[i]->idx[j]]->T[k+1]->v[5];                                       
                PD1->sbfb[i]->sbf[j]->Q[id]->v[6] = P1P->sbf[PD1->sbfb[i]->idx[j]]->T[k]->v[3];                                         
                PD1->sbfb[i]->sbf[j]->Q[id]->v[7] = P1P->sbf[PD1->sbfb[i]->idx[j]]->T[k]->v[4];                                         
                PD1->sbfb[i]->sbf[j]->Q[id]->v[8] = P1P->sbf[PD1->sbfb[i]->idx[j]]->T[k]->v[5];                                         
                PD1->sbfb[i]->sbf[j]->Q[id]->v[9] = P1P->sbf[PD1->sbfb[i]->idx[j]]->T[k]->v[6];                                         
                PD1->sbfb[i]->sbf[j]->Q[id]->v[10] = P1P->sbf[PD1->sbfb[i]->idx[j]]->T[k]->v[7];                                        
                PD1->sbfb[i]->sbf[j]->Q[id]->v[11] = P1P->sbf[PD1->sbfb[i]->idx[j]]->T[k]->v[8];
                PD1->sbfb[i]->sbf[j]->Q[id]->n = gsl_vector_alloc(3);
                gsl_vector_memcpy(PD1->sbfb[i]->sbf[j]->Q[id]->n,normal(NULL,PD1->sbfb[i]->sbf[j]->Q[id]->v));
                PD1->sbfb[i]->sbf[j]->Q[id]->c[0] = gsl_matrix_alloc(1,3);
                //c = gsl_matrix_submatrix(PD1->sbfb[i]->sbf[j]->Q[id]->c,0,0,1,3); 
                gsl_matrix_memcpy(PD1->sbfb[i]->sbf[j]->Q[id]->c[0], curl_bfq(PD1->sbfb[i]->sbf[j]->Q[id]->v, eta1, chi1, sz1[0],0, PD1->sbfb[i]->sbf[j]->Q[id]->n));
                PD1->sbfb[i]->sbf[j]->Q[id]->c[1] = gsl_matrix_alloc(4,3);
                //c = gsl_matrix_submatrix(PD1->sbfb[i]->sbf[j]->Q[id]->c,1,0,4,3);                 
                gsl_matrix_memcpy(PD1->sbfb[i]->sbf[j]->Q[id]->c[1], curl_bfq(PD1->sbfb[i]->sbf[j]->Q[id]->v, eta2, chi2, sz2[0],0, PD1->sbfb[i]->sbf[j]->Q[id]->n));
                //c = gsl_matrix_submatrix(PD1->sbfb[i]->sbf[j]->Q[id]->c,5,0,9,3);                 
                PD1->sbfb[i]->sbf[j]->Q[id]->c[2] = gsl_matrix_alloc(9,3);
                gsl_matrix_memcpy(PD1->sbfb[i]->sbf[j]->Q[id]->c[2], curl_bfq(PD1->sbfb[i]->sbf[j]->Q[id]->v, eta3, chi3, sz3[0],0, PD1->sbfb[i]->sbf[j]->Q[id]->n)); 
                id++;                                                                                             
            }                                                                                                     
        }
        PD1->sbfb[i]->sbf[3] = (psPD1b2) malloc(sizeof(sPD1b2)); 
        PD1->sbfb[i]->sbf[3]->len = P1P->sbf[PD1->sbfb[i]->idx[j]]->len/2;   
        PD1->sbfb[i]->sbf[3]-> T = (pTriangle*) malloc(PD1->sbfb[i]->sbf[j]->len*sizeof(pTriangle));
        id = 0;
        for(k = 0; k< P1P->sbf[PD1->sbfb[i]->idx[j]]->len; k = k+2)                                       
        {
            PD1->sbfb[i]->sbf[j]->T[id] = (pTriangle) malloc(sizeof(Triangle)); 
            PD1->sbfb[i]->sbf[j]->T[id]->v[0] = PD1->T[i]->T[k]->v[6];
            PD1->sbfb[i]->sbf[j]->T[id]->v[1] = PD1->T[i]->T[k]->v[7];  
            PD1->sbfb[i]->sbf[j]->T[id]->v[2] = PD1->T[i]->T[k]->v[8];  
            PD1->sbfb[i]->sbf[j]->T[id]->v[3] = PD1->T[i]->T[k]->v[0];  
            PD1->sbfb[i]->sbf[j]->T[id]->v[4] = PD1->T[i]->T[k]->v[1];  
            PD1->sbfb[i]->sbf[j]->T[id]->v[5] = PD1->T[i]->T[k]->v[2];  
            PD1->sbfb[i]->sbf[j]->T[id]->v[6] = PD1->T[i]->T[k+1]->v[0];  
            PD1->sbfb[i]->sbf[j]->T[id]->v[7] = PD1->T[i]->T[k+1]->v[1];  
            PD1->sbfb[i]->sbf[j]->T[id]->v[8] = PD1->T[i]->T[k+1]->v[2];  
            PD1->sbfb[i]->sbf[j]->T[id]->n = normal(NULL,PD1->sbfb[i]->sbf[j]->T[id]->v);                 
            PD1->sbfb[i]->sbf[j]->T[id]->c = curl_bf(PD1->sbfb[i]->sbf[j]->T[id]->v, 0, PD1->sbfb[i]->sbf[j]->T[id]->n);
            id++;
        }
    }
}

pMeshes new_Mesh(double **Vp, double **Vb, int **Dp, int **Db, int **Dq, double **Tp, int **Ep, int **Eb, int size_of_Dp, int size_of_Db, int size_of_Vp, int size_of_Vb, int size_of_Edgep, int size_of_Edgeb, int *opt)
{
    pMeshes M = (pMeshes) malloc(sizeof(Meshes));
    int type = 0;
    int i,j;


    if(opt[0] == 0)
    {
        M->P0P = new_P0Pmesh(Tp, Dp, size_of_Dp);
        type++;
        M->P0P->type2 = type;
    }
    if(opt[0] == 1)
    {
        M->P1P = new_P1Pmesh(Vp, Dp, size_of_Vp, size_of_Dp);
        type++;
        M->P1P->type2 = type;
    }
    if(opt[0] == 2)
    {

        M->PD0 = new_PD0mesh(Vp, Vb, Dp, Db, Dq, size_of_Vp, size_of_Db, size_of_Dp);
        type++;
        M->PD0->type2 = type;
    }
    if(opt[0] == 3)
    {
        M->P1P = new_P1Pmesh(Vb, Db, size_of_Vb, size_of_Db);
        M->PD1 = new_PD1mesh(Vp, Vb, Dp, Db, size_of_Vp, size_of_Dp, size_of_Vb, size_of_Db);
        complete_P1D(M->PD1, M->P1P);
        type++;
        M->PD1->type2 = type;
    }
    if(opt[0] == 4)                                                                                           
    {                                                                                                         
        M->P0P = new_P0Pmesh(Tp, Db, size_of_Db);                                                             
        type++;                                                                                               
        M->P0P->type2 = type;
    }
    if(opt[0] == 5)                                                                                           
    {                        
        M->P1P = new_P1Pmesh(Vb, Db, size_of_Vb, size_of_Db);                                                 
        type++;                                                                                               
        M->P1P->type2 = type;                                                                                                                                                                  
    }
    if(opt[0] == 6)                                                                                           
    {                                                                                                         
        M->RWG = new_RWGmesh(Tp, Ep, Dp, size_of_Dp, size_of_Edgep);                                                 
        type++;                                                                                               
        M->RWG->type2 = type;  

    }
    if(opt[0] == 7)
    {

    }
    if(opt[1] == 0)
    {
        M->P0P = new_P0Pmesh(Tp, Dp, size_of_Dp);
        type = type+2;
        M->P0P->type2 = type;
    }
    if(opt[1] == 1)
    {
        M->P1P = new_P1Pmesh(Vp, Dp, size_of_Vp, size_of_Dp);
        type = type+2;
        M->P1P->type2 = type;
    }
    if(opt[1] == 2)
    {
        
        M->PD0 = new_PD0mesh(Vp, Vb, Dp, Db, Dq, size_of_Vp, size_of_Db, size_of_Dp);
        type = type+2;
        M->PD0->type2 = type;
    }
    if(opt[1] == 3)
    {   
        M->P1P = new_P1Pmesh(Vb, Db, size_of_Vb, size_of_Db);
        M->PD1 = new_PD1mesh(Vp, Vb, Dp, Db, size_of_Vp, size_of_Dp, size_of_Vb, size_of_Db);
        complete_P1D(M->PD1, M->P1P);
        type = type+2;
        M->PD1->type2 = type;
    }
    if(opt[1] == 4)                                                                                           
    {                                                                                                         
        M->P0P = new_P0Pmesh(Tp, Db, size_of_Db);                                                             
        type++;                                                                                               
        M->P0P->type2 = type;                                                                                 
    }
    if(opt[1] == 5)                                                                                           
    {                                                                                                         
        M->P1P = new_P1Pmesh(Vb, Db, size_of_Vb, size_of_Db);                                                 
        type++;                                                                                               
        M->P1P->type2 = type;                                                                                 
    }
    if(opt[1] == 6)                                                                                           
    {                                                                                                         
        M->RWG = new_RWGmesh(Tp, Ep, Dp, size_of_Dp, size_of_Edgep);                                             
        type++;                                                                                               
        M->RWG->type2 = type;                                                                                 
    }
    if(opt[1] == 7)                                                                                           
    {                                                                                                         



    }    

    return M;
}

pSpace new_Space(double **Vp, double **Vb, int **Dp, int **Db, int **Dq, double **Tp, double **Tb, int **Ep, int **Eb, int **Ep2, int **Eb2, int size_of_Dp, int size_of_Db, int size_of_Vp, int size_of_Vb, int size_of_Edgep, int size_of_Edgeb, int opt)
{                                                                                                             
    pSpace S = (pSpace) malloc(sizeof(Space));                                                             
    int type = 0;                                                                                             
    int i,j;                                                                                                  
                                                                                                              
    S->space = opt;

    if(opt == 0)                                                                                           
    {                                                                                                         
        S->P0P = new_P0Pmesh(Tp, Dp, size_of_Dp);                                                                                          
    }                                                                                                         
    if(opt == 1)                                                                                           
    {                                                                                                         
        S->P1P = new_P1Pmesh(Vp, Dp, size_of_Vp, size_of_Dp);                                                 
    }                                                                                                         
    if(opt == 2)                                                                                           
    {                                                                                                         
                                                                                                              
        S->PD0 = new_PD0mesh(Vp, Vb, Dp, Db, Dq, size_of_Vp, size_of_Db, size_of_Dp);                         
    }                                                                                                         
    if(opt == 3)                                                                                           
    {                                                                                                         
        S->P1P = new_P1Pmesh(Vb, Db, size_of_Vb, size_of_Db);                                                 
        S->PD1 = new_PD1mesh(Vp, Vb, Dp, Db, size_of_Vp, size_of_Dp, size_of_Vb, size_of_Db);                 
        complete_P1D(S->PD1, S->P1P);                                                                         
    }                                                                                                         
    if(opt == 4)                                                                                           
    {                                                                                                         
        S->P0P = new_P0Pmesh(Tb, Db, size_of_Db);                                                             
    }                                                                                                         
    if(opt == 5)                                                                                           
    {                                                                                                         
        S->P1P = new_P1Pmesh(Vb, Db, size_of_Vb, size_of_Db);                                                 
    }                   
    if(opt == 6)                                                                                           
    {                                                                                                         
        S->RWG = new_RWGmesh(Tp, Ep, Dp, size_of_Dp, size_of_Edgep);                                          
    } 
    if(opt == 7)
    {                 
        pTriangle *T = (pTrianglep*) malloc(size_of_Db*sizeof(pTriangle));
        for(i = 0; i<size_of_Db; i++)
        {
            T[i] = new_Triangle(Tb[i], Db[i]);
            complete_Triangle( T[i], size_of_Edgeb, Eb2 );
        }
        complete_Triangle( T, size_of_Edgeb, Eb2 );
        S->PD0 = new_PD0mesh(Vp, Vb, Dp, Db, Dq, size_of_Vp, size_of_Db, size_of_Dp);
        S->BC = new_BCmesh(T, Ep2, size_of_Edgep, S->PD0);  
    }

    return S;                                                                                                     
}
int inlist(int e, int *list, int len)
{
    int i;
    
    for(i = 0; i< len; i++)
    {
   
        if(e == list[i])
        {
            return 1;
        }
    }
    return 0;
}

int ntriangles(pPD1mesh PD1, int *notriangles, int len, int lnt)
{
    int i, j;
    int l = 0;
    int k = 0;
    PD1->triangles = (int*) malloc(sizeof(int));
    int *list;
   

    if(lnt ==0)
    {
        for(i = 0; i< len; i++)                                                                                   
        {                                                                                                         
            for( j= 0; j<PD1->sbf[i]->len; j++)                                                                   
            {                                                                                                     
                if(inlist(PD1->sbf[i]->num[j], PD1->triangles, l)==0)                                            
                {                                                                                                 
                    PD1->triangles = (int *)realloc(PD1->triangles,sizeof(int)*(l+1));                            
                    PD1->triangles[l] = PD1->sbf[i]->num[j];                                                     
                    l++;                                                                                          
                }                                                                                                 
            }                                                                                                     
        }
        return l;
    }
    else
    {
    
        list = (int*) malloc(sizeof(int)); 
        for(i = 0; i< len; i++)                                                                                   
        {                                                                                                         
            for( j = 0; j<PD1->sbf[i]->len; j++)                                                                   
            {                                                                                                     
                if(inlist(PD1->sbf[i]->num[j], list, l)==0)                                            
                {                                                                                                 
                    list = (int *)realloc(list,sizeof(int)*(l+1));                            
                    list[l] = PD1->sbf[i]->num[j];                                                     
                    l++;                                                                                          
                }                                                                                                 
            }                                                                                                     
        }
        for(i = 0; i< l; i++)
        {
           if(inlist(list[i], notriangles, lnt)==0)
           {
               PD1->triangles = (int *)realloc(PD1->triangles,sizeof(int)*(k+1));
               PD1->triangles[k] = list[i];
               k++;
            }
        }
        return k;

    }
    
 
}

pMeshes sub_Mesh_P0D(pMeshes Mesh, char *type, int *indexp, int startp, int sizep)
{
    pMeshes M = (pMeshes) malloc(sizeof(Meshes));
    int i;

    pPD0mesh PD0 = (pPD0mesh) malloc(sizeof(PD0mesh));
    M->trialbf = "P0d";
    M->testbf = "P0d"; 
    PD0 -> type = "P0d";
    PD0 -> type2 = 3;
    PD0 -> ndof = sizep;
    PD0 -> sbf = (psPD0*) malloc(sizep*sizeof(psPD0));
    for( i = 0; i<sizep; i++)
    {
        PD0 -> sbf[i] = Mesh->PD0->sbf[indexp[startp+i]];
    }
    M -> PD0 = PD0;
    return M;

}

pMeshes sub_Mesh_P1D(pMeshes Mesh, char *type, int *indexp, int *notriangle, int startp, int sizep, int lnt)      
{
    pMeshes M = (pMeshes) malloc(sizeof(Meshes));                                                             
    int i, j, k, lenbf, lendof; 
    
    M->trialbf = "P1d";                                                                                   
    M->testbf = "P1d"; 
    pPD1mesh PD1 = (pPD1mesh) malloc(sizeof(PD1mesh));
    PD1 -> type = "P1d";
    PD1 -> type2 = 3;
    PD1 -> ndof = sizep;
    PD1 -> sbf = (psPD1*) malloc(sizep*sizeof(psPD1));
    // PD1->globaldofs = (int*) malloc(Mesh->ndp*sizeof(int)); 

    for( i = 0; i<sizep; i++)
    {
        PD1 -> sbf[i] = Mesh->PD1->sbf[indexp[startp+i]];
            
    }
    int idx;
    PD1->nt = ntriangles(PD1, notriangle, sizep, lnt);
    PD1->T  = (pTrianglep*) malloc(PD1->nt*sizeof(pTrianglep));
    PD1->lensbf = (int**) malloc(PD1->nt*sizeof(int*));
        
    for( i = 0; i<PD1->nt; i++)
    {
        PD1->T[i] = Mesh->PD1->T[PD1->triangles[i]];
    }

       /* for( i = 0; i<sizep; i++)
        {

            PD1 -> globaldofs[PD1 ->sbf[i]->T[0]->dofs[0]] = i; 

        }*/
        
        /*for( i = 0; i<sizep; i++)
        {
            lenbf = 0;
            PD1 -> globaldofs[PD1 ->sbf[i]->T[0]->dofs[0]] = i;
            for( j = 0; j< PD1->sbf[i]->len; j++)
            {
                printf("NUM: %d\n", PD1->sbf[i]->num[j]);
                if(inlist(PD1->sbf[i]->num[j], PD1->triangles, PD1->nt)==1)
                {
                    lenbf++;
                }
            }
       
      
     
    
   
            PD1->lenloc[i] = lenbf;
        }*/

    PD1->dofsb = (int*) malloc(sizeof(int));
    lendof = 0;
    for(i = 0; i<PD1->nt; i++)
    {
        for(j = 0; j<6; j++)
        {
            for(k = 0; k<3; k++)
            {
                if(inlist(PD1->T[i]->T[j]->dofs[k], PD1->dofsb, lendof)==0)
                {
                    PD1->dofsb = (int *)realloc(PD1->dofsb,sizeof(int)*(lendof+1));
                    PD1->dofsb[lendof] = PD1->T[i]->T[j]->dofs[k];
                    lendof++;
                    PD1->ndofb = lendof;
                }
            }
        }
        PD1->lensbf[i] = (int*) malloc(3*sizeof(int));
        lenbf = 0;
        
        for(j = 0; j<3; j++)
        {
            /*for( k = 0; k< Mesh->PD1->sbf[PD1->T[i]->dofs[0]]->len; k++)                                                             
            {
                if(inlist(Mesh->PD1->sbf[PD1->T[i]->dofs[0]]->num[k], PD1->triangles, PD1->nt)==1)
                {
                    lenbf++;
                }
            }*/
            //PD1->lensbf[i][j] = lenbf;
          
            PD1->lensbf[i][j] = Mesh->PD1->sbf[PD1->T[i]->dofs[j]]->len;
        }
    }


    pP1Pmesh P1 = (pP1Pmesh) malloc(sizeof(P1Pmesh));
    P1 -> type = "P1";                                                                                  
    P1 -> type2 = 3;
    P1 -> sbf = (psP1P*) malloc(PD1->ndofb*sizeof(psP1P));
    P1 -> ndof = PD1->ndofb;
    P1->globaldofs = (int*) malloc(Mesh->ndb*sizeof(int));
    for( i = 0; i< PD1->ndofb; i++)
    {
        P1->sbf[i] = Mesh->P1P->sbf[PD1->dofsb[i]];
    }

    for( i = 0; i< P1->ndof; i++)                                                                       
    {                                                                                                                                                          
        P1->globaldofs[P1->sbf[i]->dv] = i;  
    }
    M -> PD1 = PD1;
    M ->P1P = P1;

    return M;
    
}


pMeshes sub_Mesh_P1D2(pMeshes Mesh, char *type, int *indexp, int *indexb, int *notriangle, int startp, int sizep, int startb, int sizeb, int lnt)  
{
    pMeshes M = (pMeshes) malloc(sizeof(Meshes));                                                             
    int i, j, k, lenbf, lendof;                                                                               
                                                                                                              
    M->trialbf = "P1d";                                                                                       
    M->testbf = "P1d";                                                                                        
    pPD1mesh PD1 = (pPD1mesh) malloc(sizeof(PD1mesh));                                                        
    PD1 -> type = "P1d";
    PD1 -> ndof = sizep;                                                                                       
    PD1 -> type2 = 3;                                                                                                                                                                               
    PD1 -> sbf = (psPD1*) malloc(sizep*sizeof(psPD1)); 
    

    for( i = 0; i<sizep; i++)                                                                                 
    {                                                                                                         
        PD1 -> sbf[i] = Mesh->PD1->sbf[indexp[startp+i]];                                                     
                                                                                                              
    }                               
    PD1->sbfb = (psPD1b*) malloc(sizeof(psPD1b)); 
    PD1->T  = (pTrianglep*) malloc(sizeof(pTrianglep));                                                                                                 
    PD1->triangles = (int*) malloc(sizeof(int)); 
    int nt = 0;
    int stop;
    for( i = 0; i<Mesh->ntp; i++)                                                                               
    {    
        stop = 0;
        for(j=0;j<sizeb & stop==0 ;j++)
        {
            if(Mesh->PD1->T[i]->dofs[0]==indexb[j+startb]|| Mesh->PD1->T[i]->dofs[1]==indexb[j+startb] || Mesh->PD1->T[i]->dofs[2]==indexb[j+startb] || Mesh->PD1->T[i]->T[0]->dofs[1]==indexb[j+startb] || Mesh->PD1->T[i]->T[5]->dofs[2]==indexb[j+startb] || Mesh->PD1->T[i]->T[3]->dofs[2]==indexb[j+startb]|| Mesh->PD1->T[i]->T[0]->dofs[2]==indexb[j+startb])
            {
             //   printf("%d %d %d %d %d %d %d\n",Mesh->PD1->T[i]->dofs[0], Mesh->PD1->T[i]->dofs[1], Mesh->PD1->T[i]->dofs[2],Mesh->PD1->T[i]->T[0]->dofs[1],Mesh->PD1->T[i]->T[5]->dofs[2], Mesh->PD1->T[i]->T[3]->dofs[2],Mesh->PD1->T[i]->T[0]->dofs[2]);
               // printf("%d %d %d %d %d %d %d\n",Mesh->PD1->sbfb[i]->idx[0],Mesh->PD1->sbfb[i]->idx[1],Mesh->PD1->sbfb[i]->idx[2],Mesh->PD1->sbfb[i]->idx[3],Mesh->PD1->sbfb[i]->idx[4],Mesh->PD1->sbfb[i]->idx[5],Mesh->PD1->sbfb[i]->idx[6]);
                PD1->T =(pTrianglep*)realloc(PD1->T,sizeof(pTrianglep)*(nt+1));
                PD1->triangles = (int *)realloc(PD1->triangles,sizeof(int)*(nt+1));  
                PD1->sbfb = (psPD1b*) realloc(PD1->sbfb,sizeof(psPD1b)*(nt+1)); 
                PD1->T[nt] = Mesh->PD1->T[i]; 
                PD1->triangles[nt] = Mesh->PD1->triangles[i]; 
                PD1->sbfb[nt] = Mesh->PD1->sbfb[i];
                nt++;
                stop =1;
            }
            
        }
    }


    PD1->nt = nt;
    PD1->lensbf = (int**) malloc(PD1->nt*sizeof(int*)); 
    PD1->dofsb = (int*) malloc(sizeof(int));                                                                  
    lendof = 0;                                                                                               
    for(i = 0; i<PD1->nt; i++)                                                                                
    {                                                                                                         
        for(j = 0; j<6; j++)                                                                                  
        {                                                                                                     
            for(k = 0; k<3; k++)                                                                              
            {                                                                                                 
                if(inlist(PD1->T[i]->T[j]->dofs[k], PD1->dofsb, lendof)==0)                                   
                {                                                                                             
                    PD1->dofsb = (int *)realloc(PD1->dofsb,sizeof(int)*(lendof+1));                           
                    PD1->dofsb[lendof] = PD1->T[i]->T[j]->dofs[k];                                            
                    lendof++;                                                                                 
                    PD1->ndofb = lendof;                                                                      
                }                                                                                             
            }                                                                                                 
        } 
    }

    for(i = 0; i<PD1->nt; i++)                                                                                
    {
        PD1->lensbf[i] = (int*) malloc(3*sizeof(int));    
        for(j = 0; j<3; j++)                                                                                  
        {
            PD1->lensbf[i][j] = Mesh->PD1->sbf[PD1->T[i]->dofs[j]]->len;                                      
        }
    }
    pP1Pmesh P1 = (pP1Pmesh) malloc(sizeof(P1Pmesh));                                                         
    P1 -> type = "P1";                                                                                        
    P1 -> type2 = 3;                                                                                          
    P1 -> sbf = (psP1P*) malloc(PD1->ndofb*sizeof(psP1P));                                                    
    P1 -> ndof = PD1->ndofb;                                                                                  
    P1-> globaldofs = (int*) malloc(Mesh->ndb*sizeof(int)); 
                                                   
    for( i = 0; i< PD1->ndofb; i++)                                                                           
    {                                                                                                         
        P1->sbf[i] = (psP1P) malloc(sizeof(sP1P));
        stop=0;
        for(j=0;j<sizeb && stop==0;j++)
        {
            if(PD1->dofsb[i]==indexb[j+startb])
            {
                P1->sbf[i]->len = Mesh->P1P->sbf[PD1->dofsb[i]]->len;
                stop=1;
            }
        }
        if(stop==0)
        {
              P1->sbf[i]->len = 0;
        }
     
        P1-> sbf[i] -> T = (pTriangle*) malloc(P1->sbf[i]->len*sizeof(pTriangle));                                                      
        P1->sbf[i] -> num = (int*) malloc(P1->sbf[i]->len*sizeof(int));                                                                                                                                     
        P1->sbf[i]->dv = Mesh->P1P->sbf[PD1->dofsb[i]]->dv;                         

        for(j = 0; j<P1->sbf[i]->len; j++)                                                                                   
        {              
            P1->sbf[i] -> num[j] = Mesh->P1P->sbf[PD1->dofsb[i]]->num[j];                                                                                 
            P1->sbf[i] -> T[j] = Mesh->P1P->sbf[PD1->dofsb[i]]->T[j];                                                               
        } 
    }                                                                                                         
                                                                                                              
    for( i = 0; i< PD1->ndofb; i++)                                                                             
    {  
       P1->globaldofs[PD1->dofsb[i]] = i;
    }

   
    M -> PD1 = PD1;                                                                                           
    M ->P1P = P1; 
    
    return M; 
}
