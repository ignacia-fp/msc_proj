#include "assembly_hmat.h"


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
                  ){

    double weights[16], weightsp[4], weightsp2[4], eta[4], chi[4], etap[4], chip[4], etap2[4], chip2[4], di[16], dj[16], vert1[9], vert2[9];
    double P0_x[3*order[0]*order[0]], P0_y[3*order[0]*order[0]], IE_x[order[0]*order[0]], IE_y[order[0]*order[0]], cres[2];
    double phii[16] = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
    double phij[16] = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
    int sz[2], szp[2], szp2[2], i, j, k, l, length_arr;
    //int order2 = 3;
    int start1 = row->start;
    int start2 = col->start;
    pMeshes Mesh1 = row->M;
    pMeshes Mesh2 = col->M;
     
    
    if( strcmp(testbf, "P0d" ) == 0 && strcmp(trialbf, "P0d") == 0 )
    {
    
        quad_rules_DUAL(eta, chi, sz, order[0]);
        weights_DUAL(weights, order[0], sz[0]);
        quad_rules_PRIMAL(etap, chip, szp, order[0]);
        weights_PRIMAL(weightsp, order[0], szp[0]);
        quad_rules_PRIMAL(etap2, chip2, szp2, order[1]);
        weights_PRIMAL(weightsp2, order[1], szp2[0]);
        length_arr = order[0]*order[0];
        gsl_matrix_complex *M = gsl_matrix_complex_alloc(Mesh1->PD0->ndof, Mesh2->PD0->ndof);
        gsl_matrix_complex_set_zero(M);
        
        for( i = 0; i < Mesh1->PD0->ndof; i++ )
        {
            for( j = 0; j < Mesh2->PD0->ndof; j++ )
            {
        
                if( index[i+start1] != index[j+start2] )
                {
                    
                    for( k = 0; k <Mesh1->PD0->sbf[i]->len; k++ )
                    {
                        
                        integration_elements_DUAL( sz[0],sz[1], eta, chi,Mesh1->PD0->sbf[i]->Q[k]->v, IE_x, P0_x );
                        ev_basis_function( length_arr, 0, phii, eta, chi, trialbf );
                        vec_prod( length_arr, phii, weights, IE_x, di );
                        
                        for(l = 0; l < Mesh2->PD0->sbf[j]->len; l++)
                        {
                            integration_elements_DUAL( sz[0], sz[1], eta, chi, Mesh2->PD0->sbf[j]->Q[l]->v, IE_y, P0_y );
                            ev_basis_function( length_arr, 0, phij, eta, chi, testbf );
                            vec_prod( length_arr, phij, weights, IE_y, dj );
                            evaluate_kernel( P0_x,P0_y, kappa, length_arr, length_arr, di, dj, cres );
                            gsl_matrix_complex_set(M, i, j, gsl_complex_add(gsl_matrix_complex_get(M,i,j),gsl_complex_rect(cres[0],cres[1])));
                        }
                    }
                }
                else
                {   
                    for( k = 0; k <Mesh1->PD0->sbf[i]->len; k++ )
                    {                                                                                                                                                                 
                        for(l = 0; l < Mesh2->PD0->sbf[j]->len; l++)
                        {      
                            if(l!=k)
                            {                                                                               
                                 vert1[0] = Mesh1->PD0->sbf[i]->Q[k]->v[0];                                   
                                 vert1[1] = Mesh1->PD0->sbf[i]->Q[k]->v[1];                                   
                                 vert1[2] = Mesh1->PD0->sbf[i]->Q[k]->v[2];                                   
                                 vert1[3] = Mesh1->PD0->sbf[i]->Q[k]->v[6];                                   
                                 vert1[4] = Mesh1->PD0->sbf[i]->Q[k]->v[7];                                   
                                 vert1[5] = Mesh1->PD0->sbf[i]->Q[k]->v[8];                                   
                                 vert1[6] = Mesh1->PD0->sbf[i]->Q[k]->v[9];                                   
                                 vert1[7] = Mesh1->PD0->sbf[i]->Q[k]->v[10];                                  
                                 vert1[8] = Mesh1->PD0->sbf[i]->Q[k]->v[11];                                  
                                 integration_elements_PRIMAL( szp[0],szp[1], etap, chip, vert1, IE_x, P0_x );
                                 ev_basis_function( length_arr, 0, phii, etap, chip, "P0" ); 
                                 vec_prod( length_arr, phii, weightsp, IE_x, di );                            
                                 vert2[0] = Mesh2->PD0->sbf[j]->Q[l]->v[0];                                   
                                 vert2[1] = Mesh2->PD0->sbf[j]->Q[l]->v[1];                                   
                                 vert2[2] = Mesh2->PD0->sbf[j]->Q[l]->v[2];                                   
                                 vert2[3] = Mesh2->PD0->sbf[j]->Q[l]->v[6];                                   
                                 vert2[4] = Mesh2->PD0->sbf[j]->Q[l]->v[7];                                   
                                 vert2[5] = Mesh2->PD0->sbf[j]->Q[l]->v[8];                                   
                                 vert2[6] = Mesh2->PD0->sbf[j]->Q[l]->v[9];                                   
                                 vert2[7] = Mesh2->PD0->sbf[j]->Q[l]->v[10];                                  
                                 vert2[8] = Mesh2->PD0->sbf[j]->Q[l]->v[11];                                  
                                 integration_elements_PRIMAL( szp[0],szp[1], etap, chip , vert2, IE_y, P0_y );
                                 ev_basis_function( length_arr, 0, phij, etap, chip, "P0" ); 
                                 vec_prod( length_arr, phij, weightsp, IE_y, dj );                            
                                 evaluate_kernel( P0_x,P0_y, kappa, length_arr, length_arr, di, dj, cres );  
                                 gsl_matrix_complex_set(M, i, j, gsl_complex_add(gsl_matrix_complex_get(M,i,j),gsl_complex_rect(cres[0],cres[1])));                                                        
                                 vert2[0] = Mesh2->PD0->sbf[j]->Q[l]->v[0];                                   
                                 vert2[1] = Mesh2->PD0->sbf[j]->Q[l]->v[1];                                   
                                 vert2[2] = Mesh2->PD0->sbf[j]->Q[l]->v[2];                                   
                                 vert2[3] = Mesh2->PD0->sbf[j]->Q[l]->v[3];                                   
                                 vert2[4] = Mesh2->PD0->sbf[j]->Q[l]->v[4];                                   
                                 vert2[5] = Mesh2->PD0->sbf[j]->Q[l]->v[5];                                   
                                 vert2[6] = Mesh2->PD0->sbf[j]->Q[l]->v[6];                                   
                                 vert2[7] = Mesh2->PD0->sbf[j]->Q[l]->v[7];                                   
                                 vert2[8] = Mesh2->PD0->sbf[j]->Q[l]->v[8];                                   
                                 integration_elements_PRIMAL( szp[0],szp[1], etap, chip , vert2, IE_y, P0_y );
                                 ev_basis_function( length_arr, 0, phij, etap, chip, "P0" ); 
                                 vec_prod( length_arr, phij, weightsp, IE_y, dj );                            
                                 evaluate_kernel( P0_x,P0_y, kappa, length_arr, length_arr, di, dj, cres );   
                                 gsl_matrix_complex_set(M, i, j, gsl_complex_add(gsl_matrix_complex_get(M,i,j),gsl_complex_rect(cres[0],cres[1])));                                     
                                 vert1[0] = Mesh1->PD0->sbf[i]->Q[k]->v[0];                                   
                                 vert1[1] = Mesh1->PD0->sbf[i]->Q[k]->v[1];                                   
                                 vert1[2] = Mesh1->PD0->sbf[i]->Q[k]->v[2];                                   
                                 vert1[3] = Mesh1->PD0->sbf[i]->Q[k]->v[3];                                   
                                 vert1[4] = Mesh1->PD0->sbf[i]->Q[k]->v[4];                                   
                                 vert1[5] = Mesh1->PD0->sbf[i]->Q[k]->v[5];                                   
                                 vert1[6] = Mesh1->PD0->sbf[i]->Q[k]->v[6];                                   
                                 vert1[7] = Mesh1->PD0->sbf[i]->Q[k]->v[7];                                   
                                 vert1[8] = Mesh1->PD0->sbf[i]->Q[k]->v[8];
                                  integration_elements_PRIMAL( szp[0],szp[1], etap, chip, vert1, IE_x, P0_x ); 
                                  ev_basis_function( length_arr, 0, phii, etap, chip, "P0" ); 
                                 vec_prod( length_arr, phii, weightsp, IE_x, di );                            
                                 vert2[0] = Mesh2->PD0->sbf[j]->Q[l]->v[0];                                   
                                 vert2[1] = Mesh2->PD0->sbf[j]->Q[l]->v[1];                                   
                                 vert2[2] = Mesh2->PD0->sbf[j]->Q[l]->v[2];                                   
                                 vert2[3] = Mesh2->PD0->sbf[j]->Q[l]->v[6];                                   
                                 vert2[4] = Mesh2->PD0->sbf[j]->Q[l]->v[7];                                   
                                 vert2[5] = Mesh2->PD0->sbf[j]->Q[l]->v[8];                                   
                                 vert2[6] = Mesh2->PD0->sbf[j]->Q[l]->v[9];                                   
                                 vert2[7] = Mesh2->PD0->sbf[j]->Q[l]->v[10];                                   
                                 vert2[8] = Mesh2->PD0->sbf[j]->Q[l]->v[11];                                  
                                 integration_elements_PRIMAL( szp[0],szp[1], etap, chip , vert2, IE_y, P0_y );
                                 ev_basis_function( length_arr, 0, phij, etap, chip, "P0" ); 
                                 vec_prod( length_arr, phij, weightsp, IE_y, dj );                            
                                 evaluate_kernel( P0_x,P0_y, kappa, length_arr, length_arr, di, dj, cres );   
                                 gsl_matrix_complex_set(M, i, j, gsl_complex_add(gsl_matrix_complex_get(M,i,j),gsl_complex_rect(cres[0],cres[1])));                                                      
                                 vert2[0] = Mesh2->PD0->sbf[j]->Q[l]->v[0];                                   
                                 vert2[1] = Mesh2->PD0->sbf[j]->Q[l]->v[1];                                   
                                 vert2[2] = Mesh2->PD0->sbf[j]->Q[l]->v[2];                                   
                                 vert2[3] = Mesh2->PD0->sbf[j]->Q[l]->v[3];                                   
                                 vert2[4] = Mesh2->PD0->sbf[j]->Q[l]->v[4];                                   
                                 vert2[5] = Mesh2->PD0->sbf[j]->Q[l]->v[5];                                   
                                 vert2[6] = Mesh2->PD0->sbf[j]->Q[l]->v[6];                                   
                                 vert2[7] = Mesh2->PD0->sbf[j]->Q[l]->v[7];                                   
                                 vert2[8] = Mesh2->PD0->sbf[j]->Q[l]->v[8];                                   
                                 integration_elements_PRIMAL( szp[0],szp[1], etap, chip , vert2, IE_y, P0_y );
                                 ev_basis_function( length_arr, 0, phij, etap, chip, "P0" ); 
                                 vec_prod( length_arr, phij, weightsp, IE_y, dj );                            
                                 evaluate_kernel( P0_x,P0_y, kappa, length_arr, length_arr, di, dj, cres );   
                                 gsl_matrix_complex_set(M, i, j, gsl_complex_add(gsl_matrix_complex_get(M,i,j),gsl_complex_rect(cres[0],cres[1])));
                             
                            }
                            else
                            {
                                 vert1[0] = Mesh1->PD0->sbf[i]->Q[k]->v[0];                                   
                                 vert1[1] = Mesh1->PD0->sbf[i]->Q[k]->v[1];                                   
                                 vert1[2] = Mesh1->PD0->sbf[i]->Q[k]->v[2];                                   
                                 vert1[3] = Mesh1->PD0->sbf[i]->Q[k]->v[6];                                   
                                 vert1[4] = Mesh1->PD0->sbf[i]->Q[k]->v[7];                                   
                                 vert1[5] = Mesh1->PD0->sbf[i]->Q[k]->v[8];                                   
                                 vert1[6] = Mesh1->PD0->sbf[i]->Q[k]->v[9];                                   
                                 vert1[7] = Mesh1->PD0->sbf[i]->Q[k]->v[10];                                  
                                 vert1[8] = Mesh1->PD0->sbf[i]->Q[k]->v[11];                                  
                                 vert2[0] = Mesh1->PD0->sbf[i]->Q[k]->v[0];                                   
                                 vert2[1] = Mesh1->PD0->sbf[i]->Q[k]->v[1];                                   
                                 vert2[2] = Mesh1->PD0->sbf[i]->Q[k]->v[2];                                   
                                 vert2[3] = Mesh1->PD0->sbf[i]->Q[k]->v[3];                                   
                                 vert2[4] = Mesh1->PD0->sbf[i]->Q[k]->v[4];                                   
                                 vert2[5] = Mesh1->PD0->sbf[i]->Q[k]->v[5];                                   
                                 vert2[6] = Mesh1->PD0->sbf[i]->Q[k]->v[6];                                   
                                 vert2[7] = Mesh1->PD0->sbf[i]->Q[k]->v[7];                                   
                                 vert2[8] = Mesh1->PD0->sbf[i]->Q[k]->v[8];                                   
                                 integration_elements_PRIMAL( szp[0],szp[1], etap, chip, vert1, IE_x, P0_x );
                                 ev_basis_function( length_arr, 0, phii, etap, chip, "P0" ); 
                                 vec_prod( length_arr, phii, weightsp, IE_x, di );                            
                                 integration_elements_PRIMAL( szp[0],szp[1], etap, chip , vert2, IE_y, P0_y );
                                 ev_basis_function( length_arr, 0, phij, etap, chip, "P0" ); 
                                 vec_prod( length_arr, phij, weightsp, IE_y, dj );                            
                                 evaluate_kernel( P0_x,P0_y, kappa, szp[0], szp[0], di, dj, cres );
                                 gsl_matrix_complex_set(M, i, j, gsl_complex_add(gsl_matrix_complex_get(M,i,j),gsl_complex_rect(cres[0],cres[1])));
                                 integration_elements_PRIMAL( szp[0],szp[1], etap, chip, vert2, IE_x, P0_x );
                                 ev_basis_function( length_arr, 0, phii, etap, chip, "P0" );
                                 vec_prod( length_arr, phii, weightsp, IE_x, di );
                                 integration_elements_PRIMAL( szp[0],szp[1], etap, chip , vert1, IE_y, P0_y );
                                 ev_basis_function( length_arr, 0, phij, etap, chip, "P0" );
                                 vec_prod( length_arr, phij, weightsp, IE_y, dj );
                                 evaluate_kernel( P0_x,P0_y, kappa, szp[0], szp[0], di, dj, cres );
                                 gsl_matrix_complex_set(M, i, j, gsl_complex_add(gsl_matrix_complex_get(M,i,j),gsl_complex_rect(cres[0],cres[1])));
                                 singular(P0_x, vert2, weightsp2, cres, IE_x, szp, order[1], 1, kappa, "P0", 1.0);
                                 gsl_matrix_complex_set(M, i, j, gsl_complex_add(gsl_matrix_complex_get(M,i,j),gsl_complex_rect(cres[0],cres[1])));
                                 singular(P0_y, vert1, weightsp2, cres, IE_y, szp, order[1], 1, kappa, "P0", 1.0);
                                 gsl_matrix_complex_set(M, i, j, gsl_complex_add(gsl_matrix_complex_get(M,i,j),gsl_complex_rect(cres[0],cres[1])));
                            }
                        }                                                                                     
                    }
                 }
            }
        }
        return M;
    }
    else if( strcmp(testbf, "P1d" ) == 0 && strcmp(trialbf, "P1d") == 0 )
    {
        double eta2[4], chi2[4], weights2[4];
        int sz2[2];
        //int order3 = 2;
        //order2 =1;
        quad_rules_PRIMAL(eta, chi, sz, order[0]);
        quad_rules_PRIMAL(eta2, chi2, sz2, order[1]);
        weights_PRIMAL(weights, order[0], sz[0]);
        weights_PRIMAL(weights2, order[1], sz2[0]);
        gsl_matrix_complex *M = gsl_matrix_complex_alloc(Mesh1->P1P->ndof, Mesh2->P1P->ndof);                 
        gsl_matrix_complex_set_zero(M);
        length_arr = sz[0];
        for( i = 0; i < Mesh1->P1P->ndof; i++ )
        {
            for( j = 0; j < Mesh2->P1P->ndof; j++ )
            {
                
                for( k = 0; k <Mesh1->P1P->sbf[i]->len; k++ )
                {
                    integration_elements_PRIMAL( sz[0],sz[1], eta, chi, Mesh1->P1P->sbf[i]->T[k]->v, IE_x, P0_x );
    
                    ev_basis_function( length_arr, 0, phii, eta, chi, trialbf );
                    vec_prod( length_arr, phii, weights, IE_x, di );
                        
                    for(l = 0; l < Mesh2->P1P->sbf[j]->len; l++)
                    {
                        if( Mesh1->P1P->sbf[i]->num[k] != Mesh2->P1P->sbf[j]->num[l] )
                        {
                            integration_elements_PRIMAL( sz[0],sz[1], eta, chi , Mesh2->P1P->sbf[j]->T[l]->v, IE_y, P0_y );
                            ev_basis_function( length_arr, 0, phij, eta, chi, testbf );
                            vec_prod( length_arr, phij, weights, IE_y, dj );
                            evaluate_kernel( P0_x,P0_y, kappa, length_arr, length_arr, di, dj, cres );
                            gsl_matrix_complex_set(M, i, j, gsl_complex_add(gsl_matrix_complex_get(M,i,j),gsl_complex_rect(cres[0],cres[1])));
                        }
                        else
                        {
                            integration_elements_PRIMAL( sz2[0],sz2[1], eta2, chi2, Mesh2->P1P->sbf[j]->T[l]->v, IE_x, P0_x );
                            singular(P0_x, Mesh2->P1P->sbf[j]->T[l]->v, weights2, cres, IE_x,sz2, order[2], 0, kappa, "P1", 1.0);
                            gsl_matrix_complex_set(M, i, j, gsl_complex_add(gsl_matrix_complex_get(M,i,j),gsl_complex_rect(cres[0],cres[1])));
                        }
                    }
                }
            }
        }
        
        gsl_complex alpha = gsl_complex_rect(1,0);
        gsl_complex beta = gsl_complex_rect(0,0);
        gsl_matrix_complex *P1 = average_matrix(Mesh1,0);
        gsl_matrix_complex *P2 = average_matrix(Mesh2,0);
        gsl_matrix_complex *res = gsl_matrix_complex_alloc(P1->size1, M->size2);
        gsl_matrix_complex *res2 = gsl_matrix_complex_alloc(P1->size1, P2->size1);
        gsl_blas_zgemm(CblasNoTrans, CblasNoTrans, alpha, P1, M, beta, res);
        gsl_matrix_complex_free(P1);
        gsl_matrix_complex_free(M); 
        gsl_blas_zgemm(CblasNoTrans, CblasTrans, alpha, res, P2, beta, res2);
        gsl_matrix_complex_free(P2);                                                                          
        gsl_matrix_complex_free(res);
        return res2;
    }

}

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
                  ) {
    int i, j, k, l;                                                                                       
        int ndofb;                                                                                            
                                                                                                              
                                                                                                              
                                                                                                              
        gsl_complex aux1, aux2, aux3;                                                                         
        aux1 = gsl_complex_rect(0.0,0.0);                                                                     
        aux2 = gsl_complex_rect(0.0,0.0);                                                                     
        aux3 = gsl_complex_rect(0.0,0.0);                                                                     
        double sing[2];                                                                                       
        int idx1[7], idx2[7];                                                                                 
        

        if( strcmp(mode, "row" ) == 0)                                                                        
        { 
            i = urow;
            gsl_vector_complex_view auxm = gsl_matrix_complex_row(m, it);                                     
            gsl_vector_complex_set_zero(&auxm.vector);
            gsl_vector_complex *f1 = gsl_vector_complex_alloc(7);                                                 
            gsl_vector_complex *f2 = gsl_vector_complex_alloc(7);                                                 
            gsl_vector_complex *f3 = gsl_vector_complex_alloc(7);                                                 
            gsl_vector_complex *f4 = gsl_vector_complex_alloc(7);                                                 
            gsl_vector_complex_set(f1, 3, gsl_complex_rect(1.0,0));                                               
            gsl_vector_complex_set(f1, 4, gsl_complex_rect(0.5,0));                                               
            gsl_vector_complex_set(f1, 5, gsl_complex_rect(0.5,0));                                               
            gsl_vector_complex_set(f1, 6, gsl_complex_rect(0.5,0));                                               
            gsl_vector_complex_set(f2, 3, gsl_complex_rect(1.0,0));                                               
            gsl_vector_complex_set(f2, 4, gsl_complex_rect(0.5,0));                                               
            gsl_vector_complex_set(f2, 5, gsl_complex_rect(0.5,0));                                               
            gsl_vector_complex_set(f2, 6, gsl_complex_rect(0.5,0));                                               
                                                                                                                                                                                                       
            idx1[0] = Mesh1->P1P->globaldofs[Mesh1->PD1->T[i]->dofs[0]];                                      
            idx1[1] = Mesh1->P1P->globaldofs[Mesh1->PD1->T[i]->dofs[1]];                                      
            idx1[2] = Mesh1->P1P->globaldofs[Mesh1->PD1->T[i]->dofs[2]];                                      
            idx1[3] = Mesh1->P1P->globaldofs[Mesh1->PD1->T[i]->T[0]->dofs[2]];                                
            idx1[4] = Mesh1->P1P->globaldofs[Mesh1->PD1->T[i]->T[0]->dofs[1]];                                
            idx1[5] = Mesh1->P1P->globaldofs[Mesh1->PD1->T[i]->T[1]->dofs[2]];                                                                                                                                            idx1[6] = Mesh1->P1P->globaldofs[Mesh1->PD1->T[i]->T[3]->dofs[2]];
            gsl_vector_complex_set(f1, 0, gsl_complex_rect(1.0/(Mesh1->PD1->lensbf[i][0]),0));                
            gsl_vector_complex_set(f1, 1, gsl_complex_rect(1.0/(Mesh1->PD1->lensbf[i][1]),0));                
            gsl_vector_complex_set(f1, 2, gsl_complex_rect(1.0/(Mesh1->PD1->lensbf[i][2]),0));                
                                                                                                              
            for( j = 0; j < Mesh2->PD1->nt; j++ )                                                             
            {                                                                                                 
                idx2[0] = Mesh2->P1P->globaldofs[Mesh2->PD1->T[j]->dofs[0]];                                  
                idx2[1] = Mesh2->P1P->globaldofs[Mesh2->PD1->T[j]->dofs[1]];                                  
                idx2[2] = Mesh2->P1P->globaldofs[Mesh2->PD1->T[j]->dofs[2]];                                  
                idx2[3] = Mesh2->P1P->globaldofs[Mesh2->PD1->T[j]->T[0]->dofs[2]];                            
                idx2[4] = Mesh2->P1P->globaldofs[Mesh2->PD1->T[j]->T[0]->dofs[1]];                            
                idx2[5] = Mesh2->P1P->globaldofs[Mesh2->PD1->T[j]->T[1]->dofs[2]];                            
                idx2[6] = Mesh2->P1P->globaldofs[Mesh2->PD1->T[j]->T[3]->dofs[2]];                            
                gsl_vector_complex_set(f2, 0, gsl_complex_rect(1.0/(Mesh2->PD1->lensbf[j][0]),0));            
                gsl_vector_complex_set(f2, 1, gsl_complex_rect(1.0/(Mesh2->PD1->lensbf[j][1]),0));            
                gsl_vector_complex_set(f2, 2, gsl_complex_rect(1.0/(Mesh2->PD1->lensbf[j][2]),0));   
                
                for(k=0;k<7;k++)                                                                              
                {                                                                                             
                    for(l=0;l<7;l++)                                                                          
                    {                                                                                         
                        if(GSL_REAL(gsl_matrix_complex_get(Mat,idx1[k], idx2[l]))==0 & GSL_IMAG(gsl_matrix_complex_get(Mat,idx1[k], idx2[l]))==0)
                        {                                                                                     
                                sing[0] = 0;                                                                  
                                sing[1] = 0;                                                                  
                                aux3 = assemble_P1_HP_entries2(order, kappa, i, j, idx1[k], idx2[l],sing, Mesh1,Mesh2, k, l);
                                aux3 = gsl_complex_add(aux3,gsl_complex_rect(sing[0],sing[1]));               
                                gsl_vector_complex_set(f3, l, aux3);                                          
                                gsl_matrix_complex_set(Mat,idx1[k], idx2[l],aux3);                            
                         }                                                                                    
                        else                                                                                  
                        {                                                                                     
                           gsl_vector_complex_set(f3, l, gsl_matrix_complex_get(Mat,idx1[k], idx2[l]));       
                        }                                                                                     
                    }                                                                                         
                    gsl_blas_zdotu(f3,f2,&aux1);                                                              
                    gsl_vector_complex_set(f4, k, aux1);                                                      
                }                                                                                             
                                                                                                              
                gsl_blas_zdotu(f4,f1,&aux2);
                gsl_vector_complex_set(&auxm.vector, j,aux2);                                                                   
                                                           
            }                                                                                                 
                                                                                                             
        }
        else
        {
            j = ucol;                                                                                         
            gsl_vector_complex_view auxm = gsl_matrix_complex_row(m, it);                                     
            gsl_vector_complex_set_zero(&auxm.vector);                                                        
            gsl_vector_complex *f1 = gsl_vector_complex_alloc(7);                                             
            gsl_vector_complex *f2 = gsl_vector_complex_alloc(7);                                             
            gsl_vector_complex *f3 = gsl_vector_complex_alloc(7);                                             
            gsl_vector_complex *f4 = gsl_vector_complex_alloc(7);                                             
            gsl_vector_complex_set(f1, 3, gsl_complex_rect(1.0,0));                                           
            gsl_vector_complex_set(f1, 4, gsl_complex_rect(0.5,0));                                           
            gsl_vector_complex_set(f1, 5, gsl_complex_rect(0.5,0));                                           
            gsl_vector_complex_set(f1, 6, gsl_complex_rect(0.5,0));                                           
            gsl_vector_complex_set(f2, 3, gsl_complex_rect(1.0,0));                                           
            gsl_vector_complex_set(f2, 4, gsl_complex_rect(0.5,0));                                           
            gsl_vector_complex_set(f2, 5, gsl_complex_rect(0.5,0));                                           
            gsl_vector_complex_set(f2, 6, gsl_complex_rect(0.5,0));                                           
                                                                                                              
            idx2[0] = Mesh2->P1P->globaldofs[Mesh2->PD1->T[j]->dofs[0]];                                      
            idx2[1] = Mesh2->P1P->globaldofs[Mesh2->PD1->T[j]->dofs[1]];                                      
            idx2[2] = Mesh2->P1P->globaldofs[Mesh2->PD1->T[j]->dofs[2]];                                      
            idx2[3] = Mesh2->P1P->globaldofs[Mesh2->PD1->T[j]->T[0]->dofs[2]];                                
            idx2[4] = Mesh2->P1P->globaldofs[Mesh2->PD1->T[j]->T[0]->dofs[1]];                                
            idx2[5] = Mesh2->P1P->globaldofs[Mesh2->PD1->T[j]->T[1]->dofs[2]];                                                                                                                                            idx2[6] = Mesh2->P1P->globaldofs[Mesh2->PD1->T[j]->T[3]->dofs[2]];
            gsl_vector_complex_set(f2, 0, gsl_complex_rect(1.0/(Mesh2->PD1->lensbf[j][0]),0));                
            gsl_vector_complex_set(f2, 1, gsl_complex_rect(1.0/(Mesh2->PD1->lensbf[j][1]),0));                
            gsl_vector_complex_set(f2, 2, gsl_complex_rect(1.0/(Mesh2->PD1->lensbf[j][2]),0));                
                                                                                                              
            for( i = 0; i < Mesh1->PD1->nt; i++ )                                                             
            {                                                                                                 
                idx1[0] = Mesh1->P1P->globaldofs[Mesh1->PD1->T[i]->dofs[0]];                                  
                idx1[1] = Mesh1->P1P->globaldofs[Mesh1->PD1->T[i]->dofs[1]];                                  
                idx1[2] = Mesh1->P1P->globaldofs[Mesh1->PD1->T[i]->dofs[2]];                                  
                idx1[3] = Mesh1->P1P->globaldofs[Mesh1->PD1->T[i]->T[0]->dofs[2]];                            
                idx1[4] = Mesh1->P1P->globaldofs[Mesh1->PD1->T[i]->T[0]->dofs[1]];                            
                idx1[5] = Mesh1->P1P->globaldofs[Mesh1->PD1->T[i]->T[1]->dofs[2]];                            
                idx1[6] = Mesh1->P1P->globaldofs[Mesh1->PD1->T[i]->T[3]->dofs[2]];                            
                gsl_vector_complex_set(f2, 0, gsl_complex_rect(1.0/(Mesh1->PD1->lensbf[i][0]),0));            
                gsl_vector_complex_set(f2, 1, gsl_complex_rect(1.0/(Mesh1->PD1->lensbf[i][1]),0));            
                gsl_vector_complex_set(f2, 2, gsl_complex_rect(1.0/(Mesh1->PD1->lensbf[i][2]),0));            
                                                                                                              
                for(k=0;k<7;k++)                                                                              
                {                                                                                             
                    for(l=0;l<7;l++)                                                                          
                    {                                                                                         
                        if(GSL_REAL(gsl_matrix_complex_get(Mat,idx1[k], idx2[l]))==0 & GSL_IMAG(gsl_matrix_complex_get(Mat,idx1[k], idx2[l]))==0)
                        {                                                                                     
                                sing[0] = 0;                                                                  
                                sing[1] = 0;                                                                  
                                aux3 = assemble_P1_HP_entries2(order, kappa, i, j, idx1[k], idx2[l],sing, Mesh1,Mesh2, k, l);
                                aux3 = gsl_complex_add(aux3,gsl_complex_rect(sing[0],sing[1]));               
                                gsl_vector_complex_set(f3, l, aux3);                                          
                                gsl_matrix_complex_set(Mat,idx1[k], idx2[l],aux3);                            
                         }                                                                                    
                        else                                                                                  
                        {                                                                                     
                           gsl_vector_complex_set(f3, l, gsl_matrix_complex_get(Mat,idx1[k], idx2[l]));       
                        }                                                                                     
                    }                                                                                         
                    gsl_blas_zdotu(f3,f2,&aux1);
                    gsl_vector_complex_set(f4, k, aux1);                                                      
                }                                                                                             
                                                                                                              
                gsl_blas_zdotu(f4,f1,&aux2);                                                                  
                gsl_vector_complex_set(&auxm.vector, i,aux2);                                                  
            }  
        }                                                                                                                                                                                            
                                                                                                              
}  

gsl_matrix_complex *assemble_HLP_f(gsl_matrix_complex* Mat,                                                                  
                   pMeshes Mesh1,                                                                             
                   pMeshes Mesh2,                                                                             
                   int *order,                                                                                 
                   double kappa,                                                                              
                   // Order of quadrature                                                                     
                   // Output                                                                                  
                   int opt                                                                        
                  )
                {                                                                                   
        int i, j, k, l;                                                                                           
        int ndofb;                                                                                                
                                                                                                              
      
        gsl_matrix_complex *M = gsl_matrix_complex_alloc( Mesh1->PD1->nt, Mesh2->PD1->nt);                        
        gsl_matrix_complex_set_zero(M);   
        gsl_complex aux1, aux2, aux3;                                                                         
        aux1 = gsl_complex_rect(0.0,0.0);                                                                     
        aux2 = gsl_complex_rect(0.0,0.0);                                                                     
        aux3 = gsl_complex_rect(0.0,0.0);                                                                     
        double sing[2];                                                                                       
        int idx1[7], idx2[7], idx3[7], idx4[7];                                                                                 
        gsl_vector_complex *f1 = gsl_vector_complex_alloc(7);                                                 
        gsl_vector_complex *f2 = gsl_vector_complex_alloc(7);                                                 
        gsl_vector_complex *f3 = gsl_vector_complex_alloc(7);                                                 
        gsl_vector_complex *f4 = gsl_vector_complex_alloc(7);                                                 
        gsl_vector_complex_set(f1, 3, gsl_complex_rect(1.0,0));                                               
        gsl_vector_complex_set(f1, 4, gsl_complex_rect(0.5,0));                                               
        gsl_vector_complex_set(f1, 5, gsl_complex_rect(0.5,0));                                               
        gsl_vector_complex_set(f1, 6, gsl_complex_rect(0.5,0));                                               
        gsl_vector_complex_set(f2, 3, gsl_complex_rect(1.0,0));                                               
        gsl_vector_complex_set(f2, 4, gsl_complex_rect(0.5,0));                                               
        gsl_vector_complex_set(f2, 5, gsl_complex_rect(0.5,0));                                               
        gsl_vector_complex_set(f2, 6, gsl_complex_rect(0.5,0));                                               
    
        for( i = 0; i < Mesh1->PD1->nt; i++ )                                                                  
        {  
            idx1[0] = Mesh1->P1P->globaldofs[Mesh1->PD1->T[i]->dofs[0]];                                        
            idx1[1] = Mesh1->P1P->globaldofs[Mesh1->PD1->T[i]->dofs[1]];                                        
            idx1[2] = Mesh1->P1P->globaldofs[Mesh1->PD1->T[i]->dofs[2]];                                        
            idx1[3] = Mesh1->P1P->globaldofs[Mesh1->PD1->T[i]->T[0]->dofs[2]];                                  
            idx1[4] = Mesh1->P1P->globaldofs[Mesh1->PD1->T[i]->T[0]->dofs[1]];                                  
            idx1[5] = Mesh1->P1P->globaldofs[Mesh1->PD1->T[i]->T[5]->dofs[2]];
            idx1[6] = Mesh1->P1P->globaldofs[Mesh1->PD1->T[i]->T[3]->dofs[2]];
            idx3[0] = Mesh1->P1P->sbf[idx1[0]]->dv;                                      
            idx3[1] = Mesh1->P1P->sbf[idx1[1]]->dv;                                      
            idx3[2] = Mesh1->P1P->sbf[idx1[2]]->dv;                                      
            idx3[3] = Mesh1->P1P->sbf[idx1[3]]->dv;                                
            idx3[4] = Mesh1->P1P->sbf[idx1[4]]->dv;                                
            idx3[5] = Mesh1->P1P->sbf[idx1[5]]->dv;
            idx3[6] = Mesh1->P1P->sbf[idx1[6]]->dv;
            gsl_vector_complex_set(f1, 0, gsl_complex_rect(1.0/(Mesh1->PD1->lensbf[i][0]),0));                 
            gsl_vector_complex_set(f1, 1, gsl_complex_rect(1.0/(Mesh1->PD1->lensbf[i][1]),0));                 
            gsl_vector_complex_set(f1, 2, gsl_complex_rect(1.0/(Mesh1->PD1->lensbf[i][2]),0));                                                                                                       
            for( j = 0; j < Mesh2->PD1->nt; j++ )                                                              
            {                                                                                                 
                idx2[0] = Mesh2->P1P->globaldofs[Mesh2->PD1->T[j]->dofs[0]];                                    
                idx2[1] = Mesh2->P1P->globaldofs[Mesh2->PD1->T[j]->dofs[1]];                                    
                idx2[2] = Mesh2->P1P->globaldofs[Mesh2->PD1->T[j]->dofs[2]];                                    
                idx2[3] = Mesh2->P1P->globaldofs[Mesh2->PD1->T[j]->T[0]->dofs[2]];                              
                idx2[4] = Mesh2->P1P->globaldofs[Mesh2->PD1->T[j]->T[0]->dofs[1]];                              
                idx2[5] = Mesh2->P1P->globaldofs[Mesh2->PD1->T[j]->T[5]->dofs[2]];                              
                idx2[6] = Mesh2->P1P->globaldofs[Mesh2->PD1->T[j]->T[3]->dofs[2]];                              
                idx4[0] = Mesh2->P1P->sbf[idx2[0]]->dv;                                                           
                idx4[1] = Mesh2->P1P->sbf[idx2[1]]->dv;                                                           
                idx4[2] = Mesh2->P1P->sbf[idx2[2]]->dv;                                                           
                idx4[3] = Mesh2->P1P->sbf[idx2[3]]->dv;                                                         
                idx4[4] = Mesh2->P1P->sbf[idx2[4]]->dv;                                                           
                idx4[5] = Mesh2->P1P->sbf[idx2[5]]->dv;                                                           
                idx4[6] = Mesh2->P1P->sbf[idx2[6]]->dv;
                gsl_vector_complex_set(f2, 0, gsl_complex_rect(1.0/(Mesh2->PD1->lensbf[j][0]),0));             
                gsl_vector_complex_set(f2, 1, gsl_complex_rect(1.0/(Mesh2->PD1->lensbf[j][1]),0));             
                gsl_vector_complex_set(f2, 2, gsl_complex_rect(1.0/(Mesh2->PD1->lensbf[j][2]),0));             
               
               for(k=0;k<7;k++)                                                                              
                {                                                                                             
                    for(l=0;l<7;l++)                                                                          
                    {   
                        if(GSL_REAL(gsl_matrix_complex_get(Mat,idx3[k], idx4[l]))==0 && GSL_IMAG(gsl_matrix_complex_get(Mat,idx3[k], idx4[l]))==0 && (Mesh1->P1P->sbf[idx1[k]]->len >0 &&Mesh2->P1P->sbf[idx2[l]]->len>0))
                        {                                                                                     
                                sing[0] = 0;                                                                  
                                sing[1] = 0;                                                                  
                                aux3 = assemble_P1_HP_entries2(order, kappa, i, j, idx1[k], idx2[l],sing, Mesh1,Mesh2, k, l);
                                aux3 = gsl_complex_add(aux3,gsl_complex_rect(sing[0],sing[1]));               
                                gsl_vector_complex_set(f3, l, aux3);                                          
                                gsl_matrix_complex_set(Mat,idx3[k], idx4[l],aux3);
                         }                                                                                     
                        else if (GSL_REAL(gsl_matrix_complex_get(Mat,idx3[k], idx4[l]))!=0 && GSL_IMAG(gsl_matrix_complex_get(Mat,idx3[k], idx4[l]))!=0 && (Mesh1->P1P->sbf[idx1[k]]->len >0 &&Mesh2->P1P->sbf[idx2[l]]->len>0))                                                                                 
                        {                                                                                     
                           gsl_vector_complex_set(f3, l, gsl_matrix_complex_get(Mat,idx3[k], idx4[l]));       
                        }                                                             
                        else
                        {
                            gsl_vector_complex_set(f3, l, gsl_complex_rect(0.0,0.0));
                        }
                    }                                                                   
                    gsl_blas_zdotu(f3,f2,&aux1);
                    gsl_vector_complex_set(f4, k, aux1);                                                      
                }                                                                                             
                                                                                                              
                gsl_blas_zdotu(f4,f1,&aux2); 
                gsl_matrix_complex_set(M, i, j,aux2);    
            }                                                                                                 
        }                                                                                                                                                    
                                                                                                    
    return M;                                                                                                      

}


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
                    ){
    
    double weights[16], eta[4], chi[4], phij[16], phii[16], di[16], dj[16];
    double P0_x[3*order[0]*order[0]], P0_y[3*order[0]*order[0]], IE_x[order[0]*order[0]], IE_y[order[0]*order[0]], cres[2];
    int sz[2], i, j, k, l, length_arr;
    int start1 = row->start;
    int start2 = col->start;
    pMeshes Mesh1 = row->M;
    pMeshes Mesh2 = col->M;
     
    if( strcmp(testbf, "P0d" ) == 0 && strcmp(trialbf, "P0d") == 0 )
    {
        quad_rules_DUAL(eta, chi, sz, order[0]);                                                                     
        weights_DUAL(weights, order[0], sz[0]);                                                                      
        length_arr = order[0]*order[0];
        if( strcmp(mode, "row" ) == 0)
        {
            gsl_vector_complex_view auxm = gsl_matrix_complex_row(m, it);
            gsl_vector_complex_set_zero(&auxm.vector);
            i = urow;
            for( j= 0; j < Mesh2->PD0->ndof; j++ )
            {
                if( index[i+start1] != index[j+start2] )
                {
                    for( k = 0; k <Mesh1->PD0->sbf[i]->len; k++ )
                    {
                        integration_elements_DUAL( sz[0],sz[1], eta, chi,Mesh1->PD0->sbf[i]->Q[k]->v, IE_x, P0_x );
                        ev_basis_function( length_arr, 0, phii, eta, chi, trialbf );
                        vec_prod( length_arr, phii, weights, IE_x, di );
                        
                        for(l = 0; l < Mesh2->PD0->sbf[j]->len; l++)
                        {
                            integration_elements_DUAL( sz[0], sz[1], eta, chi, Mesh2->PD0->sbf[j]->Q[l]->v, IE_y, P0_y );
                            ev_basis_function( length_arr, 0, phij, eta, chi, testbf );
                            vec_prod( length_arr, phij, weights, IE_y, dj );
                            evaluate_kernel( P0_x,P0_y, kappa, length_arr, length_arr, di, dj, cres ); 
                            gsl_vector_complex_set(&auxm.vector, j, gsl_complex_add(gsl_vector_complex_get(&auxm.vector,j),gsl_complex_rect(cres[0],cres[1])));
        
                        }
                    }
                }
            }

        }
        else
        {
            gsl_vector_complex_view auxm = gsl_matrix_complex_row(m, it);                                     
            gsl_vector_complex_set_zero(&auxm.vector); 
            j = ucol;
            for( i = 0; i < Mesh1->PD0->ndof; i++ )
            {
                if( index[i+start1] != index[j+start2] )
                {
                    for( k = 0; k <Mesh1->PD0->sbf[i]->len; k++ )
                    {
                        integration_elements_DUAL( sz[0],sz[1], eta, chi,Mesh1->PD0->sbf[i]->Q[k]->v, IE_x, P0_x );
                        ev_basis_function( length_arr, 0, phii, eta, chi, trialbf );
                        vec_prod( length_arr, phii, weights, IE_x, di );
                        
                        for(l = 0; l < Mesh2->PD0->sbf[j]->len; l++)
                        {
                            integration_elements_DUAL( sz[0], sz[1], eta, chi, Mesh2->PD0->sbf[j]->Q[l]->v, IE_y, P0_y );
                            ev_basis_function( length_arr, 0, phij, eta, chi, testbf );
                            vec_prod( length_arr, phij, weights, IE_y, dj );
                            evaluate_kernel( P0_x,P0_y, kappa, length_arr, length_arr, di, dj, cres );
                            gsl_vector_complex_set(&auxm.vector, i, gsl_complex_add(gsl_vector_complex_get(&auxm.vector,i),gsl_complex_rect(cres[0],cres[1])));
                        }
                    }
                }
            }
     
        }

    }
    /*else if (strcmp(testbf, "P1d" ) == 0 && strcmp(trialbf, "P1d") == 0 )
    {
        double eta2[4], chi2[4], weights2[4];                                                                 
        int sz2[2];                                                                                           
        int order3 = 2;                                                                                       
        int order2 =1;                                                                                            
        quad_rules_PRIMAL(eta, chi, sz, order);                                                               
        quad_rules_PRIMAL(eta2, chi2, sz2, order2);                                                           
        weights_PRIMAL(weights, order, sz[0]);                                                                
        weights_PRIMAL(weights2, order2, sz2[0]);                                                             
        length_arr = sz[0];
        if( strcmp(mode, "row" ) == 0)                                                                        
        {
            gsl_vector_complex_view auxm = gsl_matrix_complex_row(m, it);                                     
            gsl_vector_complex_set_zero(&auxm.vector);                                                        
            i = urow;
                                                                                                    
            for( j = 0; j < Mesh2->P1P->ndof; j++ )                                                           
            {                                                                                                 
                                                                                                              
                for( k = 0; k <Mesh1->P1P->sbf[i]->len; k++ )                                                 
                {                                                                                             
                    integration_elements_PRIMAL( sz[0],sz[1], eta, chi, Mesh1->P1P->sbf[i]->T[k]->v, IE_x, P0_x );
                                                                                                              
                    ev_basis_function( length_arr, 0, phii, eta, chi, trialbf ); 
            
                    vec_prod( length_arr, phii, weights, IE_x, di );                                          
                                                                                                              
                    for(l = 0; l < Mesh2->P1P->sbf[j]->len; l++)                                              
                    {   
                        if( Mesh1->P1P->sbf[i]->num[k] != Mesh2->P1P->sbf[j]->num[l] )                        
                        {
                            integration_elements_PRIMAL( sz[0],sz[1], eta, chi , Mesh2->P1P->sbf[j]->T[l]->v, IE_y, P0_y );
                            ev_basis_function( length_arr, 0, phij, eta, chi, testbf );                       
                            vec_prod( length_arr, phij, weights, IE_y, dj );                                  
                            evaluate_kernel( P0_x,P0_y, kappa, length_arr, length_arr, di, dj, cres );        
                            gsl_vector_complex_set(&auxm.vector, j, gsl_complex_add(gsl_vector_complex_get(&auxm.vector,j),gsl_complex_rect(cres[0],cres[1]))); 
                        }
                        else                                                                                  
                        {                                                                                     
                            integration_elements_PRIMAL( sz2[0],sz2[1], eta2, chi2, Mesh2->P1P->sbf[j]->T[l]->v, IE_x, P0_x );
                            singular(P0_x, Mesh2->P1P->sbf[j]->T[l]->v, cres, IE_x, weights2,sz2, area(Mesh2->P1P->sbf[j]->T[l]->v), order3, 0, kappa, "P1");
                            gsl_vector_complex_set(&auxm.vector, j, gsl_complex_add(gsl_vector_complex_get(&auxm.vector,j),gsl_complex_rect(cres[0],cres[1])));
                        }
                   
                    }                         
                }                                                                                             
            }                                             

        }
        else
        {
            gsl_vector_complex_view auxm = gsl_matrix_complex_row(m, it);                                     
            gsl_vector_complex_set_zero(&auxm.vector);                                                        
            j = ucol;                                                                                            
            for( i = 0; i < Mesh1->P1P->ndof; i++ )                                                           
            {                                                                                                 
                                                                                                              
                for( k = 0; k <Mesh1->P1P->sbf[i]->len; k++ )                                                 
                {                                                                                             
                    integration_elements_PRIMAL( sz[0],sz[1], eta, chi, Mesh1->P1P->sbf[i]->T[k]->v, IE_x, P0_x );
                                                                                                              
                    ev_basis_function( length_arr, 0, phii, eta, chi, trialbf );                              
                    vec_prod( length_arr, phii, weights, IE_x, di );                                                                                  
                    for(l = 0; l < Mesh2->P1P->sbf[j]->len; l++)                                              
                    {    
                        if( Mesh1->P1P->sbf[i]->num[k] != Mesh2->P1P->sbf[j]->num[l] )                        
                        {                                                                                                                                                                       
                            integration_elements_PRIMAL( sz[0],sz[1], eta, chi , Mesh2->P1P->sbf[j]->T[l]->v, IE_y, P0_y );
                            ev_basis_function( length_arr, 0, phij, eta, chi, testbf );                       
                            vec_prod( length_arr, phij, weights, IE_y, dj );                                  
                            evaluate_kernel( P0_x,P0_y, kappa, length_arr, length_arr, di, dj, cres );        
                            gsl_vector_complex_set(&auxm.vector, i, gsl_complex_add(gsl_vector_complex_get(&auxm.vector,i),gsl_complex_rect(cres[0],cres[1])));    
                        }
                        else                                                                                  
                        {                                                                                     
                            integration_elements_PRIMAL( sz2[0],sz2[1], eta2, chi2, Mesh2->P1P->sbf[j]->T[l]->v, IE_x, P0_x );
                            singular(P0_x, Mesh2->P1P->sbf[j]->T[l]->v, cres, IE_x, weights2,sz2, area(Mesh2->P1P->sbf[j]->T[l]->v), order3, 0, kappa, "P1");
                            gsl_vector_complex_set(&auxm.vector, i, gsl_complex_add(gsl_vector_complex_get(&auxm.vector,i),gsl_complex_rect(cres[0],cres[1])));
                        }
                    }
                }                                                                                             
            }                                                                                                 
        }
    
    }*/
    else if (strcmp(testbf, "P1d" ) == 0 && strcmp(trialbf, "P1d") == 0 )                                   
    {                                                                                                         
        double eta2[4], chi2[4], weights2[4];                                                                 
        int sz2[2];                                                                                           
        int order3 = 2;                                                                                       
        int order2 =1;                                                                                            
        quad_rules_PRIMAL(eta, chi, sz, order);                                                               
        quad_rules_PRIMAL(eta2, chi2, sz2, order2);                                                           
        weights_PRIMAL(weights, order, sz[0]);                                                                
        weights_PRIMAL(weights2, order2, sz2[0]);                                                             
        length_arr = sz[0];  
        gsl_complex aux1, aux2;
        int idx1[7], idx2[7];                                                                                     
        gsl_vector_complex *f1 = gsl_vector_complex_alloc(7);                                                     
        gsl_vector_complex *f2 = gsl_vector_complex_alloc(7);                                                     
        gsl_vector_complex *f3 = gsl_vector_complex_alloc(7);                                                     
        gsl_vector_complex *f4 = gsl_vector_complex_alloc(7);
        gsl_vector_complex_set(f1, 3, gsl_complex_rect(1.0,0));
        gsl_vector_complex_set(f1, 4, gsl_complex_rect(0.5,0));                                       
        gsl_vector_complex_set(f1, 5, gsl_complex_rect(0.5,0));                                       
        gsl_vector_complex_set(f1, 6, gsl_complex_rect(0.5,0));
        gsl_vector_complex_set(f2, 3, gsl_complex_rect(1.0,0));                                               
        gsl_vector_complex_set(f2, 4, gsl_complex_rect(0.5,0));                                               
        gsl_vector_complex_set(f2, 5, gsl_complex_rect(0.5,0));                                               
        gsl_vector_complex_set(f2, 6, gsl_complex_rect(0.5,0)); 
        if( strcmp(mode, "row" ) == 0)                                                                        
        {                                                                                                     
            gsl_vector_complex_view auxm = gsl_matrix_complex_row(m, it);                                     
            gsl_vector_complex_set_zero(&auxm.vector);                                                        
            i = urow; 
            idx1[0] = Mesh1->P1P->globaldofs[Mesh1->PD1->T[i]->dofs[0]];                                  
            idx1[1] = Mesh1->P1P->globaldofs[Mesh1->PD1->T[i]->dofs[1]];                                  
            idx1[2] = Mesh1->P1P->globaldofs[Mesh1->PD1->T[i]->dofs[2]];                                  
            idx1[3] = Mesh1->P1P->globaldofs[Mesh1->PD1->T[i]->T[0]->dofs[1]];                            
            idx1[4] = Mesh1->P1P->globaldofs[Mesh1->PD1->T[i]->T[0]->dofs[2]];                            
            idx1[5] = Mesh1->P1P->globaldofs[Mesh1->PD1->T[i]->T[1]->dofs[1]];                            
            idx1[6] = Mesh1->P1P->globaldofs[Mesh1->PD1->T[i]->T[3]->dofs[1]];
            gsl_vector_complex_set(f1, 0, gsl_complex_rect(1.0/(2*Mesh1->PD1->lensbf[i][0]),0));          
            gsl_vector_complex_set(f1, 1, gsl_complex_rect(1.0/(2*Mesh1->PD1->lensbf[i][1]),0));          
            gsl_vector_complex_set(f1, 2, gsl_complex_rect(1.0/(2*Mesh1->PD1->lensbf[i][2]),0));
            for( j = 0; j < Mesh2->PD1->nt; j++ )                                                                 
            {     
                idx2[0] = Mesh2->P1P->globaldofs[Mesh2->PD1->T[j]->dofs[0]];                                  
                idx2[1] = Mesh2->P1P->globaldofs[Mesh2->PD1->T[j]->dofs[1]];                                  
                idx2[2] = Mesh2->P1P->globaldofs[Mesh2->PD1->T[j]->dofs[2]];                                  
                idx2[3] = Mesh2->P1P->globaldofs[Mesh2->PD1->T[j]->T[0]->dofs[1]];                            
                idx2[4] = Mesh2->P1P->globaldofs[Mesh2->PD1->T[j]->T[0]->dofs[2]];                            
                idx2[5] = Mesh2->P1P->globaldofs[Mesh2->PD1->T[j]->T[1]->dofs[1]];                            
                idx2[6] = Mesh2->P1P->globaldofs[Mesh2->PD1->T[j]->T[3]->dofs[1]];                             
                gsl_vector_complex_set(f2, 0, gsl_complex_rect(1.0/(2*Mesh2->PD1->lensbf[j][0]),0));          
                gsl_vector_complex_set(f2, 1, gsl_complex_rect(1.0/(2*Mesh2->PD1->lensbf[j][1]),0));          
                gsl_vector_complex_set(f2, 2, gsl_complex_rect(1.0/(2*Mesh2->PD1->lensbf[j][2]),0));          
                for(k=0;k<7;k++)                                                                              
                {                                                                                             
                    for(l=0;l<7;l++)                                                                          
                    {   
                        if(GSL_REAL(gsl_matrix_complex_get(M,idx1[k], idx2[l]))==0 & GSL_IMAG(gsl_matrix_complex_get(M,idx1[k], idx2[l]))==0)
                        {
                            gsl_vector_complex_set(f3, l,assemble_P1_entries2(order,kappa,idx1[k], idx2[l], Mesh1,Mesh2));
                            gsl_matrix_complex_set(M,idx1[k], idx2[l],gsl_vector_complex_get(f3,l));
                        }
                        else
                        {   
                            gsl_vector_complex_set(f3, l, gsl_matrix_complex_get(M,idx1[k], idx2[l]));
                        }
                    }                                                                                         
                    gsl_blas_zdotu(f3,f2,&aux1);                                                              
                    gsl_vector_complex_set(f4, k, aux1);                                                      
                }                                                                                             
                                                                                                              
                gsl_blas_zdotu(f4,f1,&aux2);
                gsl_vector_complex_set(&auxm.vector, j,aux2);    
            }
        }
        else                                                                      
        {                                                                                                     
            gsl_vector_complex_view auxm = gsl_matrix_complex_row(m, it);                                     
            gsl_vector_complex_set_zero(&auxm.vector);                                                        
            j = ucol; 
            idx2[0] = Mesh2->P1P->globaldofs[Mesh2->PD1->T[j]->dofs[0]];                                  
            idx2[1] = Mesh2->P1P->globaldofs[Mesh2->PD1->T[j]->dofs[1]];                                  
            idx2[2] = Mesh2->P1P->globaldofs[Mesh2->PD1->T[j]->dofs[2]];                                  
            idx2[3] = Mesh2->P1P->globaldofs[Mesh2->PD1->T[j]->T[0]->dofs[1]];                            
            idx2[4] = Mesh2->P1P->globaldofs[Mesh2->PD1->T[j]->T[0]->dofs[2]];                            
            idx2[5] = Mesh2->P1P->globaldofs[Mesh2->PD1->T[j]->T[1]->dofs[1]];                            
            idx2[6] = Mesh2->P1P->globaldofs[Mesh2->PD1->T[j]->T[3]->dofs[1]];                                
            gsl_vector_complex_set(f2, 0, gsl_complex_rect(1.0/(2*Mesh2->PD1->lensbf[j][0]),0));          
            gsl_vector_complex_set(f2, 1, gsl_complex_rect(1.0/(2*Mesh2->PD1->lensbf[j][1]),0));          
            gsl_vector_complex_set(f2, 2, gsl_complex_rect(1.0/(2*Mesh2->PD1->lensbf[j][2]),0));                                              
            for( i = 0; i < Mesh1->PD1->nt; i++ )                                                             
            {                             
                idx1[0] = Mesh1->P1P->globaldofs[Mesh1->PD1->T[i]->dofs[0]];                                       
                idx1[1] = Mesh1->P1P->globaldofs[Mesh1->PD1->T[i]->dofs[1]];                                       
                idx1[2] = Mesh1->P1P->globaldofs[Mesh1->PD1->T[i]->dofs[2]];                                       
                idx1[3] = Mesh1->P1P->globaldofs[Mesh1->PD1->T[i]->T[0]->dofs[1]];                                 
                idx1[4] = Mesh1->P1P->globaldofs[Mesh1->PD1->T[i]->T[0]->dofs[2]];                                 
                idx1[5] = Mesh1->P1P->globaldofs[Mesh1->PD1->T[i]->T[1]->dofs[1]];                                 
                idx1[6] = Mesh1->P1P->globaldofs[Mesh1->PD1->T[i]->T[3]->dofs[1]];                                 
                gsl_vector_complex_set(f1, 0, gsl_complex_rect(1.0/(2*Mesh1->PD1->lensbf[i][0]),0));               
                gsl_vector_complex_set(f1, 1, gsl_complex_rect(1.0/(2*Mesh1->PD1->lensbf[i][1]),0));               
                gsl_vector_complex_set(f1, 2, gsl_complex_rect(1.0/(2*Mesh1->PD1->lensbf[i][2]),0));                                                          
                                                                                                  
                for(k=0;k<7;k++)                                                                                   
                {                                                                                                  
                    for(l=0;l<7;l++)                                                                              
                    {                                                                                             
                        if(GSL_REAL(gsl_matrix_complex_get(M,idx1[k], idx2[l]))==0 & GSL_IMAG(gsl_matrix_complex_get(M,idx1[k], idx2[l]))==0)
                        {                                                                                     
                            gsl_vector_complex_set(f3, l,assemble_P1_entries2(order,kappa,idx1[k], idx2[l], Mesh1,Mesh2));
                            gsl_matrix_complex_set(M,idx1[k], idx2[l],gsl_vector_complex_get(f3,l));                  
                        }                                                                                     
                        else                                                                                  
                        {                                                                                     
                            gsl_vector_complex_set(f3, l, gsl_matrix_complex_get(M,idx1[k], idx2[l]));        
                        } 
                    }                                                                                             
                    gsl_blas_zdotu(f3,f2,&aux1);                                                                  
                    gsl_vector_complex_set(f4, k, aux1);                                                          
                }                                                                                                      
                                                                                                              
                gsl_blas_zdotu(f4,f1,&aux2); 
                gsl_vector_complex_set(&auxm.vector, i,aux2);                                                                                                
            }
        }

    }
}

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
                    ){
    
        int sz[2], i, j, k, l, length_arr;                                                                        
        pMeshes Mesh1 = row->M;                                                                                   
        pMeshes Mesh2 = col->M;                                                                                          
        double eta2[4], chi2[4], weights2[4];                                                                 
        int sz2[2];                                                                                           
        gsl_complex aux1, aux2;                                                                               
        int idx1[7], idx2[7], idx3[7], idx4[7];
        double sing[2];
        gsl_complex aux3 = gsl_complex_rect(0.0,0.0);                                                                                 
        gsl_vector_complex *f1 = gsl_vector_complex_alloc(7);                                                 
        gsl_vector_complex *f2 = gsl_vector_complex_alloc(7);                                                 
        gsl_vector_complex *f3 = gsl_vector_complex_alloc(7);                                                 
        gsl_vector_complex *f4 = gsl_vector_complex_alloc(7);                                                 
        gsl_vector_complex_set(f1, 3, gsl_complex_rect(1.0,0));                                               
        gsl_vector_complex_set(f1, 4, gsl_complex_rect(0.5,0));                                               
        gsl_vector_complex_set(f1, 5, gsl_complex_rect(0.5,0));                                               
        gsl_vector_complex_set(f1, 6, gsl_complex_rect(0.5,0));                                               
        gsl_vector_complex_set(f2, 3, gsl_complex_rect(1.0,0));                                               
        gsl_vector_complex_set(f2, 4, gsl_complex_rect(0.5,0));                                               
        gsl_vector_complex_set(f2, 5, gsl_complex_rect(0.5,0));                                               
        gsl_vector_complex_set(f2, 6, gsl_complex_rect(0.5,0));                                               
        if( strcmp(mode, "row" ) == 0)                                                                        
        {                                                                                                     
            gsl_vector_complex_view auxm = gsl_matrix_complex_row(m, it);                                     
            gsl_vector_complex_set_zero(&auxm.vector);                                                        
            i = urow;                                                                                         
            idx1[0] = Mesh1->P1P->globaldofs[Mesh1->PD1->T[i]->dofs[0]];                                      
            idx1[1] = Mesh1->P1P->globaldofs[Mesh1->PD1->T[i]->dofs[1]];                                      
            idx1[2] = Mesh1->P1P->globaldofs[Mesh1->PD1->T[i]->dofs[2]];                                      
            idx1[3] = Mesh1->P1P->globaldofs[Mesh1->PD1->T[i]->T[0]->dofs[2]];                                
            idx1[4] = Mesh1->P1P->globaldofs[Mesh1->PD1->T[i]->T[0]->dofs[1]];                                
            idx1[5] = Mesh1->P1P->globaldofs[Mesh1->PD1->T[i]->T[5]->dofs[2]];                                
            idx1[6] = Mesh1->P1P->globaldofs[Mesh1->PD1->T[i]->T[3]->dofs[2]];  
            idx3[0] = Mesh1->P1P->sbf[idx1[0]]->dv;                                                           
            idx3[1] = Mesh1->P1P->sbf[idx1[1]]->dv;                                                           
            idx3[2] = Mesh1->P1P->sbf[idx1[2]]->dv;                                                           
            idx3[3] = Mesh1->P1P->sbf[idx1[3]]->dv;                                                           
            idx3[4] = Mesh1->P1P->sbf[idx1[4]]->dv;                                                           
            idx3[5] = Mesh1->P1P->sbf[idx1[5]]->dv;                                                           
            idx3[6] = Mesh1->P1P->sbf[idx1[6]]->dv;                               
            gsl_vector_complex_set(f1, 0, gsl_complex_rect(1.0/(Mesh1->PD1->lensbf[i][0]),0));              
            gsl_vector_complex_set(f1, 1, gsl_complex_rect(1.0/(Mesh1->PD1->lensbf[i][1]),0));              
            gsl_vector_complex_set(f1, 2, gsl_complex_rect(1.0/(Mesh1->PD1->lensbf[i][2]),0));              
            for( j = 0; j < Mesh2->PD1->nt; j++ )                                                             
            {                                                                                                 
                idx2[0] = Mesh2->P1P->globaldofs[Mesh2->PD1->T[j]->dofs[0]];                                  
                idx2[1] = Mesh2->P1P->globaldofs[Mesh2->PD1->T[j]->dofs[1]];                                  
                idx2[2] = Mesh2->P1P->globaldofs[Mesh2->PD1->T[j]->dofs[2]];                                  
                idx2[3] = Mesh2->P1P->globaldofs[Mesh2->PD1->T[j]->T[0]->dofs[2]];                            
                idx2[4] = Mesh2->P1P->globaldofs[Mesh2->PD1->T[j]->T[0]->dofs[1]];                            
                idx2[5] = Mesh2->P1P->globaldofs[Mesh2->PD1->T[j]->T[5]->dofs[2]];                            
                idx2[6] = Mesh2->P1P->globaldofs[Mesh2->PD1->T[j]->T[3]->dofs[2]]; 
                idx4[0] = Mesh2->P1P->sbf[idx2[0]]->dv;                                                       
                idx4[1] = Mesh2->P1P->sbf[idx2[1]]->dv;                                                       
                idx4[2] = Mesh2->P1P->sbf[idx2[2]]->dv;                                                       
                idx4[3] = Mesh2->P1P->sbf[idx2[3]]->dv;                                                       
                idx4[4] = Mesh2->P1P->sbf[idx2[4]]->dv;                                                       
                idx4[5] = Mesh2->P1P->sbf[idx2[5]]->dv;                                                       
                idx4[6] = Mesh2->P1P->sbf[idx2[6]]->dv;                            
                gsl_vector_complex_set(f2, 0, gsl_complex_rect(1.0/(Mesh2->PD1->lensbf[j][0]),0));          
                gsl_vector_complex_set(f2, 1, gsl_complex_rect(1.0/(Mesh2->PD1->lensbf[j][1]),0));          
                gsl_vector_complex_set(f2, 2, gsl_complex_rect(1.0/(Mesh2->PD1->lensbf[j][2]),0));          
                for(k=0;k<7;k++)                                                                              
                {                                                                                             
                    for(l=0;l<7;l++)                                                                          
                    {                                                                                         
                        if(GSL_REAL(gsl_matrix_complex_get(M,idx3[k], idx4[l]))==0 && GSL_IMAG(gsl_matrix_complex_get(M,idx3[k], idx4[l]))==0 && (Mesh1->P1P->sbf[idx1[k]]->len >0 &&Mesh2->P1P->sbf[idx2[l]]->len>0))
                        {                                                                                     
                                sing[0] = 0;                                                                  
                                sing[1] = 0;                                                                  
                                aux3 = assemble_P1_HP_entries2(order, kappa, i, j, idx1[k], idx2[l],sing, Mesh1,Mesh2, k, l);
                                aux3 = gsl_complex_add(aux3,gsl_complex_rect(sing[0],sing[1]));               
                                gsl_vector_complex_set(f3, l, aux3);                                          
                                gsl_matrix_complex_set(M,idx3[k], idx4[l],aux3);                            
                         }                                                                                    
                        else if (GSL_REAL(gsl_matrix_complex_get(M,idx3[k], idx4[l]))!=0 && GSL_IMAG(gsl_matrix_complex_get(M,idx3[k], idx4[l]))!=0 && (Mesh1->P1P->sbf[idx1[k]]->len >0 &&Mesh2->P1P->sbf[idx2[l]]->len>0))
                        {                                                                                     
                           gsl_vector_complex_set(f3, l, gsl_matrix_complex_get(M,idx3[k], idx4[l]));       
                        }                                                                                     
                        else                                                                                  
                        {                                                                                     
                            gsl_vector_complex_set(f3, l, gsl_complex_rect(0.0,0.0));                         
                        }                                                                                    
                    }                                                                                         
                    gsl_blas_zdotu(f3,f2,&aux1);                                                              
                    gsl_vector_complex_set(f4, k, aux1);                                                      
                }                                                                                             
                                                                                                              
                gsl_blas_zdotu(f4,f1,&aux2);                                                                  
                gsl_vector_complex_set(&auxm.vector, j,aux2);                                                 
            }                                                                                                 
        }                                                                                                     
        else                                                                                                  
        {                                                                                                     
            gsl_vector_complex_view auxm = gsl_matrix_complex_row(m, it);                                     
            gsl_vector_complex_set_zero(&auxm.vector);                                                        
            j = ucol;                                                                                         
            idx2[0] = Mesh2->P1P->globaldofs[Mesh2->PD1->T[j]->dofs[0]];                                      
            idx2[1] = Mesh2->P1P->globaldofs[Mesh2->PD1->T[j]->dofs[1]];                                      
            idx2[2] = Mesh2->P1P->globaldofs[Mesh2->PD1->T[j]->dofs[2]];                                      
            idx2[3] = Mesh2->P1P->globaldofs[Mesh2->PD1->T[j]->T[0]->dofs[2]];                                
            idx2[4] = Mesh2->P1P->globaldofs[Mesh2->PD1->T[j]->T[0]->dofs[1]];                                
            idx2[5] = Mesh2->P1P->globaldofs[Mesh2->PD1->T[j]->T[5]->dofs[2]];                                
            idx2[6] = Mesh2->P1P->globaldofs[Mesh2->PD1->T[j]->T[3]->dofs[2]];
            idx4[0] = Mesh2->P1P->sbf[idx2[0]]->dv;                                                       
            idx4[1] = Mesh2->P1P->sbf[idx2[1]]->dv;                                                       
            idx4[2] = Mesh2->P1P->sbf[idx2[2]]->dv;                                                       
            idx4[3] = Mesh2->P1P->sbf[idx2[3]]->dv;                                                       
            idx4[4] = Mesh2->P1P->sbf[idx2[4]]->dv;                                                       
            idx4[5] = Mesh2->P1P->sbf[idx2[5]]->dv;                                                       
            idx4[6] = Mesh2->P1P->sbf[idx2[6]]->dv;                                  
            gsl_vector_complex_set(f2, 0, gsl_complex_rect(1.0/(Mesh2->PD1->lensbf[j][0]),0));              
            gsl_vector_complex_set(f2, 1, gsl_complex_rect(1.0/(Mesh2->PD1->lensbf[j][1]),0));              
            gsl_vector_complex_set(f2, 2, gsl_complex_rect(1.0/(Mesh2->PD1->lensbf[j][2]),0));              
            for( i = 0; i < Mesh1->PD1->nt; i++ )                                                             
            {                                                                                                 
                idx1[0] = Mesh1->P1P->globaldofs[Mesh1->PD1->T[i]->dofs[0]];                                  
                idx1[1] = Mesh1->P1P->globaldofs[Mesh1->PD1->T[i]->dofs[1]];                                  
                idx1[2] = Mesh1->P1P->globaldofs[Mesh1->PD1->T[i]->dofs[2]];                                  
                idx1[3] = Mesh1->P1P->globaldofs[Mesh1->PD1->T[i]->T[0]->dofs[2]];                            
                idx1[4] = Mesh1->P1P->globaldofs[Mesh1->PD1->T[i]->T[0]->dofs[1]];                            
                idx1[5] = Mesh1->P1P->globaldofs[Mesh1->PD1->T[i]->T[5]->dofs[2]];                            
                idx1[6] = Mesh1->P1P->globaldofs[Mesh1->PD1->T[i]->T[3]->dofs[2]];                            
                idx3[0] = Mesh1->P1P->sbf[idx1[0]]->dv;                                                           
                idx3[1] = Mesh1->P1P->sbf[idx1[1]]->dv;                                                           
                idx3[2] = Mesh1->P1P->sbf[idx1[2]]->dv;                                                           
                idx3[3] = Mesh1->P1P->sbf[idx1[3]]->dv;                                                           
                idx3[4] = Mesh1->P1P->sbf[idx1[4]]->dv;                                                           
                idx3[5] = Mesh1->P1P->sbf[idx1[5]]->dv;                                                           
                idx3[6] = Mesh1->P1P->sbf[idx1[6]]->dv; 
                gsl_vector_complex_set(f1, 0, gsl_complex_rect(1.0/(Mesh1->PD1->lensbf[i][0]),0));          
                gsl_vector_complex_set(f1, 1, gsl_complex_rect(1.0/(Mesh1->PD1->lensbf[i][1]),0));          
                gsl_vector_complex_set(f1, 2, gsl_complex_rect(1.0/(Mesh1->PD1->lensbf[i][2]),0));          
                                                                                                              
                for(k=0;k<7;k++)                                                                              
                {                                                                                             
                    for(l=0;l<7;l++)                                                                          
                    {                                                                                         
                    
                        if(GSL_REAL(gsl_matrix_complex_get(M,idx3[k], idx4[l]))==0 && GSL_IMAG(gsl_matrix_complex_get(M,idx3[k], idx4[l]))==0 && (Mesh1->P1P->sbf[idx1[k]]->len >0 &&Mesh2->P1P->sbf[idx2[l]]->len>0))
                        {                                                                                     
                                sing[0] = 0;                                                                  
                                sing[1] = 0;                                                                  
                                aux3 = assemble_P1_HP_entries2(order, kappa, i, j, idx1[k], idx2[l],sing, Mesh1,Mesh2, k, l);
                                aux3 = gsl_complex_add(aux3,gsl_complex_rect(sing[0],sing[1]));               
                                gsl_vector_complex_set(f3, l, aux3);                                          
                                gsl_matrix_complex_set(M,idx3[k], idx4[l],aux3);                            
                         }                                                                                    
                        else if (GSL_REAL(gsl_matrix_complex_get(M,idx3[k], idx4[l]))!=0 && GSL_IMAG(gsl_matrix_complex_get(M,idx3[k], idx4[l]))!=0 && (Mesh1->P1P->sbf[idx1[k]]->len >0 &&Mesh2->P1P->sbf[idx2[l]]->len>0))
                        {                                                                                     
                           gsl_vector_complex_set(f3, l, gsl_matrix_complex_get(M,idx3[k], idx4[l]));       
                        }                                                                                     
                        else                                                                                  
                        {                                                                                     
                            gsl_vector_complex_set(f3, l, gsl_complex_rect(0.0,0.0));                         
                        }                                                                                    
                    }                                                                                         
                    gsl_blas_zdotu(f3,f2,&aux1);                                                              
                    gsl_vector_complex_set(f4, k, aux1);                                                      
                }                                                                                             
                                                                                                              
                gsl_blas_zdotu(f4,f1,&aux2);
                gsl_vector_complex_set(&auxm.vector, i,aux2);                                                 
            }                                                                                                 
        }                                                                                                     
                                                                                                              
                                                                                                             
} 

gsl_matrix_complex *average_matrix(                                                                           
                   pMeshes Mesh,                                                                               
                   // Order of quadrature                                                                     
                   // Output                                                                                  
                   int opt                                                                                    
                   ){                                                                                         
                                                                                                              
    int i,idx;                                                                                        
    int ndofb;                                                                                                
    if(opt==0)                                                                                                
    {                                                                                                         
        gsl_matrix_complex *M = gsl_matrix_complex_alloc( Mesh->PD1->nt,Mesh->P1P->ndof);                     
        gsl_matrix_complex_set_zero(M);                                                                       
        ndofb = Mesh->P1P->ndof;                                                                              
        for( i = 0; i < Mesh->PD1->nt; i++ )                                                                  
        {                                                                                                     
          gsl_matrix_complex_set(M, i, Mesh->P1P->globaldofs[Mesh->PD1->T[i]->dofs[0]], gsl_complex_rect(1.0/(2*Mesh->PD1->lensbf[i][0]),0));
          gsl_matrix_complex_set(M, i, Mesh->P1P->globaldofs[Mesh->PD1->T[i]->dofs[1]], gsl_complex_rect(1.0/(2*Mesh->PD1->lensbf[i][1]),0));
          gsl_matrix_complex_set(M, i, Mesh->P1P->globaldofs[Mesh->PD1->T[i]->dofs[2]], gsl_complex_rect(1.0/(2*Mesh->PD1->lensbf[i][2]),0));
          gsl_matrix_complex_set(M, i, Mesh->P1P->globaldofs[Mesh->PD1->T[i]->T[0]->dofs[1]], gsl_complex_rect(1.0,0));
          gsl_matrix_complex_set(M, i, Mesh->P1P->globaldofs[Mesh->PD1->T[i]->T[0]->dofs[2]], gsl_complex_rect(0.5,0));
          gsl_matrix_complex_set(M, i, Mesh->P1P->globaldofs[Mesh->PD1->T[i]->T[1]->dofs[1]], gsl_complex_rect(0.5,0));
          gsl_matrix_complex_set(M, i, Mesh->P1P->globaldofs[Mesh->PD1->T[i]->T[3]->dofs[1]], gsl_complex_rect(0.5,0));
        }                                                                                                     
        return M;                                                                                             
                                                                                                              
                                                                                                              
    }                                                                                                         
}


gsl_complex assemble_P1_entries(int *order, double kappa, int i, int j, pMeshes Mesh1, pMeshes Mesh2)
{

   double weights[16], eta[4], chi[4], phij[16], phii[16], di[16], dj[16];                                   
   double P0_x[3*order[0]*order[0]], P0_y[3*order[0]*order[0]], IE_x[order[0]*order[0]], IE_y[order[0]*order[0]], cres[2]; 
   double eta2[4], chi2[4], weights2[4];                                                                 
   int sz2[2], sz[2], length_arr;                                                                                           
   //int order3 = 2;                                                                                       
   //int order2 =1;
   gsl_complex aux = gsl_complex_rect(0.0,0.0);
   quad_rules_PRIMAL(eta, chi, sz, order[0]);                                                               
   quad_rules_PRIMAL(eta2, chi2, sz2, order[1]);                                                           
   weights_PRIMAL(weights, order, sz[0]);                                                                
   weights_PRIMAL(weights2, order[1], sz2[0]);                                                             
   length_arr = sz[0]; 
   int k, l;
   for( k = 0; k <Mesh1->P1P->sbf[i]->len; k++ )                                                 
   {                                                                                             
        integration_elements_PRIMAL( sz[0],sz[1], eta, chi, Mesh1->P1P->sbf[i]->T[k]->v, IE_x, P0_x );
                                                                                                              
        ev_basis_function( length_arr, 0, phii, eta, chi, "P1" );                              
        vec_prod( length_arr, phii, weights, IE_x, di );                                          
        for(l = 0; l < Mesh2->P1P->sbf[j]->len; l++)                                              
        {                                                                                         
           if( Mesh1->P1P->sbf[i]->num[k] != Mesh2->P1P->sbf[j]->num[l] )                        
           {                                                                                     
                integration_elements_PRIMAL( sz[0],sz[1], eta, chi , Mesh2->P1P->sbf[j]->T[l]->v, IE_y, P0_y );
                ev_basis_function( length_arr, 0, phij, eta, chi, "P1" );                       
                vec_prod( length_arr, phij, weights, IE_y, dj );                                  
                evaluate_kernel( P0_x,P0_y, kappa, length_arr, length_arr, di, dj, cres );        
                aux = gsl_complex_add(aux, gsl_complex_rect(cres[0],cres[1]));
           }                                                                                     
           else                                                                                  
           {                                                                                     
                integration_elements_PRIMAL( sz2[0],sz2[1], eta2, chi2, Mesh2->P1P->sbf[j]->T[l]->v, IE_x, P0_x );
                singular(P0_x, Mesh2->P1P->sbf[j]->T[l]->v, weights, cres, IE_x, sz2,  order[2], 0, kappa, "P1", 1.0);
                aux = gsl_complex_add(aux, gsl_complex_rect(cres[0],cres[1]));
           }                                                                                     
        }                                                                                         
     } 
     
     return aux;
}

gsl_complex assemble_P1_entries2(int *order, double kappa, int i, int j, pMeshes Mesh1, pMeshes Mesh2)
{
    
    double weights[16], eta[4], chi[4], phij[16], phii[16], di[16], dj[16];
    double P0_x[3*order[0]*order[0]], P0_y[3*order[0]*order[0]], IE_x[order[0]*order[0]], IE_y[order[0]*order[0]], cres[2];
    double eta2[4], chi2[4],eta3[4], chi3[4], weights2[4], weights3[4];
    double vert[12];
    int sz[2], sz2[2], sz3[2],length_arr, length_arr2;
    //int order3 = 2;
    //int order2 =1;
    gsl_complex aux = gsl_complex_rect(0.0,0.0);
    quad_rules_PRIMAL(eta, chi, sz, order[0]);
    quad_rules_PRIMAL(eta2, chi2, sz2, order[1]);
    quad_rules_DUAL(eta3, chi3, sz3, order[0]);
    
    weights_PRIMAL(weights, order[0], sz[0]);
    weights_PRIMAL(weights2, order[1], sz2[0]);
    weights_DUAL(weights3, order[0], sz[0]);
    
    length_arr = sz[0];
    length_arr2 = sz3[0]*sz3[0];
    int k, l;
    
    for( k = 0; k <Mesh1->P1P->sbf[i]->len; k=k+2 )
    {
        vert[0] = Mesh1->P1P->sbf[i]->T[k]->v[0];
        vert[1] = Mesh1->P1P->sbf[i]->T[k]->v[1];
        vert[2] = Mesh1->P1P->sbf[i]->T[k]->v[2];
        vert[3] = Mesh1->P1P->sbf[i]->T[k+1]->v[3];
        vert[4] = Mesh1->P1P->sbf[i]->T[k+1]->v[4];
        vert[5] = Mesh1->P1P->sbf[i]->T[k+1]->v[5];
        vert[6] = Mesh1->P1P->sbf[i]->T[k]->v[3];
        vert[7] = Mesh1->P1P->sbf[i]->T[k]->v[4];
        vert[8] = Mesh1->P1P->sbf[i]->T[k]->v[5];
        vert[9] = Mesh1->P1P->sbf[i]->T[k]->v[6];
        vert[10] = Mesh1->P1P->sbf[i]->T[k]->v[7];
        vert[11] = Mesh1->P1P->sbf[i]->T[k]->v[8];
        integration_elements_DUAL( sz3[0],sz3[1], eta3, chi3,vert, IE_x, P0_x );
        ev_basis_function( length_arr2, 0, phii, eta3, chi3, "P1di" );
        vec_prod( length_arr2, phii, weights3, IE_x, di );
        for(l = 0; l < Mesh2->P1P->sbf[j]->len; l=l+2)
        {
            
            if( Mesh1->P1P->sbf[i]->num[k] != Mesh2->P1P->sbf[j]->num[l] )
            {
                vert[0] = Mesh2->P1P->sbf[j]->T[l]->v[0];
                vert[1] = Mesh2->P1P->sbf[j]->T[l]->v[1];
                vert[2] = Mesh2->P1P->sbf[j]->T[l]->v[2];
                vert[3] = Mesh2->P1P->sbf[j]->T[l+1]->v[3];
                vert[4] = Mesh2->P1P->sbf[j]->T[l+1]->v[4];
                vert[5] = Mesh2->P1P->sbf[j]->T[l+1]->v[5];
                vert[6] = Mesh2->P1P->sbf[j]->T[l]->v[3];
                vert[7] = Mesh2->P1P->sbf[j]->T[l]->v[4];
                vert[8] = Mesh2->P1P->sbf[j]->T[l]->v[5];
                vert[9] = Mesh2->P1P->sbf[j]->T[l]->v[6];
                vert[10] = Mesh2->P1P->sbf[j]->T[l]->v[7];
                vert[11] = Mesh2->P1P->sbf[j]->T[l]->v[8];
                integration_elements_DUAL( sz3[0],sz3[1], eta3, chi3,vert, IE_y, P0_y );
                ev_basis_function( length_arr2, 0, phij, eta3, chi3, "P1di" );
                vec_prod( length_arr2, phij, weights3, IE_y, dj );
                evaluate_kernel( P0_x,P0_y, kappa, length_arr2, length_arr2, di, dj, cres );
                aux = gsl_complex_add(aux, gsl_complex_rect(cres[0],cres[1]));
            }
            else
            {
                integration_elements_PRIMAL( sz2[0],sz2[1], eta2, chi2, Mesh2->P1P->sbf[j]->T[l]->v, IE_x, P0_x );
                singular(P0_x, Mesh2->P1P->sbf[j]->T[l]->v, weights2, cres, IE_x, sz2,  order[2], 0, kappa, "P1", 1.0);
                aux = gsl_complex_add(aux, gsl_complex_rect(cres[0],cres[1]));
                
                integration_elements_PRIMAL( sz2[0],sz2[1], eta2, chi2, Mesh2->P1P->sbf[j]->T[l+1]->v, IE_x, P0_x );
                singular(P0_x, Mesh2->P1P->sbf[j]->T[l+1]->v, weights2, cres, IE_x, sz2,  order[2], 0, kappa, "P1", 1.0);
                aux = gsl_complex_add(aux, gsl_complex_rect(cres[0],cres[1]));
                
                integration_elements_PRIMAL( sz[0],sz[1], eta, chi, Mesh1->P1P->sbf[i]->T[k]->v, IE_x, P0_x );
                
                ev_basis_function( length_arr, 0, phii, eta, chi, "P1" );
                vec_prod( length_arr, phii, weights, IE_x, di );
                integration_elements_PRIMAL( sz[0],sz[1], eta, chi , Mesh2->P1P->sbf[j]->T[l+1]->v, IE_y, P0_y );
                ev_basis_function( length_arr, 0, phij, eta, chi, "P1" );
                vec_prod( length_arr, phij, weights, IE_y, dj );
                evaluate_kernel( P0_x,P0_y, kappa, length_arr, length_arr, di, dj, cres );
                aux = gsl_complex_add(aux, gsl_complex_rect(cres[0],cres[1]));
                
                integration_elements_PRIMAL( sz[0],sz[1], eta, chi, Mesh1->P1P->sbf[i]->T[k+1]->v, IE_x, P0_x );
                
                ev_basis_function( length_arr, 0, phii, eta, chi, "P1" );
                vec_prod( length_arr, phii, weights, IE_x, di );
                integration_elements_PRIMAL( sz[0],sz[1], eta, chi , Mesh2->P1P->sbf[j]->T[l]->v, IE_y, P0_y );
                ev_basis_function( length_arr, 0, phij, eta, chi, "P1" );
                vec_prod( length_arr, phij, weights, IE_y, dj );
                evaluate_kernel( P0_x,P0_y, kappa, length_arr, length_arr, di, dj, cres );
                aux = gsl_complex_add(aux, gsl_complex_rect(cres[0],cres[1]));
            }
        }
        
    }
    
    return aux;
}
gsl_matrix_complex *assemble_SLD_P1_f(
                    gsl_matrix_complex* Mat,                         
                   pMeshes Mesh1,
                   pMeshes Mesh2,
                   int *order,
                   double kappa,
                   // Order of quadrature                                                                     
                   // Output                                                                                  
                   int opt                                                                                    
                   ){                                                                                         
                                                                                                              
    int i, j, k, l;                                                                                        
    int ndofb;                                                                                                
                                    
    gsl_matrix_complex *M = gsl_matrix_complex_alloc( Mesh1->PD1->nt, Mesh2->PD1->nt);                     
    gsl_matrix_complex_set_zero(M);  
    gsl_complex aux1, aux2;                                                                     
    aux1 = gsl_complex_rect(0.0,0.0);
    aux2 = gsl_complex_rect(0.0,0.0); 

    int idx1[7], idx2[7];
    gsl_vector_complex *f1 = gsl_vector_complex_alloc(7);
    gsl_vector_complex *f2 = gsl_vector_complex_alloc(7);
    gsl_vector_complex *f3 = gsl_vector_complex_alloc(7);
    gsl_vector_complex *f4 = gsl_vector_complex_alloc(7);
    gsl_vector_complex_set(f1, 3, gsl_complex_rect(1.0,0));
    gsl_vector_complex_set(f1, 4, gsl_complex_rect(0.5,0));                                            
    gsl_vector_complex_set(f1, 5, gsl_complex_rect(0.5,0));                                            
    gsl_vector_complex_set(f1, 6, gsl_complex_rect(0.5,0));
    gsl_vector_complex_set(f2, 3, gsl_complex_rect(1.0,0));                                                   
    gsl_vector_complex_set(f2, 4, gsl_complex_rect(0.5,0));                                                   
    gsl_vector_complex_set(f2, 5, gsl_complex_rect(0.5,0));                                                   
    gsl_vector_complex_set(f2, 6, gsl_complex_rect(0.5,0));  
    for( i = 0; i < Mesh1->PD1->nt; i++ )                                          
    {   
        idx1[0] = Mesh1->P1P->globaldofs[Mesh1->PD1->T[i]->dofs[0]];
        idx1[1] = Mesh1->P1P->globaldofs[Mesh1->PD1->T[i]->dofs[1]];                                       
        idx1[2] = Mesh1->P1P->globaldofs[Mesh1->PD1->T[i]->dofs[2]];                                       
        idx1[3] = Mesh1->P1P->globaldofs[Mesh1->PD1->T[i]->T[0]->dofs[1]];                                 
        idx1[4] = Mesh1->P1P->globaldofs[Mesh1->PD1->T[i]->T[0]->dofs[2]];                                 
        idx1[5] = Mesh1->P1P->globaldofs[Mesh1->PD1->T[i]->T[1]->dofs[1]];                                 
        idx1[6] = Mesh1->P1P->globaldofs[Mesh1->PD1->T[i]->T[3]->dofs[1]]; 
        gsl_vector_complex_set(f1, 0, gsl_complex_rect(1.0/(2*Mesh1->PD1->lensbf[i][0]),0));               
        gsl_vector_complex_set(f1, 1, gsl_complex_rect(1.0/(2*Mesh1->PD1->lensbf[i][1]),0));               
        gsl_vector_complex_set(f1, 2, gsl_complex_rect(1.0/(2*Mesh1->PD1->lensbf[i][2]),0));
        for( j = 0; j < Mesh2->PD1->nt; j++ )                                                                     
        {
           idx2[0] = Mesh2->P1P->globaldofs[Mesh2->PD1->T[j]->dofs[0]];                                     
           idx2[1] = Mesh2->P1P->globaldofs[Mesh2->PD1->T[j]->dofs[1]];                                     
           idx2[2] = Mesh2->P1P->globaldofs[Mesh2->PD1->T[j]->dofs[2]];                                     
           idx2[3] = Mesh2->P1P->globaldofs[Mesh2->PD1->T[j]->T[0]->dofs[1]];                               
           idx2[4] = Mesh2->P1P->globaldofs[Mesh2->PD1->T[j]->T[0]->dofs[2]];                               
           idx2[5] = Mesh2->P1P->globaldofs[Mesh2->PD1->T[j]->T[1]->dofs[1]];                               
           idx2[6] = Mesh2->P1P->globaldofs[Mesh2->PD1->T[j]->T[3]->dofs[1]]; 
           gsl_vector_complex_set(f2, 0, gsl_complex_rect(1.0/(2*Mesh2->PD1->lensbf[j][0]),0));               
           gsl_vector_complex_set(f2, 1, gsl_complex_rect(1.0/(2*Mesh2->PD1->lensbf[j][1]),0));               
           gsl_vector_complex_set(f2, 2, gsl_complex_rect(1.0/(2*Mesh2->PD1->lensbf[j][2]),0));               
           for(k=0;k<7;k++)
           {
                for(l=0;l<7;l++)                                                                                   
                {            
                    if(GSL_REAL(gsl_matrix_complex_get(Mat,idx1[k], idx2[l]))==0 & GSL_IMAG(gsl_matrix_complex_get(Mat,idx1[k], idx2[l]))==0)
                    {                                                                                     
                        gsl_vector_complex_set(f3, l,assemble_P1_entries2(order,kappa,idx1[k], idx2[l], Mesh1,Mesh2));
                        gsl_matrix_complex_set(Mat,idx1[k], idx2[l],gsl_vector_complex_get(f3,l));                  
                    }                                                                                     
                    else                                                                                  
                    {                                                                                     
                        gsl_vector_complex_set(f3, l, gsl_matrix_complex_get(Mat,idx1[k], idx2[l]));        
                    }                                                                                 
                }
                gsl_blas_zdotu(f3,f2,&aux1);
                gsl_vector_complex_set(f4, k, aux1);
           }

           gsl_blas_zdotu(f4,f1,&aux2);      
           gsl_matrix_complex_set(M, i, j,aux2);

        }
    }                                                                                                
    return M;                                                                                             
                                                                                                             
}

gsl_complex assemble_P1_HP_entries2(int *order, double kappa, int i, int j, int idx, int idy, double *sing, pMeshes Mesh1, pMeshes Mesh2, int dofx, int dofy)
{
                                                                                                                 
    double weights_x[16], weights_y[16], etax[4], etay[4], eta2[4], chi2[4], chix[4], chiy[4];                                    
    int sz_x[2], sz_y[2], k, l, m;                                                                            
    double phii[16] = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};                                       
    double phij[16] = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};                                       
    gsl_complex res;                                                                                          
    gsl_vector_complex *di;                                                                                   
    gsl_vector_complex *dj;                                                                                   
    gsl_vector *n_x = gsl_vector_complex_alloc(3);                                                            
    gsl_vector *n_y = gsl_vector_complex_alloc(3);                                                            
    gsl_matrix *curl_x;                                                                                       
    gsl_matrix *curl_y;                                                                                       
    double cres[2];                                                                                           
    gsl_complex aux = gsl_complex_rect(0.0,0.0);                                                              
    gsl_complex singaux = gsl_complex_rect(0.0,0.0);                                                          
    gsl_complex singaux2 = gsl_complex_rect(0.0,0.0);                                                         
    int length_arr = order[0]*order[0];                                                                       
    double P0_x[3*length_arr], P0_y[3*length_arr], IE_x[length_arr], IE_y[length_arr];                        
    
    if((dofx==0 || dofx==1 || dofx==2 ) &&(dofy==0 || dofy==1 || dofy==2 ))                                   
    {     
        quad_rules_DUAL(etax, chix, sz_x, order[0]); 
        quad_rules_DUAL(etay, chiy, sz_y, order[0]);
        weights_DUAL(weights_x, order[0], sz_x[0]);                                                           
        weights_DUAL(weights_y, order[0], sz_y[0]);
        ev_basis_function( 2*sz_x[0], 0, phii, etax, chiy, "P1di" );                                            
        ev_basis_function( 2*sz_x[0], 0, phij, etax, chiy, "P1di" ); 
        int stop =0;
        for( k = 0; k <Mesh1->P1P->sbf[idx]->len && stop==0; k=k+2 )                                           
        {                                                                                                     
            assembly_P1d1(Mesh1, i, dofx, k/2, sz_x, P0_x, IE_x, etax, chix); 
            for(l = 0; l < Mesh2->P1P->sbf[idy]->len &&stop ==0; l=l+2)                                        
            {   
                if( Mesh1->P1P->sbf[idx]->num[k] != Mesh2->P1P->sbf[idy]->num[l] )                              
                {                                                                                             
                    assembly_P1d1(Mesh2, j, dofy, l/2, sz_y, P0_y, IE_y, etay, chiy);      
                    res = HP_kern2(kappa, weights_x, weights_y, phii, phij, P0_x, P0_y, IE_x, IE_y, Mesh1->PD1->sbfb[i]->sbf[dofx]->Q[k/2]->n, Mesh2->PD1->sbfb[j]->sbf[dofy]->Q[l/2]->n,  Mesh1->PD1->sbfb[i]->sbf[dofx]->Q[k/2]->c[order[0]-1],  Mesh2->PD1->sbfb[j]->sbf[dofy]->Q[l/2]->c[order[0]-1], length_arr, length_arr );

                    aux = gsl_complex_add(aux, res);                                                          
                }                                                                                             
                else                                                                                          
                {                                                                                             
                    aux = assemble_P1_HP_entriesd2(order, kappa, idx, idy, sing, Mesh1, Mesh2);                         
                    stop = 1;                                                                                 
                } 
            }                                                                                                 
        }                                                                                                     
    }                                                                                                         
    if(dofx==3 &&(dofy==0 || dofy==1 || dofy==2 ))                                                            
    {                                                                                                         
                                                                                                              
        int stop=0;
        quad_rules_PRIMAL(etax,chix,sz_x,order[0]);                                                             
        quad_rules_DUAL(etay, chiy, sz_y, order[0]);
        weights_PRIMAL(weights_x, order[0], sz_x[0]);                                                         
        weights_DUAL(weights_y, order[0], sz_y[0]);                                                             
        ev_basis_function( sz_x[0], 0, phii, etax, chix, "P1" );                                                
        ev_basis_function( 2*sz_y[0], 0, phij, etay, chiy, "P1di" );
        for( k = 0; k <6 && stop==0; k=k+2 )                                                                  
        {                                                                                                     
            assembly_P1d3(Mesh1, i, dofx, k, sz_x, P0_x, IE_x, etax, chix);                
            for(l = 0; l < Mesh2->P1P->sbf[idy]->len && stop==0; l=l+2)                                        
            {                                                                                                 
                if( Mesh1->PD1->T[i]->num[k] != Mesh2->P1P->sbf[idy]->num[l+1] && Mesh1->PD1->T[i]->num[k+1] != Mesh2->P1P->sbf[idy]->num[l] )
                {                                                                                             
                    assembly_P1d1(Mesh2, j, dofy, l/2, sz_y, P0_y, IE_y, etay, chiy);      
                    res = HP_kern3(kappa, weights_x, weights_y, phii, phij, P0_x, P0_y, IE_x, IE_y, Mesh1->PD1->sbfb[i]->sbf[dofx]->T[k/2]->n, Mesh2->PD1->sbfb[j]->sbf[dofy]->Q[l/2]->n, Mesh1->PD1->sbfb[i]->sbf[dofx]->T[k/2]->c, Mesh2->PD1->sbfb[j]->sbf[dofy]->Q[l/2]->c[order[0]-1], sz_x[0], length_arr);
                    aux = gsl_complex_add(aux, res);                                                          
                }                                                                                             
                else                                                                                          
                {                                                                                             
                    aux = assemble_P1_HP_entriesd2(order, kappa, idx, idy, sing, Mesh1, Mesh2);                         
                    stop = 1;
                 } 
            }  
            
        }                                                                                                     
                                                                                                              
    }                                                                                                         
    else if((dofx==4 || dofx==5 || dofx==6)&&(dofy==0 || dofy==1 || dofy==2))                                 
    {                                                                                                         
                                                                                                              
        int stop = 0;
        quad_rules_PRIMAL(etax,chix,sz_x,order[0]);                                                             
        quad_rules_DUAL(etay, chiy, sz_y, order[0]);
        weights_PRIMAL(weights_x, order[0], sz_x[0]);                                                         
        weights_DUAL(weights_y, order[0], sz_y[0]);                                                              
        ev_basis_function( sz_x[0], 0, phii, etax, chix, "P1" );                                                
        ev_basis_function( 2*sz_y[0], 0, phij, etay, chiy, "P1di" ); 

        for( k = 0; k <Mesh1->P1P->sbf[idx]->len && stop==0; k=k+1 )                                           
        {                                                                                                     
            assembly_P1d2(Mesh1, idx, k, sz_x, P0_x, IE_x, etax, chix);              
            for(l = 0; l < Mesh2->P1P->sbf[idy]->len && stop==0; l=l+2)                                        
            {                                                                    
                if( Mesh1->P1P->sbf[idx]->num[k] != Mesh2->P1P->sbf[idy]->num[l] && Mesh1->P1P->sbf[idx]->num[k] != Mesh2->P1P->sbf[idy]->num[l+1] )
                {                                                                                             
                    assembly_P1d1(Mesh2, j, dofy, l/2, sz_y, P0_y, IE_y, etay, chiy);      
                    res = HP_kern3(kappa, weights_x, weights_y, phii, phij, P0_x, P0_y, IE_x, IE_y, Mesh1->P1P->sbf[idx]->T[k]->n, Mesh2->PD1->sbfb[j]->sbf[dofy]->Q[l/2]->n, Mesh1->P1P->sbf[idx]->T[k]->c, Mesh2->PD1->sbfb[j]->sbf[dofy]->Q[l/2]->c[order[0]-1], sz_x[0], length_arr );
                    aux = gsl_complex_add(aux, res);                                                          
                }                                                                                             
                else if(Mesh1->P1P->sbf[idx]->num[k] == Mesh2->P1P->sbf[idy]->num[l] || Mesh1->P1P->sbf[idx]->num[k] == Mesh2->P1P->sbf[idy]->num[l+1] )
                {                                                                                             
                    aux = assemble_P1_HP_entriesd2(order, kappa, idx, idy, sing, Mesh1, Mesh2);                         
                    stop = 1;                                                                                 
                }                                                    
            }                                                                                                 
        }                                                                                                     
                                                                                                              
    }                                                                                                         
    else if((dofx==0 || dofx==1 || dofx==2) &&dofy==3)                                                        
    {                                                                                                         
        int stop = 0;
        quad_rules_DUAL(etax, chix, sz_x, order[0]);                                                            
        quad_rules_PRIMAL(etay,chiy,sz_y,order[0]);  
        weights_DUAL(weights_x, order[0], sz_x[0]);                                                           
        weights_PRIMAL(weights_y, order[0], sz_y[0]);                                                             
        ev_basis_function( 2*sz_x[0], 0, phii, etax, chix, "P1di" );                                            
        ev_basis_function( sz_y[0], 0, phij, etay, chiy, "P1" ); 

        for( k = 0; k <Mesh1->P1P->sbf[idx]->len && stop==0; k=k+2 )                                           
        {                                                                                                     
            assembly_P1d1(Mesh1, i, dofx, k/2, sz_x, P0_x, IE_x, etax, chix);              
            for(l = 0; l < 6 && stop==0; l=l+2)                                                               
            {                                                                                                 
                if( Mesh2->PD1->T[j]->num[l] != Mesh1->P1P->sbf[idx]->num[k+1] && Mesh2->PD1->T[j]->num[l+1] != Mesh1->P1P->sbf[idx]->num[k] )
                {                                                                                             
                    assembly_P1d3(Mesh2, j, dofy, l, sz_y, P0_y, IE_y, etay, chiy);        
                    res = HP_kern4(kappa, weights_x, weights_y, phii, phij, P0_x, P0_y, IE_x, IE_y, Mesh1->PD1->sbfb[i]->sbf[dofx]->Q[k/2]->n, Mesh2->PD1->sbfb[j]->sbf[dofy]->T[l/2]->n, Mesh1->PD1->sbfb[i]->sbf[dofx]->Q[k/2]->c[order[0]-1], Mesh2->PD1->sbfb[j]->sbf[dofy]->T[l/2]->c, length_arr, sz_y[0]);  
                    aux = gsl_complex_add(aux, res);                                                          
                }                                                                                             
                 else                                                                                         
                {                                                                                             
                    aux = assemble_P1_HP_entriesd2(order, kappa, idx, idy, sing, Mesh1, Mesh2);                         
                    stop = 1;                                                                                 
                                                                                                              
                }                                                                                             
            }                                                                                                 
        }                                                                                                     
                                                                                                              
    }
    else if(dofx==3 &&dofy==3)                                                                                
    {            
         
        quad_rules_PRIMAL(etax,chix,sz_x,order[0]);                                                             
        quad_rules_PRIMAL(etay,chiy,sz_y,order[0]);    
        weights_PRIMAL(weights_x, order[0], sz_x[0]);                                                         
        weights_PRIMAL(weights_y, order[0], sz_y[0]);                                                            
        ev_basis_function( sz_x[0], 0, phii, etax, chix, "P1" );                                                
        ev_basis_function( sz_y[0], 0, phij, etay, chiy, "P1" );  
        if(i!=j)                                                                                              
        {                                                                                                      
            for( k = 0; k <6; k=k+2 )                                                                         
            {                                                                                                 
                assembly_P1d3(Mesh1, i, dofx, k, sz_x, P0_x, IE_x, etax, chix);            
                for(l = 0; l < 6; l=l+2)                                                                      
                {                                                                                             
                    assembly_P1d3(Mesh2, j, dofy, l, sz_y, P0_y, IE_y, etay, chiy);        
                    res = HP_kern(kappa, weights_x, weights_y, phii, phij, P0_x, P0_y, IE_x, IE_y, Mesh1->PD1->sbfb[i]->sbf[dofx]->T[k/2]->n, Mesh2->PD1->sbfb[j]->sbf[dofy]->T[l/2]->n, Mesh1->PD1->sbfb[i]->sbf[dofx]->T[k/2]->c, Mesh2->PD1->sbfb[j]->sbf[dofy]->T[l/2]->c, sz_x[0], sz_y[0] );
                    aux = gsl_complex_add(aux, res);                                                          
                }                                                                                             
            }                                                                                                 
        }                                                                                                     
        else                                                                                                  
        {  
            res = assemble_P1_HP_entriesd2(order, kappa, idx, idy, sing, Mesh1, Mesh2);                                 
            aux = gsl_complex_add(aux, res);                                                                  
        }                                                                                                     
    }                                                                                                         
    else if((dofx==4 || dofx==5 || dofx==6)&&dofy==3)                                                         
    {                                                                                                         
        int stop = 0; 
    
        quad_rules_PRIMAL(etax,chix,sz_x,order[0]);                                                             
        quad_rules_PRIMAL(etay,chiy,sz_y,order[0]);     
        weights_PRIMAL(weights_x, order[0], sz_x[0]);                                                         
        weights_PRIMAL(weights_y, order[0], sz_y[0]);                                                        
        ev_basis_function( sz_x[0], 0, phii, etax, chix, "P1" );                                                
        ev_basis_function( sz_y[0], 0, phij, etay, chiy, "P1" );                                               
        for( k = 0; k < Mesh1->P1P->sbf[idx]->len && stop==0; k=k+1 )                                            
        {                                                                                                     
            assembly_P1d2(Mesh1, idx, k, sz_x, P0_x, IE_x, etax, chix);              
                                                                                                              
            for(l = 0; l < 6 && stop==0; l=l+2)                                                               
            {                                                                                                 
                if( Mesh1->P1P->sbf[idx]->num[k] != Mesh2->PD1->T[j]->num[l+1] && Mesh1->P1P->sbf[idx]->num[k] != Mesh2->PD1->T[j]->num[l] )
                {                                                                                             
                    assembly_P1d3(Mesh2, j, dofy, l, sz_y, P0_y, IE_y, etay, chiy);        
                    res = HP_kern(kappa, weights_x, weights_y, phii, phij, P0_x, P0_y, IE_x, IE_y, Mesh1->P1P->sbf[idx]->T[k]->n, Mesh2->PD1->sbfb[j]->sbf[dofy]->T[l/2]->n, Mesh1->P1P->sbf[idx]->T[k]->c, Mesh2->PD1->sbfb[j]->sbf[dofy]->T[l/2]->c, sz_x[0], sz_y[0] );
                    aux = gsl_complex_add(aux, res);                                                          
                }                                                                                             
                else if( Mesh1->P1P->sbf[idx]->num[k] == Mesh2->PD1->T[j]->num[l] || Mesh1->P1P->sbf[idx]->num[k] == Mesh2->PD1->T[j]->num[l+1])
                {                                                                                             
                    aux = assemble_P1_HP_entriesd2(order, kappa, idx, idy, sing, Mesh1, Mesh2);                         
                    stop = 1;                                                                                 
                                                                                                              
                }                                                                                             
            }                                                                                                 
        }                                                                                                     
    }
    else if((dofx==0 || dofx==1 || dofx==2) &&(dofy==4 || dofy==5 || dofy==6))                                
    {                                                                                                         
        int stop = 0;
    
        quad_rules_DUAL(etax, chix, sz_x, order[0]);                                                            
        quad_rules_PRIMAL(etay, chiy, sz_y, order[0]);
        weights_DUAL(weights_x, order[0], sz_x[0]);                                                           
        weights_PRIMAL(weights_y, order[0], sz_y[0]);                                                              
        ev_basis_function( 2*sz_x[0], 0, phii, etax, chix, "P1di" );                                            
        ev_basis_function( sz_y[0], 0, phij, etay, chiy, "P1" );       
        for( k = 0; k <Mesh1->P1P->sbf[idx]->len && stop==0; k=k+2 )                                           
        {                                                                                                     
            assembly_P1d1(Mesh1, i, dofx, k/2, sz_x, P0_x, IE_x, etax, chix);                    
            for(l = 0; l < Mesh2->P1P->sbf[idy]->len && stop==0; l=l+1)                                        
            {                                                                                                 
                if( Mesh1->P1P->sbf[idx]->num[k] != Mesh2->P1P->sbf[idy]->num[l] && Mesh1->P1P->sbf[idx]->num[k+1] != Mesh2->P1P->sbf[idy]->num[l] )
                {                                                                                             
                    assembly_P1d2(Mesh2, idy, l, sz_y, P0_y, IE_y, etay, chiy);      
                    res = HP_kern4(kappa, weights_x, weights_y, phii, phij, P0_x, P0_y, IE_x, IE_y, Mesh1->PD1->sbfb[i]->sbf[dofx]->Q[k/2]->n, Mesh2->P1P->sbf[idy]->T[l]->n, Mesh1->PD1->sbfb[i]->sbf[dofx]->Q[k/2]->c[order[0]-1], Mesh2->P1P->sbf[idy]->T[l]->c, length_arr, sz_y[0] );
                    aux = gsl_complex_add(aux, res);                                                          
                }                                                                                             
                else if(Mesh1->P1P->sbf[idx]->num[k] == Mesh2->P1P->sbf[idy]->num[l] || Mesh1->P1P->sbf[idx]->num[k+1] == Mesh2->P1P->sbf[idy]->num[l] )
                {                                                                                             
                    aux = assemble_P1_HP_entriesd2(order, kappa, idx, idy, sing, Mesh1, Mesh2);                         
                    stop = 1;                                                                                 
                }                                                                                             
            }                                                                                                 
        }                                                                                                     
    }                                                                                                         
    else if(dofx==3 &&(dofy==4 || dofy==5 || dofy==6))                                                        
    {                                                                                                         
        int stop = 0; 
     
        quad_rules_PRIMAL(etax,chix,sz_x,order[0]);                                                             
        quad_rules_PRIMAL(etay,chiy,sz_y,order[0]);   
        weights_PRIMAL(weights_x, order[0], sz_x[0]);                                                         
        weights_PRIMAL(weights_y, order[0], sz_y[0]);                                                            
        ev_basis_function( sz_x[0], 0, phii, etax, chix, "P1" );                                                
        ev_basis_function( sz_y[0], 0, phij, etay, chiy, "P1" );                                                
        for( k = 0; k < 6 && stop==0; k=k+2 )                                                                 
        {                                                                                                     
            assembly_P1d3(Mesh1, i, dofx, k, sz_x, P0_x, IE_x, etax, chix);                
                                                                                                              
            for(l = 0; l < Mesh2->P1P->sbf[idy]->len && stop==0; l=l+1)                                        
            {                                                                                                 
                if( Mesh2->P1P->sbf[idy]->num[l] != Mesh1->PD1->T[i]->num[k+1] && Mesh2->P1P->sbf[idy]->num[l] != Mesh1->PD1->T[i]->num[k] )
                {                                                                                             
                    assembly_P1d2(Mesh2, idy, l, sz_y, P0_y, IE_y, etay, chiy);      
                    res = HP_kern(kappa, weights_x, weights_y, phii, phij, P0_x, P0_y, IE_x, IE_y, Mesh1->PD1->sbfb[i]->sbf[dofx]->T[k/2]->n, Mesh2->P1P->sbf[idy]->T[l]->n, Mesh1->PD1->sbfb[i]->sbf[dofx]->T[k/2]->c, Mesh2->P1P->sbf[idy]->T[l]->c, sz_x[0], sz_y[0] );
                    aux = gsl_complex_add(aux, res);                                                          
                }                                                                                             
                else if(Mesh2->P1P->sbf[idy]->num[l] == Mesh1->PD1->T[i]->num[k+1] || Mesh2->P1P->sbf[idy]->num[l] == Mesh1->PD1->T[i]->num[k])
                {                                                                                             
                    aux = assemble_P1_HP_entriesd2(order, kappa, idx, idy, sing, Mesh1, Mesh2);                         
                    stop = 1;                                                                                 
                }                                                                                             
            }                                                                                                 
        }                                                                                                     
    }                                                                                                         
    else if((dofx==4 || dofx==5 || dofx==6)&&(dofy==4 || dofy==5 || dofy==6))                                 
    {            
        
        res = assemble_P1_HP_entriesd2(order, kappa, idx, idy, sing, Mesh1, Mesh2);                                      
        aux = gsl_complex_add(aux, res);                                                                       
    }                                                                                                                                                            
    return aux;                                                                                               
}


gsl_complex assemble_P1_HP_entriesd2(int *order, double kappa, int i, int j, double *sing, pMeshes Mesh1, pMeshes Mesh2)        
{                                                                                                             
                                                                                                              
    double weights[4], weights2[4], eta[4], eta2[4], chi2[4], chi[4];                                                      
    int sz[2], sz2[2], k, l,m, length_arr;                                                                    
    quad_rules_PRIMAL(eta, chi, sz, order[0]);                                                                
    quad_rules_PRIMAL(eta2, chi2, sz2, order[1]);                                                             
    weights_PRIMAL(weights, order[0], sz[0]);  
    weights_PRIMAL(weights2, order[1], sz2[0]);
    length_arr = sz[0];                                                                                       
    double P0_x[3*length_arr], P0_y[3*length_arr], IE_x[length_arr], IE_y[length_arr];                        
    double phii[16] = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};                                       
    double phij[16] = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};                                       
    double phij2[16] = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};  
    ev_basis_function( length_arr, 0, phii, eta, chi, "P1" ); 
    ev_basis_function( length_arr, 0, phij, eta, chi, "P1" ); 
    ev_basis_function( 1, 0, phij2, eta2, chi2, "P1" ); 
    gsl_complex res;                                                                                          
    gsl_vector_complex *di;                                                                                   
    gsl_vector_complex *dj;                                                                                   
    gsl_matrix *n_x;                                                                                          
    gsl_matrix *n_y;                                                                                          
    gsl_vector *curl_x;                                                                                       
    gsl_vector *curl_y;                                                                                       
    double cres[2];                                                                                           
    gsl_complex aux = gsl_complex_rect(0.0,0.0);                                                              
    gsl_complex aux2 = gsl_complex_rect(0.0,0.0);                                                             
   for( k = 0; k <Mesh1->P1P->sbf[i]->len; k++ )                                                               
   {                                                                                                          
        integration_elements_PRIMAL( sz[0],sz[1], eta, chi, Mesh1->P1P->sbf[i]->T[k]->v, IE_x, P0_x );          
        for(l = 0; l < Mesh2->P1P->sbf[j]->len; l++)                                                           
        {                                                                                                     
           if( Mesh1->P1P->sbf[i]->num[k] != Mesh2->P1P->sbf[j]->num[l] )                                       
           {                                                                                                  
                integration_elements_PRIMAL( sz[0],sz[1], eta, chi , Mesh2->P1P->sbf[j]->T[l]->v, IE_y, P0_y );
                res = HP_kern(kappa, weights, weights, phii, phij, P0_x, P0_y, IE_x, IE_y, Mesh1->P1P->sbf[i]->T[k]->n, Mesh2->P1P->sbf[j]->T[l]->n, Mesh1->P1P->sbf[i]->T[k]->c, Mesh2->P1P->sbf[j]->T[l]->c, length_arr, length_arr );
                aux = gsl_complex_add(aux, res);                                                              
           }                                                                                                  
           else                                                                                               
           {                                                                                                  
                double resc, resn;                                                                            
                integration_elements_PRIMAL( sz2[0],sz2[1], eta2, chi2, Mesh1->P1P->sbf[i]->T[k]->v, IE_x, P0_x );
                integration_elements_PRIMAL( sz[0],sz[1], eta, chi , Mesh2->P1P->sbf[j]->T[l]->v, IE_y, P0_y );
                singular(P0_x, Mesh2->P1P->sbf[j]->T[l]->v, weights2, cres, IE_x,sz2, order[2], 0, kappa, "P1", phij2[0]);                         
                gsl_blas_ddot(Mesh1->P1P->sbf[i]->T[k]->n, Mesh2->P1P->sbf[j]->T[l]->n, &resn);                                                               
                res =  gsl_complex_rect(-resn*kappa*kappa*cres[0],-resn*kappa*kappa*cres[1]);                 
                aux2 = gsl_complex_add(aux2, res);                                                            
                singular(P0_x, Mesh2->P1P->sbf[j]->T[l]->v, weights, cres, IE_x, sz2, order[2], 1, kappa, "P0", 1.0);   
                gsl_blas_ddot(Mesh1->P1P->sbf[i]->T[k]->c, Mesh2->P1P->sbf[j]->T[l]->c, &resc);                                                         
                res = gsl_complex_rect(cres[0]*resc,cres[1]*resc);                                            
                aux2 = gsl_complex_add(aux2, res);                                                            
                                                                                                              
           }
        }                                                                                                     
     }                                                                                                        
                                                                                                              
     sing[0] = GSL_REAL(aux2);                                                                                
     sing[1] = GSL_IMAG(aux2);                                                                                
     return aux;                                                                                              
}      

