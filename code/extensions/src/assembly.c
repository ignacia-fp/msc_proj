#include "assembly.h"


void assemble_SLD(
                  // Input:
                  // V: vertices of each element
                  pMeshes Mesh,
                  // wavenumber
                  double kappa,
                  // Order of quadrature
                  int *order,
                  // Output: complex and real part
                  double *M1,
                  double *M2,
                  char *trialbf,
                  char *testbf
                  ){

    double weights[16], weightsp[4], weightsp2[4], etap[4], chip[4], eta[4], chi[4], phij[16], phii[16], di[16], dj[16], vert1[9], vert2[9];
    double P0_x[3*order[0]*order[0]], P0_y[3*order[0]*order[0]], IE_x[order[0]*order[0]], IE_y[order[0]*order[0]], cres[2];
    int sz[2],szp[2], i, j, k, l, length_arr, idx;
    quad_rules_DUAL(eta, chi, sz, order[0]);
    quad_rules_PRIMAL(etap, chip, szp, order[0]);  
    weights_DUAL(weights, order[0], sz[0]);
    weights_PRIMAL(weightsp, order[0], szp[0]);
    length_arr = order[0]*order[0];
    //int order2 = 3;


    if( strcmp(testbf, "P0d" ) == 0 && strcmp(trialbf, "P0d") == 0 )
    {
        for( i = 0; i < Mesh->PD0->ndof; i++ )
        {
            for( j = 0; j < Mesh->PD0->ndof; j++ )
            {
                idx = i*Mesh->PD0->ndof + j;
                if( i != j )
                {
                    
                    for( k = 0; k <Mesh->PD0->sbf[i]->len; k++ )
                    {
                        
                        integration_elements_DUAL( sz[0],sz[1], eta, chi,Mesh->PD0->sbf[i]->Q[k]->v, IE_x, P0_x );
                        ev_basis_function( length_arr, 0, phii, eta, chi, trialbf );
                        vec_prod( length_arr, phii, weights, IE_x, di );
                        
                        for(l = 0; l < Mesh->PD0->sbf[j]->len; l++)
                        {
                            integration_elements_DUAL( sz[0], sz[1], eta, chi, Mesh->PD0->sbf[j]->Q[l]->v, IE_y, P0_y );
                            ev_basis_function( length_arr, 0, phij, eta, chi, testbf );
                            vec_prod( length_arr, phij, weights, IE_y, dj );
                            evaluate_kernel( P0_x,P0_y, kappa, length_arr, length_arr, di, dj, cres );
                            M1[idx]+= cres[0];
                            M2[idx]+= cres[1];
                        }
                    }
                }
                else
                {   
                    for( k = 0; k <Mesh->PD0->sbf[i]->len; k++ )
                    {
                        for(l = 0; l < Mesh->PD0->sbf[j]->len; l++)
                        {
                            if(l!=k)
                            {
                                 vert1[0] = Mesh->PD0->sbf[i]->Q[k]->v[0];
                                 vert1[1] = Mesh->PD0->sbf[i]->Q[k]->v[1];
                                 vert1[2] = Mesh->PD0->sbf[i]->Q[k]->v[2];
                                 vert1[3] = Mesh->PD0->sbf[i]->Q[k]->v[6];                                   
                                 vert1[4] = Mesh->PD0->sbf[i]->Q[k]->v[7];                                   
                                 vert1[5] = Mesh->PD0->sbf[i]->Q[k]->v[8];
                                 vert1[6] = Mesh->PD0->sbf[i]->Q[k]->v[9];                                   
                                 vert1[7] = Mesh->PD0->sbf[i]->Q[k]->v[10];                                   
                                 vert1[8] = Mesh->PD0->sbf[i]->Q[k]->v[11];
                                 integration_elements_PRIMAL( szp[0],szp[1], etap, chip, vert1, IE_x, P0_x );
                                 ev_basis_function( length_arr, 0, phii, etap, chip, "P0" );
                                 vec_prod( length_arr, phii, weightsp, IE_x, di );
                                 vert2[0] = Mesh->PD0->sbf[j]->Q[l]->v[0];
                                 vert2[1] = Mesh->PD0->sbf[j]->Q[l]->v[1];
                                 vert2[2] = Mesh->PD0->sbf[j]->Q[l]->v[2];
                                 vert2[3] = Mesh->PD0->sbf[j]->Q[l]->v[6];                                   
                                 vert2[4] = Mesh->PD0->sbf[j]->Q[l]->v[7];                                   
                                 vert2[5] = Mesh->PD0->sbf[j]->Q[l]->v[8];
                                 vert2[6] = Mesh->PD0->sbf[j]->Q[l]->v[9];                                   
                                 vert2[7] = Mesh->PD0->sbf[j]->Q[l]->v[10];                                   
                                 vert2[8] = Mesh->PD0->sbf[j]->Q[l]->v[11];
                                 integration_elements_PRIMAL( szp[0],szp[1], etap, chip , vert2, IE_y, P0_y );
                                 ev_basis_function( length_arr, 0, phij, etap, chip, "P0" );
                                 vec_prod( length_arr, phij, weightsp, IE_y, dj );
                                 evaluate_kernel( P0_x,P0_y, kappa, length_arr, length_arr, di, dj, cres );
                                 M1[idx]+=cres[0];
                                 M2[idx]+=cres[1];
                                 vert2[0] = Mesh->PD0->sbf[j]->Q[l]->v[0];
                                 vert2[1] = Mesh->PD0->sbf[j]->Q[l]->v[1];
                                 vert2[2] = Mesh->PD0->sbf[j]->Q[l]->v[2];
                                 vert2[3] = Mesh->PD0->sbf[j]->Q[l]->v[3];                                   
                                 vert2[4] = Mesh->PD0->sbf[j]->Q[l]->v[4];                                   
                                 vert2[5] = Mesh->PD0->sbf[j]->Q[l]->v[5];
                                 vert2[6] = Mesh->PD0->sbf[j]->Q[l]->v[6];                                   
                                 vert2[7] = Mesh->PD0->sbf[j]->Q[l]->v[7];                                   
                                 vert2[8] = Mesh->PD0->sbf[j]->Q[l]->v[8];
                                 integration_elements_PRIMAL( szp[0],szp[1], etap, chip , vert2, IE_y, P0_y );
                                 ev_basis_function( length_arr, 0, phij, etap, chip, "P0" );
                                 vec_prod( length_arr, phij, weightsp, IE_y, dj );
                                 evaluate_kernel( P0_x,P0_y, kappa, length_arr, length_arr, di, dj, cres );
                                 M1[idx]+=cres[0];
                                 M2[idx]+=cres[1];
                                 vert1[0] = Mesh->PD0->sbf[i]->Q[k]->v[0];
                                 vert1[1] = Mesh->PD0->sbf[i]->Q[k]->v[1];
                                 vert1[2] = Mesh->PD0->sbf[i]->Q[k]->v[2];
                                 vert1[3] = Mesh->PD0->sbf[i]->Q[k]->v[3];                                   
                                 vert1[4] = Mesh->PD0->sbf[i]->Q[k]->v[4];                                   
                                 vert1[5] = Mesh->PD0->sbf[i]->Q[k]->v[5];
                                 vert1[6] = Mesh->PD0->sbf[i]->Q[k]->v[6];                                   
                                 vert1[7] = Mesh->PD0->sbf[i]->Q[k]->v[7];                                   
                                 vert1[8] = Mesh->PD0->sbf[i]->Q[k]->v[8];
                                 integration_elements_PRIMAL( szp[0],szp[1], etap, chip, vert1, IE_x, P0_x );
                                 ev_basis_function( length_arr, 0, phii, etap, chip, "P0" ); 
                                 vec_prod( length_arr, phii, weightsp, IE_x, di );
                                 vert2[0] = Mesh->PD0->sbf[j]->Q[l]->v[0];
                                 vert2[1] = Mesh->PD0->sbf[j]->Q[l]->v[1];
                                 vert2[2] = Mesh->PD0->sbf[j]->Q[l]->v[2];
                                 vert2[3] = Mesh->PD0->sbf[j]->Q[l]->v[6];                                   
                                 vert2[4] = Mesh->PD0->sbf[j]->Q[l]->v[7];                                   
                                 vert2[5] = Mesh->PD0->sbf[j]->Q[l]->v[8];                                   
                                 vert2[6] = Mesh->PD0->sbf[j]->Q[l]->v[9];                                   
                                 vert2[7] = Mesh->PD0->sbf[j]->Q[l]->v[10];                                   
                                 vert2[8] = Mesh->PD0->sbf[j]->Q[l]->v[11]; 
                                 integration_elements_PRIMAL( szp[0],szp[1], etap, chip , vert2, IE_y, P0_y );
                                 ev_basis_function( length_arr, 0, phij, etap, chip, "P0" ); 
                                 vec_prod( length_arr, phij, weightsp, IE_y, dj );
                                 evaluate_kernel( P0_x,P0_y, kappa, length_arr, length_arr, di, dj, cres );
                                 M1[idx]+=cres[0];
                                 M2[idx]+=cres[1]; 
                                 vert2[0] = Mesh->PD0->sbf[j]->Q[l]->v[0];
                                 vert2[1] = Mesh->PD0->sbf[j]->Q[l]->v[1];
                                 vert2[2] = Mesh->PD0->sbf[j]->Q[l]->v[2];
                                 vert2[3] = Mesh->PD0->sbf[j]->Q[l]->v[3];                                   
                                 vert2[4] = Mesh->PD0->sbf[j]->Q[l]->v[4];                                   
                                 vert2[5] = Mesh->PD0->sbf[j]->Q[l]->v[5];
                                 vert2[6] = Mesh->PD0->sbf[j]->Q[l]->v[6];                                   
                                 vert2[7] = Mesh->PD0->sbf[j]->Q[l]->v[7];                                   
                                 vert2[8] = Mesh->PD0->sbf[j]->Q[l]->v[8];
                                 integration_elements_PRIMAL( szp[0],szp[1], etap, chip , vert2, IE_y, P0_y );
                                 ev_basis_function( length_arr, 0, phij, etap, chip, "P0" ); 
                                 vec_prod( length_arr, phij, weightsp, IE_y, dj );
                                 evaluate_kernel( P0_x,P0_y, kappa, length_arr, length_arr, di, dj, cres );
                                 M1[idx]+=cres[0];
                                 M2[idx]+=cres[1]; 
                            }
                            else
                            {
                                vert1[0] = Mesh->PD0->sbf[i]->Q[k]->v[0];
                                vert1[1] = Mesh->PD0->sbf[i]->Q[k]->v[1];
                                vert1[2] = Mesh->PD0->sbf[i]->Q[k]->v[2];
                                vert1[3] = Mesh->PD0->sbf[i]->Q[k]->v[6];                                   
                                vert1[4] = Mesh->PD0->sbf[i]->Q[k]->v[7];                                   
                                vert1[5] = Mesh->PD0->sbf[i]->Q[k]->v[8];
                                vert1[6] = Mesh->PD0->sbf[i]->Q[k]->v[9];                                   
                                vert1[7] = Mesh->PD0->sbf[i]->Q[k]->v[10];                                   
                                vert1[8] = Mesh->PD0->sbf[i]->Q[k]->v[11];
                                vert2[0] = Mesh->PD0->sbf[i]->Q[k]->v[0];
                                vert2[1] = Mesh->PD0->sbf[i]->Q[k]->v[1];
                                vert2[2] = Mesh->PD0->sbf[i]->Q[k]->v[2];
                                vert2[3] = Mesh->PD0->sbf[i]->Q[k]->v[3];                                   
                                vert2[4] = Mesh->PD0->sbf[i]->Q[k]->v[4];                                   
                                vert2[5] = Mesh->PD0->sbf[i]->Q[k]->v[5];
                                vert2[6] = Mesh->PD0->sbf[i]->Q[k]->v[6];                                   
                                vert2[7] = Mesh->PD0->sbf[i]->Q[k]->v[7];                                   
                                vert2[8] = Mesh->PD0->sbf[i]->Q[k]->v[8];
                                integration_elements_PRIMAL( szp[0],szp[1], etap, chip, vert1, IE_x, P0_x );
                                ev_basis_function( length_arr, 0, phii, etap, chip, "P0" ); 
                                vec_prod( length_arr, phii, weightsp, IE_x, di );
                                integration_elements_PRIMAL( szp[0],szp[1], etap, chip , vert2, IE_y, P0_y );
                                ev_basis_function( length_arr, 0, phij, etap, chip, "P0" ); 
                                vec_prod( length_arr, phij, weightsp, IE_y, dj );
                                evaluate_kernel( P0_x,P0_y, kappa, szp[0], szp[0], di, dj, cres );
                                M1[idx]+=cres[0];
                                M2[idx]+=cres[1];
                                integration_elements_PRIMAL( szp[0],szp[1], etap, chip, vert2, IE_x, P0_x ); 
                                ev_basis_function( length_arr, 0, phii, etap, chip, "P0" ); 
                                vec_prod( length_arr, phii, weightsp, IE_x, di );                            
                                integration_elements_PRIMAL( szp[0],szp[1], etap, chip , vert1, IE_y, P0_y );
                                ev_basis_function( length_arr, 0, phij, etap, chip, "P0" ); 
                                vec_prod( length_arr, phij, weightsp, IE_y, dj );                            
                                evaluate_kernel( P0_x,P0_y, kappa, szp[0], szp[0], di, dj, cres );                                                
                                M1[idx]+=cres[0];                                                            
                                M2[idx]+=cres[1];
                                singular(P0_x, vert2, weightsp2, cres, IE_x, szp, order[1], 1, kappa, "P0", 1.0);
                                M1[idx]+=cres[0];
                                M2[idx]+=cres[1];
                                singular(P0_y, vert1, weightsp2, cres, IE_y, szp, order[1], 1, kappa, "P0", 1.0);
                                M1[idx]+=cres[0];
                                M2[idx]+=cres[1];
                            }
                        }
                    }
                
                }

            }
        }
    }
}


void assemble_SLP(
                  // Input:
                  // Space data
                  pMeshes Mesh,
                  // wavenumber
                  double kappa,
                  // Order of quadrature
                  int *order,
                  // Output: complex and real part
                  double *M1,
                  double *M2,
                  char *trialbf,
                  char *testbf
                  ){
    
    double weights[4], weights2[4], eta[4], eta2[4], chi2[4], chi[4], di[4], dj[4], cres[2];
    int sz[2], sz2[2], i, j, k, l, length_arr;
    quad_rules_PRIMAL(eta, chi, sz, order[0]);
    quad_rules_PRIMAL(eta2, chi2, sz2, order[1]);
    weights_PRIMAL(weights, order[0], sz[0]);
    weights_PRIMAL(weights2, order[1], sz2[0]);
    length_arr = sz[0];
    double P0_x[3*length_arr], P0_y[3*length_arr], IE_x[length_arr], IE_y[length_arr];
    double phii[16] = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};                                       
    double phij[16] = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}; 
   
    if(strcmp(testbf, "P0" ) == 0 && strcmp(trialbf, "P0") == 0)
    {
        for( i = 0; i < Mesh->P0P->nt; i++ )
        {
            integration_elements_PRIMAL( sz[0],sz[1], eta, chi, Mesh->P0P->sbf[i]->T->v, IE_x, P0_x );
            vec_prod( length_arr, phii, weights, IE_x, di );
            for( j = 0; j < Mesh->P0P->nt; j++ )
            { 
                if (i != j)
                {
                       integration_elements_PRIMAL( sz[0],sz[1], eta, chi , Mesh->P0P->sbf[j]->T->v, IE_y, P0_y );
                       vec_prod( length_arr, phij, weights, IE_y, dj );
                       evaluate_kernel( P0_x,P0_y, kappa, length_arr, length_arr, di, dj, cres );
                       M1[ i*Mesh->P0P->nt + j ]+= cres[0];
                       M2[ i*Mesh->P0P->nt + j ]+= cres[1];
                }
                else
                {
                    integration_elements_PRIMAL( sz2[0],sz2[1], eta2, chi2, Mesh->P0P->sbf[i]->T->v, IE_x, P0_x ); 
                    singular(P0_x, Mesh->P0P->sbf[j]->T->v, weights2, cres, IE_x, sz,  order[2], 1, kappa, "P0", 1.0);
                    M1[ i*Mesh->P0P->nt + j ] += cres[0];
                    M2[ i*Mesh->P0P->nt + j ] += cres[1];
                    
                }
                
            }
        }
    }
    else if( strcmp(testbf,"P1") == 0 && strcmp(trialbf,"P0") == 0 )
    {
        for( i = 0; i < Mesh->P1P->ndof; i++ )
        {
            for( j = 0; j < Mesh->P0P->nt; j++ )
            {
                for( k = 0; k < Mesh->P1P->sbf[i]->len; k++ )
                {
                    integration_elements_PRIMAL( sz[0],sz[1], eta, chi, Mesh->P1P->sbf[i]->T[k]->v, IE_x, P0_x );
                    ev_basis_function( length_arr, 0, phii, eta, chi, trialbf );
                    vec_prod( length_arr, phii, weights, IE_x, di );
                    
                    if( Mesh->P1P->sbf[i]->num[k] != j )
                    {
                        integration_elements_PRIMAL( sz[0],sz[1], eta, chi , Mesh->P0P->sbf[j]->T->v, IE_y, P0_y );
                        ev_basis_function( length_arr, 0, phij, eta, chi, testbf );
                        vec_prod( length_arr, phij, weights, IE_y, dj );
                        evaluate_kernel( P0_x,P0_y, kappa, length_arr, length_arr, di, dj, cres );
                        M1[ i*Mesh->P0P->nt + j ]+= cres[0];
                        M2[ i*Mesh->P0P->nt + j ]+= cres[1];
                    }
                }
            }
        }
    }
    else if( (strcmp(testbf,"P1") == 0 && strcmp(trialbf,"P1") == 0)||(strcmp(testbf,"P1b") == 0 && strcmp(trialbf,"P1b") == 0 ) )
    {
        for( i = 0; i < Mesh->P1P->ndof; i++ )
        {
            for( j = 0; j < Mesh->P1P->ndof; j++ )                                                
            { 
                
                for( k = 0; k < Mesh->P1P->sbf[i]->len; k++ )
                {   

                    integration_elements_PRIMAL( sz[0],sz[1], eta, chi, Mesh->P1P->sbf[i]->T[k]->v, IE_x, P0_x );
                    ev_basis_function( length_arr, 0, phii, eta, chi, trialbf );
                    vec_prod( length_arr, phii, weights, IE_x, di );
                    for( l = 0; l < Mesh->P1P->sbf[j]->len; l++ )
                    {
                        if( Mesh->P1P->sbf[i]->num[k] != Mesh->P1P->sbf[j]->num[l] )
                        {
                            integration_elements_PRIMAL( sz[0],sz[1], eta, chi , Mesh->P1P->sbf[j]->T[l]->v, IE_y, P0_y );
                            ev_basis_function( length_arr, 0, phij, eta, chi, testbf );
                            vec_prod( length_arr, phij, weights, IE_y, dj );
                            evaluate_kernel( P0_x,P0_y, kappa, length_arr, length_arr, di, dj, cres );
                            M1[ i*Mesh->P1P->ndof + j ]+= cres[0];
                            M2[ i*Mesh->P1P->ndof + j ]+= cres[1];
                        }
                        else
                        {
                            integration_elements_PRIMAL( sz2[0],sz2[1], eta2, chi2, Mesh->P1P->sbf[j]->T[l]->v, IE_x, P0_x );
                            singular(P0_x, Mesh->P1P->sbf[j]->T[l]->v, weights2, cres, IE_x, sz2, order[2], 0, kappa, "P1", 1.0);
                            M1[ i*Mesh->P1P->ndof + j ]+= cres[0];
                            M2[ i*Mesh->P1P->ndof + j ]+= cres[1];
                            
                        }

                    }
                }
            }   
        }
    }
        

}


void assemble_HP(
                 // Input:
                 // Space data
                 pMeshes Mesh,
                 // wavenumber
                 double kappa,
                 // Order of quadrature
                 int *order,
                 // Output: complex and real part
                 double *M1,
                 double *M2,
                 char *trialbf,
                 char *testbf
                 ){
    
    double weights[4], weights2[4], eta[4], eta2[4], chi2[4], chi[4];
    int sz[2], sz2[2], i, j, k, l,m, length_arr;
    quad_rules_PRIMAL(eta, chi, sz, order[0]);
    quad_rules_PRIMAL(eta2, chi2, sz2, order[1]);
    weights_PRIMAL(weights, order[0], sz[0]);
    weights_PRIMAL(weights2, order[1], sz2[0]); 
    length_arr = sz[0];
    double P0_x[3*length_arr], P0_y[3*length_arr], IE_x[length_arr], IE_y[length_arr];
    double phii[16] = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
    double phij[16] = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
    double phij2[16] = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};  
    gsl_complex res;
    double cres[2]; 
    
    if( (strcmp(testbf,"P1") == 0 && strcmp(trialbf,"P1") == 0) || (strcmp(testbf,"P1b") == 0 && strcmp(trialbf,"P1b") == 0))
    {
        ev_basis_function( length_arr, 0, phii, eta, chi, trialbf );
        ev_basis_function( length_arr, 0, phij, eta, chi, testbf );
        ev_basis_function( 1, 0, phij2, eta2, chi2, "P1" );  
        for( m = 0; m < Mesh->P1P->ndof*Mesh->P1P->ndof; m++ )
        {
            i = m/Mesh->P1P->ndof;
            j = m%Mesh->P1P->ndof;
            for( k = 0; k < Mesh->P1P->sbf[i]->len; k++ )
            {
                
                integration_elements_PRIMAL( sz[0],sz[1], eta, chi, Mesh->P1P->sbf[i]->T[k]->v, IE_x, P0_x );
                //n_x = normal(P0_x, Mesh->P1P->sbf[i]->T[k]->v);
                //ev_basis_function( length_arr, 0, phii, eta, chi, trialbf );
                //curl_x = curl_bf(Mesh->P1P->sbf[i]->T[k]->v, 0, n_x);
                for( l = 0; l < Mesh->P1P->sbf[j]->len; l++ )
                {
                    
                    if( Mesh->P1P->sbf[i]->num[k] != Mesh->P1P->sbf[j]->num[l] )
                    {
                        integration_elements_PRIMAL( sz[0],sz[1], eta, chi , Mesh->P1P->sbf[j]->T[l]->v, IE_y, P0_y );
                        //ev_basis_function( length_arr, 0, phij, eta, chi, testbf );
                        //n_y = normal(P0_y, Mesh->P1P->sbf[j]->T[l]->v);
                        //curl_y = curl_bf(Mesh->P1P->sbf[j]->T[l]->v, 0, n_y);
                        res = HP_kern(kappa, weights, weights, phii, phij, P0_x, P0_y, IE_x, IE_y, Mesh->P1P->sbf[i]->T[k]->n, Mesh->P1P->sbf[j]->T[l]->n, Mesh->P1P->sbf[i]->T[k]->c, Mesh->P1P->sbf[j]->T[l]->c, length_arr, length_arr );
                        M1[ i*Mesh->P1P->ndof + j ]+= GSL_REAL(res);
                        M2[ i*Mesh->P1P->ndof + j ]+= GSL_IMAG(res);
                        
                    }
                    else
                    {
                        double resc, resn;
                        integration_elements_PRIMAL( sz2[0],sz2[1], eta2, chi2, Mesh->P1P->sbf[i]->T[k]->v, IE_x, P0_x );
                        //ev_basis_function( 1, 0, phij, eta2, chi2, "P1" );
                        integration_elements_PRIMAL( sz[0],sz[1], eta, chi , Mesh->P1P->sbf[j]->T[l]->v, IE_y, P0_y );
                        singular(P0_x, Mesh->P1P->sbf[j]->T[l]->v, weights2, cres, IE_x,sz2, order[2], 0, kappa, "P1", phij2[0]);
                        //n_y = normal(P0_y, Mesh->P1P->sbf[j]->T[l]->v);                                       
                        //curl_y = curl_bf(Mesh->P1P->sbf[j]->T[l]->v, 0, n_y);
                        gsl_blas_ddot(Mesh->P1P->sbf[i]->T[k]->n, Mesh->P1P->sbf[j]->T[l]->n, &resn);
                        M1[ i*Mesh->P1P->ndof + j ]+= (-resn*kappa*kappa*cres[0]);
                        M2[ i*Mesh->P1P->ndof + j ]+= (-resn*kappa*kappa*cres[1]);
                        singular(P0_x, Mesh->P1P->sbf[j]->T[l]->v, weights2, cres, IE_x, sz2, order[2], 1, kappa, "P0", 1.0);
                        gsl_blas_ddot(Mesh->P1P->sbf[i]->T[k]->c, Mesh->P1P->sbf[j]->T[l]->c, &resc);
                        M1[ i*Mesh->P1P->ndof + j ]+= cres[0]*(resc);                                          
                        M2[ i*Mesh->P1P->ndof + j ]+= cres[1]*(resc);
                     
                    }
                    
                }
            }
        }
    }
    else if( strcmp(testbf,"P1d") == 0 && strcmp(trialbf,"P1d") == 0)
    {
        gsl_matrix_complex* Mat;                                                                                
        Mat = gsl_matrix_complex_calloc(Mesh->P1P->ndof, Mesh->P1P->ndof);
        gsl_complex aux1, aux2, aux3;                                                                     
        aux1 = gsl_complex_rect(0.0,0.0);
        aux2 = gsl_complex_rect(0.0,0.0);
        aux3 = gsl_complex_rect(0.0,0.0); 
        double sing[2];
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
                         
        for( i = 0; i < Mesh->PD1->nt; i++ )
        {
            idx1[0] = Mesh->P1P->globaldofs[Mesh->PD1->T[i]->dofs[0]];                                    
            idx1[1] = Mesh->P1P->globaldofs[Mesh->PD1->T[i]->dofs[1]];                                    
            idx1[2] = Mesh->P1P->globaldofs[Mesh->PD1->T[i]->dofs[2]];                                    
            idx1[3] = Mesh->P1P->globaldofs[Mesh->PD1->T[i]->T[0]->dofs[2]];                                  
            idx1[4] = Mesh->P1P->globaldofs[Mesh->PD1->T[i]->T[0]->dofs[1]];                                  
            idx1[5] = Mesh->P1P->globaldofs[Mesh->PD1->T[i]->T[5]->dofs[2]];                                  
            idx1[6] = Mesh->P1P->globaldofs[Mesh->PD1->T[i]->T[3]->dofs[2]]; 
            gsl_vector_complex_set(f1, 0, gsl_complex_rect(1.0/(Mesh->PD1->lensbf[i][0]),0));             
            gsl_vector_complex_set(f1, 1, gsl_complex_rect(1.0/(Mesh->PD1->lensbf[i][1]),0));             
            gsl_vector_complex_set(f1, 2, gsl_complex_rect(1.0/(Mesh->PD1->lensbf[i][2]),0));

            for( j = 0; j < Mesh->PD1->nt; j++ )                                                                  
            {
                idx2[0] = Mesh->P1P->globaldofs[Mesh->PD1->T[j]->dofs[0]];                                     
                idx2[1] = Mesh->P1P->globaldofs[Mesh->PD1->T[j]->dofs[1]];                                     
                idx2[2] = Mesh->P1P->globaldofs[Mesh->PD1->T[j]->dofs[2]];                                     
                idx2[3] = Mesh->P1P->globaldofs[Mesh->PD1->T[j]->T[0]->dofs[2]];                              
                idx2[4] = Mesh->P1P->globaldofs[Mesh->PD1->T[j]->T[0]->dofs[1]];                              
                idx2[5] = Mesh->P1P->globaldofs[Mesh->PD1->T[j]->T[5]->dofs[2]];                              
                idx2[6] = Mesh->P1P->globaldofs[Mesh->PD1->T[j]->T[3]->dofs[2]]; 
                gsl_vector_complex_set(f2, 0, gsl_complex_rect(1.0/(Mesh->PD1->lensbf[j][0]),0));               
                gsl_vector_complex_set(f2, 1, gsl_complex_rect(1.0/(Mesh->PD1->lensbf[j][1]),0));               
                gsl_vector_complex_set(f2, 2, gsl_complex_rect(1.0/(Mesh->PD1->lensbf[j][2]),0));              

                for(k=0;k<7;k++)
                {
                    for(l=0;l<7;l++)                                                                              
                    {   
                        if(GSL_REAL(gsl_matrix_complex_get(Mat,idx1[k], idx2[l]))==0 && GSL_IMAG(gsl_matrix_complex_get(Mat,idx1[k], idx2[l]))==0)
                        {    
                        
                            sing[0] = 0;
                            sing[1] = 0;
                            aux3 = assemble_P1_HP_entries(order, kappa, idx1[k], idx2[l],sing, Mesh);
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
                M1[i*Mesh->PD1->nt+j] = GSL_REAL(aux2);
                M2[i*Mesh->PD1->nt+j] = GSL_IMAG(aux2);   
            }
        }
        gsl_matrix_complex_free(Mat);
    }
    
}


void assemble_HP2(                                                                                             
                 // Input:                                                                                    
                 // Space data                                                                                
                 pMeshes Mesh,                                                                                
                 // wavenumber                                                                                
                 double kappa,                                                                                
                 // Order of quadrature                                                                       
                 int *order,                                                                                  
                 // Output: complex and real part                                                             
                 double *M1,                                                                                  
                 double *M2,                                                                                  
                 char *trialbf,                                                                               
                 char *testbf                                                                                 
                 ){                                                                                           
                                                                                                              
    double weights[4], weights2[4], eta[4], eta2[4], chi2[4], chi[4];                                                      
    int sz[2], sz2[2], i, j, k, l,m, length_arr;                                                              
    quad_rules_PRIMAL(eta, chi, sz, order[0]);                                                                
    quad_rules_PRIMAL(eta2, chi2, sz2, order[1]);                                                             
    weights_PRIMAL(weights, order[0], sz[0]);                                
    weights_PRIMAL(weights2, order[1], sz2[0]); 
    length_arr = sz[0];                                                                                       
    double P0_x[3*length_arr], P0_y[3*length_arr], IE_x[length_arr], IE_y[length_arr];                        
    double phii[16] = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};                                       
    double phij[16] = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};                                       
    gsl_complex res;                                                                                          
    gsl_matrix *n_x;                                                                                          
    gsl_matrix *n_y;                                                                                          
    gsl_vector *curl_x;                                                                                       
    gsl_vector *curl_y;                                                                                       
    double cres[2];                                                                                           
                                                                                                              
    if( (strcmp(testbf,"P1") == 0 && strcmp(trialbf,"P1") == 0) || (strcmp(testbf,"P1b") == 0 && strcmp(trialbf,"P1b") == 0))
    {                                                                                                         
                                                                                                              
        for( m = 0; m < Mesh->P1P->ndof*Mesh->P1P->ndof; m++ )                                                
        {                                                                                                     
            i = m/Mesh->P1P->ndof;                                                                            
            j = m%Mesh->P1P->ndof;                                                                            
            for( k = 0; k < Mesh->P1P->sbf[i]->len; k++ )                                                     
            {                                                                                                 
                                                                                                              
                integration_elements_PRIMAL( sz[0],sz[1], eta, chi, Mesh->P1P->sbf[i]->T[k]->v, IE_x, P0_x ); 
                n_x = normal(P0_x, Mesh->P1P->sbf[i]->T[k]->v);                                               
                ev_basis_function( length_arr, 0, phii, eta, chi, trialbf );                                  
                curl_x = curl_bf(Mesh->P1P->sbf[i]->T[k]->v, 0, n_x);                                         
                for( l = 0; l < Mesh->P1P->sbf[j]->len; l++ )                                                 
                {                                                                                             
                                                                                                              
                    if( Mesh->P1P->sbf[i]->num[k] != Mesh->P1P->sbf[j]->num[l] )                              
                    {                                                                                         
                        integration_elements_PRIMAL( sz[0],sz[1], eta, chi , Mesh->P1P->sbf[j]->T[l]->v, IE_y, P0_y );
                        ev_basis_function( length_arr, 0, phij, eta, chi, testbf );                           
                        n_y = normal(P0_y, Mesh->P1P->sbf[j]->T[l]->v);                                       
                        curl_y = curl_bf(Mesh->P1P->sbf[j]->T[l]->v, 0, n_y);                                 
                        res = HP_kern(kappa, weights, weights, phii, phij, P0_x, P0_y, IE_x, IE_y, n_x, n_y, curl_x, curl_y, length_arr, length_arr );
                        M1[ i*Mesh->P1P->ndof + j ]+= GSL_REAL(res);                                          
                        M2[ i*Mesh->P1P->ndof + j ]+= GSL_IMAG(res);                                         
                                                                                                              
                    }                                                                                         
                    else                                                                                      
                    {                                                                                         
                        double resc, resn;                                                                    
                        integration_elements_PRIMAL( sz2[0],sz2[1], eta2, chi2, Mesh->P1P->sbf[i]->T[k]->v, IE_x, P0_x );
                        ev_basis_function( 1, 0, phij, eta2, chi2, "P1" );                                    
                        integration_elements_PRIMAL( sz[0],sz[1], eta, chi , Mesh->P1P->sbf[j]->T[l]->v, IE_y, P0_y );
                        singular(P0_x, Mesh->P1P->sbf[j]->T[l]->v, weights2, cres, IE_x,sz2, order[2], 0, kappa, "P1", phij[0]);
                        n_y = normal(P0_y, Mesh->P1P->sbf[j]->T[l]->v);                                       
                        curl_y = curl_bf(Mesh->P1P->sbf[j]->T[l]->v, 0, n_y);                                 
                        gsl_blas_ddot(n_x, n_y, &resn);                                                       
                        M1[ i*Mesh->P1P->ndof + j ]+= (-resn*kappa*kappa*cres[0]);                            
                        M2[ i*Mesh->P1P->ndof + j ]+= (-resn*kappa*kappa*cres[1]);                            
                        singular(P0_x, Mesh->P1P->sbf[j]->T[l]->v, weights2, cres, IE_x, sz2, order[2], 1, kappa, "P0", 1.0);
                        gsl_blas_ddot(curl_x, curl_y, &resc);                                                 
                        M1[ i*Mesh->P1P->ndof + j ]+= cres[0]*(resc);                                         
                        M2[ i*Mesh->P1P->ndof + j ]+= cres[1]*(resc);                                         
                                                                                                              
                    }                                                                                         
                                                                                                              
                }                                                                                             
            }                                                                                                 
        }                                                                                                     
    }                                                                                                         
    else if( strcmp(testbf,"P1d") == 0 && strcmp(trialbf,"P1d") == 0)                                         
    {                                                                                                         
        gsl_matrix_complex* Mat;                                                                              
        Mat = gsl_matrix_complex_calloc(Mesh->P1P->ndof, Mesh->P1P->ndof);                                    
        gsl_complex aux1, aux2, aux3;                                                                               
        aux1 = gsl_complex_rect(0.0,0.0);                                                                     
        aux2 = gsl_complex_rect(0.0,0.0);
        aux3 = gsl_complex_rect(0.0,0.0);
        double sing[2];                                                                          
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
    
         for( i = 0; i < Mesh->PD1->nt; i++ )
        {                                                                                                     
            /*idx1[0] = Mesh->P1P->globaldofs[Mesh->PD1->T[i]->dofs[0]];                                        
            idx1[1] = Mesh->P1P->globaldofs[Mesh->PD1->T[i]->dofs[1]];                                        
            idx1[2] = Mesh->P1P->globaldofs[Mesh->PD1->T[i]->dofs[2]];                                        
            idx1[3] = Mesh->P1P->globaldofs[Mesh->PD1->T[i]->T[0]->dofs[2]];                                  
            idx1[4] = Mesh->P1P->globaldofs[Mesh->PD1->T[i]->T[0]->dofs[1]];                                  
            idx1[5] = Mesh->P1P->globaldofs[Mesh->PD1->T[i]->T[5]->dofs[2]];                                                                                                                                            idx1[6] = Mesh->P1P->globaldofs[Mesh->PD1->T[i]->T[3]->dofs[2]];                                  
            */gsl_vector_complex_set(f1, 0, gsl_complex_rect(1.0/(Mesh->PD1->lensbf[i][0]),0));                 
            gsl_vector_complex_set(f1, 1, gsl_complex_rect(1.0/(Mesh->PD1->lensbf[i][1]),0));                 
            gsl_vector_complex_set(f1, 2, gsl_complex_rect(1.0/(Mesh->PD1->lensbf[i][2]),0));                       
            for( j = 0; j < Mesh->PD1->nt; j++ )                                                              
            {                                                                                                 
                /*idx2[0] = Mesh->P1P->globaldofs[Mesh->PD1->T[j]->dofs[0]];                                    
                idx2[1] = Mesh->P1P->globaldofs[Mesh->PD1->T[j]->dofs[1]];                                    
                idx2[2] = Mesh->P1P->globaldofs[Mesh->PD1->T[j]->dofs[2]];                                    
                idx2[3] = Mesh->P1P->globaldofs[Mesh->PD1->T[j]->T[0]->dofs[2]];                              
                idx2[4] = Mesh->P1P->globaldofs[Mesh->PD1->T[j]->T[0]->dofs[1]];                               
                idx2[5] = Mesh->P1P->globaldofs[Mesh->PD1->T[j]->T[5]->dofs[2]];                               
                idx2[6] = Mesh->P1P->globaldofs[Mesh->PD1->T[j]->T[3]->dofs[2]];                              
                */gsl_vector_complex_set(f2, 0, gsl_complex_rect(1.0/(Mesh->PD1->lensbf[j][0]),0));             
                gsl_vector_complex_set(f2, 1, gsl_complex_rect(1.0/(Mesh->PD1->lensbf[j][1]),0));             
                gsl_vector_complex_set(f2, 2, gsl_complex_rect(1.0/(Mesh->PD1->lensbf[j][2]),0));             
              
                for(k=0;k<7;k++)                                                                              
                {                                                                                             
                    for(l=0;l<7;l++)                                                                          
                    {  
                     
                        //if(GSL_REAL(gsl_matrix_complex_get(Mat,idx1[k], idx2[l]))==0 & GSL_IMAG(gsl_matrix_complex_get(Mat,idx1[k], idx2[l]))==0)
                        if(GSL_REAL(gsl_matrix_complex_get(Mat,Mesh->PD1->sbfb[i]->idx[k], Mesh->PD1->sbfb[j]->idx[l]))==0 & GSL_IMAG(gsl_matrix_complex_get(Mat,Mesh->PD1->sbfb[i]->idx[k], Mesh->PD1->sbfb[j]->idx[l]))==0)
                        {
                                sing[0] = 0;
                                sing[1] = 0;
                                aux3 = assemble_P1_HP_entriesd(order, kappa, i, j, Mesh->PD1->sbfb[i]->idx[k], Mesh->PD1->sbfb[j]->idx[l],sing, Mesh, k, l);
                                aux3 = gsl_complex_add(aux3,gsl_complex_rect(sing[0],sing[1]));
                                gsl_vector_complex_set(f3, l, aux3);
                                gsl_matrix_complex_set(Mat,Mesh->PD1->sbfb[i]->idx[k], Mesh->PD1->sbfb[j]->idx[l],aux3);
                        }                                                                                         
                        else                                              
                        {                                                                                     
                           gsl_vector_complex_set(f3, l, gsl_matrix_complex_get(Mat,Mesh->PD1->sbfb[i]->idx[k], Mesh->PD1->sbfb[j]->idx[l])); 
                        }                                                                                         
                    }                                                                                         
                    gsl_blas_zdotu(f3,f2,&aux1);                                                              
                    gsl_vector_complex_set(f4, k, aux1); 
                }                                                                                             
                                                                                                              
                gsl_blas_zdotu(f4,f1,&aux2); 
                M1[i*Mesh->PD1->nt+j] = GSL_REAL(aux2);                                                       
                M2[i*Mesh->PD1->nt+j] = GSL_IMAG(aux2);                                                       
            }                                                                                                 
        }                                                                                                     
        gsl_matrix_complex_free(Mat);                                                                         
    }                                                                                                         
                                                                                                              
} 

void assemble_maxwell_weakly(// Input:                                                                                    
                 // Space data                                                                                
                 pMeshes Mesh,                                                                                
                 // wavenumber                                                                                
                 double kappa,                                                                                
                 // Order of quadrature                                                                       
                 int *order,                                                                                  
                 // Output: complex and real part                                                             
                 double *M1,                                                                                  
                 double *M2,                                                                                  
                 char *trialbf,                                                                               
                 char *testbf )
{
    double weights[4], weights2[4], eta[4], eta2[4], chi2[4], chi[4];                                                      
    int sz[2], sz2[2], i, j, m, length_arr;                                                              
    quad_rules_PRIMAL(eta, chi, sz, order[0]);                                                                
    quad_rules_PRIMAL(eta2, chi2, sz2, order[1]);                                                             
    weights_PRIMAL(weights, order[0], sz[0]);                                                                 
    weights_PRIMAL(weights2, order[1], sz2[0]);
    length_arr = sz[0];                                                                                       
    double P0_x[3*length_arr], P0_y[3*length_arr], IE_x[length_arr], IE_y[length_arr];                        
    double *phii = (double*)malloc(3*length_arr*sizeof(double));                                      
    double *phij = (double*)malloc(3*length_arr*sizeof(double));                                                                                                                    
    double di[3*length_arr], dj[3*length_arr];                                                                                                                                       
    double cres[2]; 

    for( m = 0; m < Mesh->RWG->sz_idx; m++ )                                                                  
    { 


        i = Mesh->RWG->idx1[m]/Mesh->RWG->nedge;                                                              
        j = Mesh->RWG->idx1[m]%Mesh->RWG->nedge;                                                              
        integration_elements_PRIMAL( sz[0],sz[1], eta, chi, Mesh->RWG->sbf[i]->T1->v, IE_x, P0_x );           
        ev_basis_function_maxwell(P0_x, phii, Mesh->RWG->sbf[i]->p1, Mesh->RWG->sbf[i]->T1->a, Mesh->RWG->sbf[i]->len, Mesh->RWG->sbf[i]->T1->n, length_arr, 0, "RWG");
        vec_prod_maxwell( length_arr, weights, IE_x, phii, di ); 
        integration_elements_PRIMAL( sz[0],sz[1], eta, chi , Mesh->RWG->sbf[j]->T1->v, IE_y, P0_y );          
        ev_basis_function_maxwell(P0_y, phij, Mesh->RWG->sbf[j]->p1, Mesh->RWG->sbf[j]->T1->a, Mesh->RWG->sbf[j]->len, Mesh->RWG->sbf[j]->T1->n, length_arr, 0, "RWG");
        vec_prod_maxwell( length_arr, weights, IE_y, phij, dj );   
        evaluate_kernel_maxwell( P0_x,P0_y, kappa, length_arr, length_arr, di, dj, cres);
        M1[ i*Mesh->RWG->nedge + j ]+= cres[0];                                                               
        M2[ i*Mesh->RWG->nedge + j ]+= cres[1];                                                               
        integration_elements_PRIMAL( sz[0],sz[1], eta, chi , Mesh->RWG->sbf[j]->T2->v, IE_y, P0_y );          
        ev_basis_function_maxwell(P0_y, phij, Mesh->RWG->sbf[j]->p2, Mesh->RWG->sbf[j]->T2->a, Mesh->RWG->sbf[j]->len, Mesh->RWG->sbf[j]->T2->n, length_arr, 1, "RWG");
        vec_prod_maxwell( length_arr, weights, IE_y, phij, dj );                                          
        evaluate_kernel_maxwell( P0_x,P0_y, kappa, length_arr, length_arr, di, dj, cres ); 
        M1[ i*Mesh->RWG->nedge + j ]+= cres[0];                                                               
        M2[ i*Mesh->RWG->nedge + j ]+= cres[1];                                                               
        integration_elements_PRIMAL( sz[0],sz[1], eta, chi, Mesh->RWG->sbf[i]->T2->v, IE_x, P0_x );           
        ev_basis_function_maxwell(P0_x, phii, Mesh->RWG->sbf[i]->p2, Mesh->RWG->sbf[i]->T2->a, Mesh->RWG->sbf[i]->len, Mesh->RWG->sbf[i]->T2->n, length_arr, 1, "RWG");
        vec_prod_maxwell( length_arr, weights, IE_x, phii, di );  
        integration_elements_PRIMAL( sz[0],sz[1], eta, chi , Mesh->RWG->sbf[j]->T1->v, IE_y, P0_y );          
        ev_basis_function_maxwell(P0_y, phij, Mesh->RWG->sbf[j]->p1, Mesh->RWG->sbf[j]->T1->a, Mesh->RWG->sbf[j]->len, Mesh->RWG->sbf[j]->T1->n, length_arr, 0, "RWG");
        vec_prod_maxwell( length_arr, weights, IE_y, phij, dj );                                          
        evaluate_kernel_maxwell( P0_x,P0_y, kappa, length_arr, length_arr, di, dj, cres );
        M1[ i*Mesh->RWG->nedge + j ]+= cres[0];                                                               
        M2[ i*Mesh->RWG->nedge + j ]+= cres[1];                                                               
        integration_elements_PRIMAL( sz[0],sz[1], eta, chi , Mesh->RWG->sbf[j]->T2->v, IE_y, P0_y );          
        ev_basis_function_maxwell(P0_y, phij, Mesh->RWG->sbf[j]->p2, Mesh->RWG->sbf[j]->T2->a, Mesh->RWG->sbf[j]->len, Mesh->RWG->sbf[j]->T2->n, length_arr, 1, "RWG");
        vec_prod_maxwell( length_arr, weights, IE_y, phij, dj );                                          
        evaluate_kernel_maxwell( P0_x,P0_y, kappa, length_arr, length_arr, di, dj, cres ); 
        M1[ i*Mesh->RWG->nedge + j ]+= cres[0];                                                               
        M2[ i*Mesh->RWG->nedge + j ]+= cres[1];

    }
     for( m=0; m< Mesh->RWG->sz_sing; m++)                                                                     
    {  

        i = Mesh->RWG->idx_sing[m]/Mesh->RWG->nedge;                                                          
        j = Mesh->RWG->idx_sing[m]%Mesh->RWG->nedge;                                                          


        if(Mesh->RWG->sing1[m]==3)                                                                            
        {
            singular_maxwell_weakly2(kappa, order, Mesh->RWG->sbf[i]->T1->v, Mesh->RWG->sbf[i]->p1, Mesh->RWG->sbf[j]->p1, cres, Mesh->RWG->sbf[i]->T1->n, Mesh->RWG->sbf[i]->T1->a, Mesh->RWG->sbf[i]->len, Mesh->RWG->sbf[j]->len, 0, 0);
            M1[ i*Mesh->RWG->nedge + j ]+= cres[0];                                                           
            M2[ i*Mesh->RWG->nedge + j ]+= cres[1]; 

        }
        else if(Mesh->RWG->sing1[m]==2)
        {
            singular_maxwell_weakly_ce(kappa, order, Mesh->RWG->sbf[i]->T1->v, Mesh->RWG->sbf[j]->T1->v, Mesh->RWG->sbf[i]->p1, Mesh->RWG->sbf[j]->p1, cres,NULL, Mesh->RWG->sbf[i]->len, Mesh->RWG->sbf[j]->len, Mesh->RWG->sbf[i]->T1->a, Mesh->RWG->sbf[j]->T1->a, Mesh->RWG->sbf[i]->T1->edges, Mesh->RWG->sbf[j]->T1->edges, 0, 0);
            M1[ i*Mesh->RWG->nedge + j ]+= cres[0];                                                           
            M2[ i*Mesh->RWG->nedge + j ]+= cres[1];

        }
        else if(Mesh->RWG->sing1[m]==1)                                                                       
        {
            singular_maxwell_weakly_cv(kappa, order, Mesh->RWG->sbf[i]->T1->v, Mesh->RWG->sbf[j]->T1->v, Mesh->RWG->sbf[i]->p1, Mesh->RWG->sbf[j]->p1, cres,NULL, Mesh->RWG->sbf[i]->len, Mesh->RWG->sbf[j]->len, Mesh->RWG->sbf[i]->T1->a, Mesh->RWG->sbf[j]->T1->a, Mesh->RWG->sbf[i]->T1->dofs, Mesh->RWG->sbf[j]->T1->dofs, 0, 0);
            M1[ i*Mesh->RWG->nedge + j ]+= cres[0];                                                           
            M2[ i*Mesh->RWG->nedge + j ]+= cres[1]; 

        }
        else if(Mesh->RWG->sing1[m]==0)
        {
            integration_elements_PRIMAL( sz[0],sz[1], eta, chi, Mesh->RWG->sbf[i]->T1->v, IE_x, P0_x );           
            ev_basis_function_maxwell(P0_x, phii, Mesh->RWG->sbf[i]->p1, Mesh->RWG->sbf[i]->T1->a, Mesh->RWG->sbf[i]->len, Mesh->RWG->sbf[i]->T1->n, length_arr, 0, "RWG");
            vec_prod_maxwell( length_arr, weights, IE_x, phii, di );                                              
            integration_elements_PRIMAL( sz[0],sz[1], eta, chi , Mesh->RWG->sbf[j]->T1->v, IE_y, P0_y );          
            ev_basis_function_maxwell(P0_y, phij, Mesh->RWG->sbf[j]->p1, Mesh->RWG->sbf[j]->T1->a, Mesh->RWG->sbf[j]->len, Mesh->RWG->sbf[j]->T1->n, length_arr, 0, "RWG");
            vec_prod_maxwell( length_arr, weights, IE_y, phij, dj );                                              
            evaluate_kernel_maxwell( P0_x,P0_y, kappa, length_arr, length_arr, di, dj, cres);                     
            M1[ i*Mesh->RWG->nedge + j ]+= cres[0];                                                               
            M2[ i*Mesh->RWG->nedge + j ]+= cres[1];

        }
        if(Mesh->RWG->sing2[m]==3)
        {
            singular_maxwell_weakly2(kappa, order, Mesh->RWG->sbf[i]->T1->v, Mesh->RWG->sbf[i]->p1, Mesh->RWG->sbf[j]->p2, cres, Mesh->RWG->sbf[i]->T1->n, Mesh->RWG->sbf[i]->T1->a, Mesh->RWG->sbf[i]->len, Mesh->RWG->sbf[j]->len, 0, 1);
            M1[ i*Mesh->RWG->nedge + j ]+= cres[0];                                                           
            M2[ i*Mesh->RWG->nedge + j ]+= cres[1]; 

        }
        else if(Mesh->RWG->sing2[m]==2)
        {
            singular_maxwell_weakly_ce(kappa, order, Mesh->RWG->sbf[i]->T1->v, Mesh->RWG->sbf[j]->T2->v, Mesh->RWG->sbf[i]->p1, Mesh->RWG->sbf[j]->p2, cres,NULL, Mesh->RWG->sbf[i]->len, Mesh->RWG->sbf[j]->len, Mesh->RWG->sbf[i]->T1->a, Mesh->RWG->sbf[j]->T2->a, Mesh->RWG->sbf[i]->T1->edges, Mesh->RWG->sbf[j]->T2->edges, 0, 1);
            M1[ i*Mesh->RWG->nedge + j ]+= cres[0];                                                           
            M2[ i*Mesh->RWG->nedge + j ]+= cres[1];

        }
        else if(Mesh->RWG->sing2[m]==1)                                                                       
        {
            singular_maxwell_weakly_cv(kappa, order, Mesh->RWG->sbf[i]->T1->v, Mesh->RWG->sbf[j]->T2->v, Mesh->RWG->sbf[i]->p1, Mesh->RWG->sbf[j]->p2, cres,NULL, Mesh->RWG->sbf[i]->len, Mesh->RWG->sbf[j]->len, Mesh->RWG->sbf[i]->T1->a, Mesh->RWG->sbf[j]->T2->a, Mesh->RWG->sbf[i]->T1->dofs, Mesh->RWG->sbf[j]->T2->dofs, 0, 1);
            M1[ i*Mesh->RWG->nedge + j ]+= cres[0];                                                           
            M2[ i*Mesh->RWG->nedge + j ]+= cres[1]; 

        }
        else if(Mesh->RWG->sing2[m]==0)                                                                       
        {   
            integration_elements_PRIMAL( sz[0],sz[1], eta, chi, Mesh->RWG->sbf[i]->T1->v, IE_x, P0_x );       
            ev_basis_function_maxwell(P0_x, phii, Mesh->RWG->sbf[i]->p1, Mesh->RWG->sbf[i]->T1->a, Mesh->RWG->sbf[i]->len, Mesh->RWG->sbf[i]->T1->n, length_arr, 0, "RWG");
            vec_prod_maxwell( length_arr, weights, IE_x, phii, di ); 
            integration_elements_PRIMAL( sz[0],sz[1], eta, chi , Mesh->RWG->sbf[j]->T2->v, IE_y, P0_y );          
            ev_basis_function_maxwell(P0_y, phij, Mesh->RWG->sbf[j]->p2, Mesh->RWG->sbf[j]->T2->a, Mesh->RWG->sbf[j]->len, Mesh->RWG->sbf[j]->T2->n, length_arr, 1, "RWG");
            vec_prod_maxwell( length_arr, weights, IE_y, phij, dj );                                              
            evaluate_kernel_maxwell( P0_x,P0_y, kappa, length_arr, length_arr, di, dj, cres );                    
            M1[ i*Mesh->RWG->nedge + j ]+= cres[0];                                                               
            M2[ i*Mesh->RWG->nedge + j ]+= cres[1];

        } 
        if(Mesh->RWG->sing3[m]==3)                                                                            
        {
            singular_maxwell_weakly2(kappa, order, Mesh->RWG->sbf[i]->T2->v, Mesh->RWG->sbf[i]->p2, Mesh->RWG->sbf[j]->p1, cres, Mesh->RWG->sbf[i]->T1->n, Mesh->RWG->sbf[i]->T2->a, Mesh->RWG->sbf[i]->len,Mesh->RWG->sbf[j]->len,  1, 0);
            M1[ i*Mesh->RWG->nedge + j ]+= cres[0];                                                           
            M2[ i*Mesh->RWG->nedge + j ]+= cres[1];

        }
        else if(Mesh->RWG->sing3[m]==2)
        {
            singular_maxwell_weakly_ce(kappa, order, Mesh->RWG->sbf[i]->T2->v, Mesh->RWG->sbf[j]->T1->v, Mesh->RWG->sbf[i]->p2, Mesh->RWG->sbf[j]->p1, cres,NULL, Mesh->RWG->sbf[i]->len, Mesh->RWG->sbf[j]->len, Mesh->RWG->sbf[i]->T2->a, Mesh->RWG->sbf[j]->T1->a, Mesh->RWG->sbf[i]->T2->edges, Mesh->RWG->sbf[j]->T1->edges, 1, 0);
            M1[ i*Mesh->RWG->nedge + j ]+= cres[0];                                                           
            M2[ i*Mesh->RWG->nedge + j ]+= cres[1];

        }
        else if(Mesh->RWG->sing3[m]==1)                                                                       
        {
            singular_maxwell_weakly_cv(kappa, order, Mesh->RWG->sbf[i]->T2->v, Mesh->RWG->sbf[j]->T1->v, Mesh->RWG->sbf[i]->p2, Mesh->RWG->sbf[j]->p1, cres,NULL, Mesh->RWG->sbf[i]->len, Mesh->RWG->sbf[j]->len, Mesh->RWG->sbf[i]->T2->a, Mesh->RWG->sbf[j]->T1->a, Mesh->RWG->sbf[i]->T2->dofs, Mesh->RWG->sbf[j]->T1->dofs, 1, 0);
            M1[ i*Mesh->RWG->nedge + j ]+= cres[0];                                                           
            M2[ i*Mesh->RWG->nedge + j ]+= cres[1];

        } 
        else if(Mesh->RWG->sing3[m]==0)                                                                       
        {   
            integration_elements_PRIMAL( sz[0],sz[1], eta, chi, Mesh->RWG->sbf[i]->T2->v, IE_x, P0_x );           
            ev_basis_function_maxwell(P0_x, phii, Mesh->RWG->sbf[i]->p2, Mesh->RWG->sbf[i]->T2->a, Mesh->RWG->sbf[i]->len, Mesh->RWG->sbf[i]->T2->n, length_arr, 1, "RWG");
            vec_prod_maxwell( length_arr, weights, IE_x, phii, di );                                              
            integration_elements_PRIMAL( sz[0],sz[1], eta, chi , Mesh->RWG->sbf[j]->T1->v, IE_y, P0_y );          
            ev_basis_function_maxwell(P0_y, phij, Mesh->RWG->sbf[j]->p1, Mesh->RWG->sbf[j]->T1->a, Mesh->RWG->sbf[j]->len, Mesh->RWG->sbf[j]->T1->n, length_arr, 0, "RWG");
            vec_prod_maxwell( length_arr, weights, IE_y, phij, dj );                                              
            evaluate_kernel_maxwell( P0_x,P0_y, kappa, length_arr, length_arr, di, dj, cres );                    
            M1[ i*Mesh->RWG->nedge + j ]+= cres[0];                                                               
            M2[ i*Mesh->RWG->nedge + j ]+= cres[1];

        }
        if(Mesh->RWG->sing4[m]==3)                                                                            
        { 
            singular_maxwell_weakly2(kappa, order, Mesh->RWG->sbf[i]->T2->v, Mesh->RWG->sbf[i]->p2, Mesh->RWG->sbf[j]->p2, cres, Mesh->RWG->sbf[i]->T1->n, Mesh->RWG->sbf[i]->T2->a, Mesh->RWG->sbf[i]->len, Mesh->RWG->sbf[j]->len, 1, 1);
            M1[ i*Mesh->RWG->nedge + j ]+= cres[0];                                                           
            M2[ i*Mesh->RWG->nedge + j ]+= cres[1]; 

        }
        else if(Mesh->RWG->sing4[m]==2)
        {
            singular_maxwell_weakly_ce(kappa, order, Mesh->RWG->sbf[i]->T2->v, Mesh->RWG->sbf[j]->T2->v, Mesh->RWG->sbf[i]->p2, Mesh->RWG->sbf[j]->p2, cres,NULL, Mesh->RWG->sbf[i]->len, Mesh->RWG->sbf[j]->len, Mesh->RWG->sbf[i]->T2->a, Mesh->RWG->sbf[j]->T2->a, Mesh->RWG->sbf[i]->T2->edges, Mesh->RWG->sbf[j]->T2->edges, 1, 1);
            M1[ i*Mesh->RWG->nedge + j ]+= cres[0];                                                           
            M2[ i*Mesh->RWG->nedge + j ]+= cres[1];

        }
        else if(Mesh->RWG->sing4[m]==1)                           
        { 
            singular_maxwell_weakly_cv(kappa, order, Mesh->RWG->sbf[i]->T2->v, Mesh->RWG->sbf[j]->T2->v, Mesh->RWG->sbf[i]->p2, Mesh->RWG->sbf[j]->p2, cres,NULL, Mesh->RWG->sbf[i]->len, Mesh->RWG->sbf[j]->len, Mesh->RWG->sbf[i]->T2->a, Mesh->RWG->sbf[j]->T2->a, Mesh->RWG->sbf[i]->T2->dofs, Mesh->RWG->sbf[j]->T2->dofs, 1, 1);
            M1[ i*Mesh->RWG->nedge + j ]+= cres[0];                                                           
            M2[ i*Mesh->RWG->nedge + j ]+= cres[1];

        }
        else if(Mesh->RWG->sing4[m]==0)                                                                       
        {   
            integration_elements_PRIMAL( sz[0],sz[1], eta, chi, Mesh->RWG->sbf[i]->T2->v, IE_x, P0_x );       
            ev_basis_function_maxwell(P0_x, phii, Mesh->RWG->sbf[i]->p2, Mesh->RWG->sbf[i]->T2->a, Mesh->RWG->sbf[i]->len, Mesh->RWG->sbf[i]->T2->n, length_arr, 1, "RWG");
            vec_prod_maxwell( length_arr, weights, IE_x, phii, di );
            integration_elements_PRIMAL( sz[0],sz[1], eta, chi , Mesh->RWG->sbf[j]->T2->v, IE_y, P0_y );          
            ev_basis_function_maxwell(P0_y, phij, Mesh->RWG->sbf[j]->p2, Mesh->RWG->sbf[j]->T2->a, Mesh->RWG->sbf[j]->len, Mesh->RWG->sbf[j]->T2->n, length_arr, 1, "RWG");
            vec_prod_maxwell( length_arr, weights, IE_y, phij, dj );                                              
            evaluate_kernel_maxwell( P0_x,P0_y, kappa, length_arr, length_arr, di, dj, cres );                    
            M1[ i*Mesh->RWG->nedge + j ]+= cres[0];                                                               
            M2[ i*Mesh->RWG->nedge + j ]+= cres[1]; 

        }

    }
}

void assemble_maxwell_T(// Input:                                                                                    
                 // Space data                                                                                
                 pSpace trial,
                 pSpace test,
                 // wavenumber                                                                                
                 double kappa,                                                                                
                 // Order of quadrature                                                                       
                 int *order,                                                                                  
                 // Output: complex and real part                                                             
                 double *M1,                                                                                  
                 double *M2, 
                 double *M3,                                                                                  
                 double *M4 )                                                                               
{                                                                                                             
    double weights[4], weights2[4],eta[4], eta2[4], chi2[4], chi[4];                                          
    int sz[2], sz2[2], i, j, m, length_arr;                                                              
    quad_rules_PRIMAL(eta, chi, sz, order[0]);                                                                
    quad_rules_PRIMAL(eta2, chi2, sz2, order[1]);                                                             
    weights_PRIMAL(weights2, order[1], sz2[0]);                                                               
    weights_PRIMAL(weights, order[0], sz[0]);                                                                 
    length_arr = sz[0];                                                                                       
    double P0_x[3*length_arr], P0_y[3*length_arr], IE_x[length_arr], IE_y[length_arr];                                                                                                                
    double di[length_arr], dj[length_arr];                      
    double di2[3*length_arr], dj2[3*length_arr];                                              
    double cres1[2], cres2[2];                                                                                           
    double divx[length_arr], divy[length_arr]; 
    double phii[length_arr], phij[length_arr];                                                                
    for( m = 0; m < trial->RWG->sz_idx; m++ )                                                                  
    {          


        i = trial->RWG->idx1[m]/trial->RWG->nedge;                                                              
        j = test->RWG->idx1[m]%test->RWG->nedge;                                                              
        integration_elements_PRIMAL( sz[0],sz[1], eta, chi, trial->RWG->sbf[i]->T1->v, IE_x, P0_x );           
        divergence(trial->RWG->sbf[i]->T1->a,trial->RWG->sbf[i]->len, trial->RWG->sbf[i]->T1->v,"RWG",0, divx, length_arr);
        vec_prod( length_arr, divx, weights, IE_x, di );
        ev_basis_function_maxwell(P0_x, phii, trial->RWG->sbf[i]->p1, trial->RWG->sbf[i]->T1->a, trial->RWG->sbf[i]->len, trial->RWG->sbf[i]->T1->n, length_arr, 0, "RWG");
        vec_prod_maxwell( length_arr, weights, IE_x, phii, di2 );                                                       
        integration_elements_PRIMAL( sz[0],sz[1], eta, chi , test->RWG->sbf[j]->T1->v, IE_y, P0_y );          
        divergence(test->RWG->sbf[j]->T1->a, test->RWG->sbf[j]->len, test->RWG->sbf[j]->T1->v,"RWG",0, divy, length_arr);
        vec_prod( length_arr, divy, weights, IE_y, dj );               
        ev_basis_function_maxwell(P0_y, phij, test->RWG->sbf[j]->p1, test->RWG->sbf[j]->T1->a, test->RWG->sbf[j]->len, test->RWG->sbf[j]->T1->n, length_arr, 0, "RWG");
        vec_prod_maxwell( length_arr, weights, IE_y, phij, dj2 );                                       
        evaluate_kernel_maxwell_T( P0_x,P0_y, kappa, length_arr, length_arr, di, dj, di2, dj2, cres1, cres2 );                            
        M1[ i*test->RWG->nedge + j ]+= cres1[0];                                                               
        M2[ i*test->RWG->nedge + j ]+= cres1[1];                                                            
        M3[ i*test->RWG->nedge + j ]+= cres2[0];                                                               
        M4[ i*test->RWG->nedge + j ]+= cres2[1];

        integration_elements_PRIMAL( sz[0],sz[1], eta, chi , test->RWG->sbf[j]->T2->v, IE_y, P0_y );          
        divergence(test->RWG->sbf[j]->T2->a, test->RWG->sbf[j]->len, test->RWG->sbf[j]->T2->v,"RWG",1, divy, length_arr);
        vec_prod( length_arr, divy, weights, IE_y, dj );            
        ev_basis_function_maxwell(P0_y, phij, test->RWG->sbf[j]->p2, test->RWG->sbf[j]->T2->a, test->RWG->sbf[j]->len, test->RWG->sbf[j]->T2->n, length_arr, 1, "RWG");
        vec_prod_maxwell( length_arr, weights, IE_y, phij, dj2 );                                           
        evaluate_kernel_maxwell_T( P0_x,P0_y, kappa, length_arr, length_arr, di, dj, di2, dj2, cres1, cres2 );                                                                                       
        M1[ i*test->RWG->nedge + j ]+= cres1[0];                                                              
        M2[ i*test->RWG->nedge + j ]+= cres1[1];                                                              
        M3[ i*test->RWG->nedge + j ]+= cres2[0];                                                              
        M4[ i*test->RWG->nedge + j ]+= cres2[1];   

        integration_elements_PRIMAL( sz[0],sz[1], eta, chi, trial->RWG->sbf[i]->T2->v, IE_x, P0_x );           
        divergence(trial->RWG->sbf[i]->T2->a, trial->RWG->sbf[i]->len, trial->RWG->sbf[i]->T2->v,"RWG",1, divx, length_arr);
        vec_prod( length_arr, divx, weights, IE_x, di );      
        ev_basis_function_maxwell(P0_x, phii, trial->RWG->sbf[i]->p2, trial->RWG->sbf[i]->T2->a, trial->RWG->sbf[i]->len, trial->RWG->sbf[i]->T2->n, length_arr, 1, "RWG");
        vec_prod_maxwell( length_arr, weights, IE_x, phii, di2 );                                                 
        integration_elements_PRIMAL( sz[0],sz[1], eta, chi , test->RWG->sbf[j]->T1->v, IE_y, P0_y );          
        divergence(test->RWG->sbf[j]->T1->a, test->RWG->sbf[j]->len, test->RWG->sbf[j]->T1->v,"RWG",0, divy, length_arr);
        vec_prod( length_arr, divy, weights, IE_y, dj );      
        ev_basis_function_maxwell(P0_y, phij, test->RWG->sbf[j]->p1, test->RWG->sbf[j]->T1->a, test->RWG->sbf[j]->len, test->RWG->sbf[j]->T1->n, length_arr, 0, "RWG");
        vec_prod_maxwell( length_arr, weights, IE_y, phij, dj2 );                                                
        evaluate_kernel_maxwell_T( P0_x,P0_y, kappa, length_arr, length_arr, di, dj, di2, dj2, cres1, cres2 );                                                                                       
        M1[ i*test->RWG->nedge + j ]+= cres1[0];                                                              
        M2[ i*test->RWG->nedge + j ]+= cres1[1];                                                              
        M3[ i*test->RWG->nedge + j ]+= cres2[0];                                                              
        M4[ i*test->RWG->nedge + j ]+= cres2[1];   

        integration_elements_PRIMAL( sz[0],sz[1], eta, chi , test->RWG->sbf[j]->T2->v, IE_y, P0_y );          
        divergence(test->RWG->sbf[j]->T2->a, test->RWG->sbf[j]->len, test->RWG->sbf[j]->T2->v,"RWG",1, divy, length_arr);
        vec_prod( length_arr, divy, weights, IE_y, dj );      
        ev_basis_function_maxwell(P0_y, phij, test->RWG->sbf[j]->p2, test->RWG->sbf[j]->T2->a, test->RWG->sbf[j]->len, test->RWG->sbf[j]->T2->n, length_arr, 1, "RWG");
        vec_prod_maxwell( length_arr, weights, IE_y, phij, dj2 );                                                
        evaluate_kernel_maxwell_T( P0_x,P0_y, kappa, length_arr, length_arr, di, dj, di2, dj2, cres1, cres2 );                            
        M1[ i*test->RWG->nedge + j ]+= cres1[0];                                                              
        M2[ i*test->RWG->nedge + j ]+= cres1[1];                                                              
        M3[ i*test->RWG->nedge + j ]+= cres2[0];                        
        M4[ i*test->RWG->nedge + j ]+= cres2[1];



    }
    for( m=0; m< trial->RWG->sz_sing; m++)                                                                     
    {                                                                                                         
        i = trial->RWG->idx_sing[m]/trial->RWG->nedge;                                                          
        j = test->RWG->idx_sing[m]%test->RWG->nedge;                                                          
                                                                                                              
                                                                                                              
        if(trial->RWG->sing1[m]==3)                                                                            
        {                                                                                                     
            singular_maxwell_hypersingular2(kappa, order, trial->RWG->sbf[i]->T1->v, trial->RWG->sbf[i]->p1, test->RWG->sbf[j]->p1, cres1, trial->RWG->sbf[i]->T1->n, trial->RWG->sbf[i]->T1->a, trial->RWG->sbf[i]->len, test->RWG->sbf[j]->len, 0, 0);
                                                                                                              
            M1[ i*test->RWG->nedge + j ]+= cres1[0];                                                          
            M2[ i*test->RWG->nedge + j ]+= cres1[1];
            singular_maxwell_weakly2(kappa, order, trial->RWG->sbf[i]->T1->v, trial->RWG->sbf[i]->p1, test->RWG->sbf[j]->p1, cres2, trial->RWG->sbf[i]->T1->n, trial->RWG->sbf[i]->T1->a, trial->RWG->sbf[i]->len, test->RWG->sbf[j]->len, 0, 0);
            M3[ i*test->RWG->nedge + j ]+= cres2[0];                                                           
            M4[ i*test->RWG->nedge + j ]+= cres2[1];                                                           
                                                                                                              
        }                                                                                                     
        else if(trial->RWG->sing1[m]==2)                                                                       
        {                                                                                                     
            singular_maxwell_hypersingular_ce(kappa, order, trial->RWG->sbf[i]->T1->v, test->RWG->sbf[j]->T1->v, trial->RWG->sbf[i]->p1, test->RWG->sbf[j]->p1, cres1,NULL, trial->RWG->sbf[i]->len, test->RWG->sbf[j]->len, trial->RWG->sbf[i]->T1->a, test->RWG->sbf[j]->T1->a, trial->RWG->sbf[i]->T1->edges, test->RWG->sbf[j]->T1->edges, 0, 0);
            M1[ i*test->RWG->nedge + j ]+= cres1[0];                                                           
            M2[ i*test->RWG->nedge + j ]+= cres1[1];
            singular_maxwell_weakly_ce(kappa, order, trial->RWG->sbf[i]->T1->v, test->RWG->sbf[j]->T1->v, trial->RWG->sbf[i]->p1, test->RWG->sbf[j]->p1, cres2,NULL, trial->RWG->sbf[i]->len, test->RWG->sbf[j]->len, trial->RWG->sbf[i]->T1->a, test->RWG->sbf[j]->T1->a, trial->RWG->sbf[i]->T1->edges, test->RWG->sbf[j]->T1->edges, 0, 0);
            M3[ i*test->RWG->nedge + j ]+= cres2[0];                                                           
            M4[ i*test->RWG->nedge + j ]+= cres2[1];                                                            
                                                                                                              
        }                                                                                                     
        else if(trial->RWG->sing1[m]==1)                                                                       
        {                                                                                                     
            singular_maxwell_hypersingular_cv(kappa, order, trial->RWG->sbf[i]->T1->v, test->RWG->sbf[j]->T1->v, trial->RWG->sbf[i]->p1, test->RWG->sbf[j]->p1, cres1,NULL, trial->RWG->sbf[i]->len, test->RWG->sbf[j]->len, trial->RWG->sbf[i]->T1->a, test->RWG->sbf[j]->T1->a, trial->RWG->sbf[i]->T1->dofs, test->RWG->sbf[j]->T1->dofs, 0, 0);
            M1[ i*test->RWG->nedge + j ]+= cres1[0];                                                           
            M2[ i*test->RWG->nedge + j ]+= cres1[1];                                                           
            singular_maxwell_weakly_cv(kappa, order, trial->RWG->sbf[i]->T1->v, test->RWG->sbf[j]->T1->v, trial->RWG->sbf[i]->p1, test->RWG->sbf[j]->p1, cres2,NULL, trial->RWG->sbf[i]->len, test->RWG->sbf[j]->len, trial->RWG->sbf[i]->T1->a, test->RWG->sbf[j]->T1->a, trial->RWG->sbf[i]->T1->dofs, test->RWG->sbf[j]->T1->dofs, 0, 0);
            M3[ i*test->RWG->nedge + j ]+= cres2[0];                                                           
            M4[ i*test->RWG->nedge + j ]+= cres2[1];                                                                                                 
        }                                                                                                     
        else if(trial->RWG->sing1[m]==0)                                                                       
        {                                                                                                     
            integration_elements_PRIMAL( sz[0],sz[1], eta, chi, trial->RWG->sbf[i]->T1->v, IE_x, P0_x );       
            divergence(trial->RWG->sbf[i]->T1->a, trial->RWG->sbf[i]->len, trial->RWG->sbf[i]->T1->v,"RWG",0, divx, length_arr);
            vec_prod( length_arr, divx, weights, IE_x, di );                                                  
            ev_basis_function_maxwell(P0_x, phii, trial->RWG->sbf[i]->p1, trial->RWG->sbf[i]->T1->a, trial->RWG->sbf[i]->len, trial->RWG->sbf[i]->T1->n, length_arr, 0, "RWG");
            vec_prod_maxwell( length_arr, weights, IE_x, phii, di2 );  
            integration_elements_PRIMAL( sz[0],sz[1], eta, chi , test->RWG->sbf[j]->T1->v, IE_y, P0_y );      
            divergence(test->RWG->sbf[j]->T1->a, test->RWG->sbf[j]->len, test->RWG->sbf[j]->T1->v,"RWG",0, divy, length_arr);
            vec_prod( length_arr, divy, weights, IE_y, dj );          
            ev_basis_function_maxwell(P0_y, phij, test->RWG->sbf[j]->p1, test->RWG->sbf[j]->T1->a, test->RWG->sbf[j]->len, test->RWG->sbf[j]->T1->n, length_arr, 0, "RWG");
            vec_prod_maxwell( length_arr, weights, IE_y, phij, dj2 );                                        
            evaluate_kernel_maxwell_T( P0_x,P0_y, kappa, length_arr, length_arr, di, dj, di2, dj2, cres1, cres2 );                        
            M1[ i*test->RWG->nedge + j ]+= cres1[0];                                                           
            M2[ i*test->RWG->nedge + j ]+= cres1[1];
            M3[ i*test->RWG->nedge + j ]+= cres2[0];                                                          
            M4[ i*test->RWG->nedge + j ]+= cres2[1];                                                           
                                                                                                              
        }                                                                                                     
        if(trial->RWG->sing2[m]==3)                                                                            
        {                                                                                                     
            singular_maxwell_hypersingular2(kappa, order, trial->RWG->sbf[i]->T1->v, trial->RWG->sbf[i]->p1, test->RWG->sbf[j]->p2, cres1, trial->RWG->sbf[i]->T1->n, trial->RWG->sbf[i]->T1->a, trial->RWG->sbf[i]->len, test->RWG->sbf[j]->len,0, 1);
                                                                                                              
             M1[ i*test->RWG->nedge + j ]+= cres1[0];                                                          
             M2[ i*test->RWG->nedge + j ]+= cres1[1];                                       
            singular_maxwell_weakly2(kappa, order, trial->RWG->sbf[i]->T1->v, trial->RWG->sbf[i]->p1, test->RWG->sbf[j]->p2, cres2, trial->RWG->sbf[i]->T1->n, trial->RWG->sbf[i]->T1->a, trial->RWG->sbf[i]->len, test->RWG->sbf[j]->len, 0, 1);
            M3[ i*test->RWG->nedge + j ]+= cres2[0];                                                           
            M4[ i*test->RWG->nedge + j ]+= cres2[1];                   
        } 
        else if(trial->RWG->sing2[m]==2)                                                                       
        {                                                                                                     
            singular_maxwell_hypersingular_ce(kappa, order, trial->RWG->sbf[i]->T1->v, test->RWG->sbf[j]->T2->v, trial->RWG->sbf[i]->p1, test->RWG->sbf[j]->p2, cres1,NULL, trial->RWG->sbf[i]->len, test->RWG->sbf[j]->len, trial->RWG->sbf[i]->T1->a, test->RWG->sbf[j]->T2->a, trial->RWG->sbf[i]->T1->edges, test->RWG->sbf[j]->T2->edges, 0, 1);
            M1[ i*test->RWG->nedge + j ]+= cres1[0];                                                           
            M2[ i*test->RWG->nedge + j ]+= cres1[1];  
            singular_maxwell_weakly_ce(kappa, order, trial->RWG->sbf[i]->T1->v, test->RWG->sbf[j]->T2->v, trial->RWG->sbf[i]->p1, test->RWG->sbf[j]->p2, cres2,NULL, trial->RWG->sbf[i]->len, test->RWG->sbf[j]->len, trial->RWG->sbf[i]->T1->a, test->RWG->sbf[j]->T2->a, trial->RWG->sbf[i]->T1->edges, test->RWG->sbf[j]->T2->edges, 0, 1);
            M3[ i*test->RWG->nedge + j ]+= cres2[0];                                                           
            M4[ i*test->RWG->nedge + j ]+= cres2[1];                                                         
        }                                                                                                     
        else if(trial->RWG->sing2[m]==1)                                                                       
        {                                                                                                     
            singular_maxwell_hypersingular_cv(kappa, order, trial->RWG->sbf[i]->T1->v, test->RWG->sbf[j]->T2->v, trial->RWG->sbf[i]->p1, test->RWG->sbf[j]->p2, cres1,NULL, trial->RWG->sbf[i]->len, test->RWG->sbf[j]->len, trial->RWG->sbf[i]->T1->a, test->RWG->sbf[j]->T2->a, trial->RWG->sbf[i]->T1->dofs, test->RWG->sbf[j]->T2->dofs, 0, 1);
            M1[ i*test->RWG->nedge + j ]+= cres1[0];                                                           
            M2[ i*test->RWG->nedge + j ]+= cres1[1];   
            singular_maxwell_weakly_cv(kappa, order, trial->RWG->sbf[i]->T1->v, test->RWG->sbf[j]->T2->v, trial->RWG->sbf[i]->p1, test->RWG->sbf[j]->p2, cres2,NULL, trial->RWG->sbf[i]->len, test->RWG->sbf[j]->len, trial->RWG->sbf[i]->T1->a, test->RWG->sbf[j]->T2->a, trial->RWG->sbf[i]->T1->dofs, test->RWG->sbf[j]->T2->dofs, 0, 1);
            M3[ i*test->RWG->nedge + j ]+= cres2[0];                                                           
            M4[ i*test->RWG->nedge + j ]+= cres2[1];                                                        
        }                                                                                                     
        else if(trial->RWG->sing2[m]==0)                                                                       
        {                                                                                                     
            integration_elements_PRIMAL( sz[0],sz[1], eta, chi, trial->RWG->sbf[i]->T1->v, IE_x, P0_x );       
            divergence(trial->RWG->sbf[i]->T1->a, trial->RWG->sbf[i]->len, trial->RWG->sbf[i]->T1->v,"RWG",0, divx, length_arr);
            vec_prod( length_arr, divx, weights, IE_x, di );
            ev_basis_function_maxwell(P0_x, phii, trial->RWG->sbf[i]->p1, trial->RWG->sbf[i]->T1->a, trial->RWG->sbf[i]->len, trial->RWG->sbf[i]->T1->n, length_arr, 0, "RWG");
            vec_prod_maxwell( length_arr, weights, IE_x, phii, di2 );                                                   
            integration_elements_PRIMAL( sz[0],sz[1], eta, chi , test->RWG->sbf[j]->T2->v, IE_y, P0_y );      
            divergence(test->RWG->sbf[j]->T2->a, test->RWG->sbf[j]->len, test->RWG->sbf[j]->T2->v,"RWG",1, divy, length_arr);
            vec_prod( length_arr, divy, weights, IE_y, dj );      
            ev_basis_function_maxwell(P0_y, phij, test->RWG->sbf[j]->p2, test->RWG->sbf[j]->T2->a, test->RWG->sbf[j]->len, test->RWG->sbf[j]->T2->n, length_arr, 1, "RWG");
            vec_prod_maxwell( length_arr, weights, IE_y, phij, dj2 ); 
            evaluate_kernel_maxwell_T( P0_x,P0_y, kappa, length_arr, length_arr, di, dj, di2, dj2, cres1, cres2 );                                           
            M1[ i*test->RWG->nedge + j ]+= cres1[0];                                                           
            M2[ i*test->RWG->nedge + j ]+= cres1[1];   
            M3[ i*test->RWG->nedge + j ]+= cres2[0];                                                          
            M4[ i*test->RWG->nedge + j ]+= cres2[1];                                                        
        }                                                                                                     
        if(trial->RWG->sing3[m]==3)                                                                            
        {                                                                                                     
            singular_maxwell_hypersingular2(kappa, order, trial->RWG->sbf[i]->T2->v, trial->RWG->sbf[i]->p2, test->RWG->sbf[j]->p1, cres1, trial->RWG->sbf[i]->T1->n, trial->RWG->sbf[i]->T2->a, trial->RWG->sbf[i]->len, test->RWG->sbf[j]->len, 1, 0);
            M1[ i*test->RWG->nedge + j ]+= cres1[0];                                                           
            M2[ i*test->RWG->nedge + j ]+= cres1[1];  
            singular_maxwell_weakly2(kappa, order, trial->RWG->sbf[i]->T2->v, trial->RWG->sbf[i]->p2, test->RWG->sbf[j]->p1, cres2, trial->RWG->sbf[i]->T1->n, trial->RWG->sbf[i]->T2->a, trial->RWG->sbf[i]->len, test->RWG->sbf[j]->len,  1, 0);
            M3[ i*test->RWG->nedge + j ]+= cres2[0];                                                           
            M4[ i*test->RWG->nedge + j ]+= cres2[1];                                                          
        }                                                                                                     
        else if(trial->RWG->sing3[m]==2)                                                                       
        {                                                                                                     
            singular_maxwell_hypersingular_ce(kappa, order, trial->RWG->sbf[i]->T2->v, test->RWG->sbf[j]->T1->v, trial->RWG->sbf[i]->p2, test->RWG->sbf[j]->p1, cres1,NULL, trial->RWG->sbf[i]->len, test->RWG->sbf[j]->len, trial->RWG->sbf[i]->T2->a, test->RWG->sbf[j]->T1->a, trial->RWG->sbf[i]->T2->edges, test->RWG->sbf[j]->T1->edges, 1, 0);
            M1[ i*test->RWG->nedge + j ]+= cres1[0];                                                           
            M2[ i*test->RWG->nedge + j ]+= cres1[1]; 
            singular_maxwell_weakly_ce(kappa, order, trial->RWG->sbf[i]->T2->v, test->RWG->sbf[j]->T1->v, trial->RWG->sbf[i]->p2, test->RWG->sbf[j]->p1, cres2,NULL, trial->RWG->sbf[i]->len, test->RWG->sbf[j]->len, trial->RWG->sbf[i]->T2->a, test->RWG->sbf[j]->T1->a, trial->RWG->sbf[i]->T2->edges, test->RWG->sbf[j]->T1->edges, 1, 0);
            M3[ i*test->RWG->nedge + j ]+= cres2[0];                                                           
            M4[ i*test->RWG->nedge + j ]+= cres2[1];                                                          
        }                                                                                                     
        else if(trial->RWG->sing3[m]==1)                                                                       
        {                                                                                                     
            singular_maxwell_hypersingular_cv(kappa, order, trial->RWG->sbf[i]->T2->v, test->RWG->sbf[j]->T1->v, trial->RWG->sbf[i]->p2, test->RWG->sbf[j]->p1, cres1,NULL, trial->RWG->sbf[i]->len, test->RWG->sbf[j]->len, trial->RWG->sbf[i]->T2->a, test->RWG->sbf[j]->T1->a, trial->RWG->sbf[i]->T2->dofs, test->RWG->sbf[j]->T1->dofs, 1, 0);
            M1[ i*test->RWG->nedge + j ]+= cres1[0];                                                           
            M2[ i*test->RWG->nedge + j ]+= cres1[1]; 
            singular_maxwell_weakly_cv(kappa, order, trial->RWG->sbf[i]->T2->v, test->RWG->sbf[j]->T1->v, trial->RWG->sbf[i]->p2, test->RWG->sbf[j]->p1, cres2,NULL, trial->RWG->sbf[i]->len, test->RWG->sbf[j]->len, trial->RWG->sbf[i]->T2->a, test->RWG->sbf[j]->T1->a, trial->RWG->sbf[i]->T2->dofs, test->RWG->sbf[j]->T1->dofs, 1, 0);
            M3[ i*test->RWG->nedge + j ]+= cres2[0];                                                           
            M4[ i*test->RWG->nedge + j ]+= cres2[1];                                                          
        }                                                                                                     
        else if(trial->RWG->sing3[m]==0)                                                                       
        {                                                                                                     
            integration_elements_PRIMAL( sz[0],sz[1], eta, chi, trial->RWG->sbf[i]->T2->v, IE_x, P0_x );       
            divergence(trial->RWG->sbf[i]->T2->a,trial->RWG->sbf[i]->len, trial->RWG->sbf[i]->T2->v,"RWG",1, divx, length_arr);
            vec_prod( length_arr, divx, weights, IE_x, di );                                                  
            ev_basis_function_maxwell(P0_x, phii, trial->RWG->sbf[i]->p2, trial->RWG->sbf[i]->T2->a, trial->RWG->sbf[i]->len, trial->RWG->sbf[i]->T2->n, length_arr, 1, "RWG");
            vec_prod_maxwell( length_arr, weights, IE_x, phii, di2 );
            integration_elements_PRIMAL( sz[0],sz[1], eta, chi , test->RWG->sbf[j]->T1->v, IE_y, P0_y );      
            divergence(test->RWG->sbf[j]->T1->a, test->RWG->sbf[j]->len, test->RWG->sbf[j]->T1->v,"RWG",0, divy, length_arr);
            vec_prod( length_arr, divy, weights, IE_y, dj );      
            ev_basis_function_maxwell(P0_y, phij, test->RWG->sbf[j]->p1, test->RWG->sbf[j]->T1->a, test->RWG->sbf[j]->len, test->RWG->sbf[j]->T1->n, length_arr, 0, "RWG");
            vec_prod_maxwell( length_arr, weights, IE_y, phij, dj2 );                                            
            evaluate_kernel_maxwell_T( P0_x,P0_y, kappa, length_arr, length_arr, di, dj, di2, dj2, cres1, cres2 );                      
            M1[ i*test->RWG->nedge + j ]+= cres1[0];                                                           
            M2[ i*test->RWG->nedge + j ]+= cres1[1];                                                           
            M3[ i*test->RWG->nedge + j ]+= cres2[0];                                                          
            M4[ i*test->RWG->nedge + j ]+= cres2[1];
        } 
        if(trial->RWG->sing4[m]==3)                                                                            
        {                                                                                                     
            singular_maxwell_hypersingular2(kappa, order, trial->RWG->sbf[i]->T2->v, trial->RWG->sbf[i]->p2, test->RWG->sbf[j]->p2, cres1, trial->RWG->sbf[i]->T1->n, trial->RWG->sbf[i]->T2->a, trial->RWG->sbf[i]->len, test->RWG->sbf[j]->len, 1, 1);
            M1[ i*test->RWG->nedge + j ]+= cres1[0];                                                           
            M2[ i*test->RWG->nedge + j ]+= cres1[1]; 
            singular_maxwell_weakly2(kappa, order, trial->RWG->sbf[i]->T2->v, trial->RWG->sbf[i]->p2, test->RWG->sbf[j]->p2, cres2, trial->RWG->sbf[i]->T1->n, trial->RWG->sbf[i]->T2->a, trial->RWG->sbf[i]->len, test->RWG->sbf[j]->len, 1, 1);
            M3[ i*test->RWG->nedge + j ]+= cres2[0];                                                           
            M4[ i*test->RWG->nedge + j ]+= cres2[1];                                                          
        }                                                                                                     
        else if(trial->RWG->sing4[m]==2)                                                                       
        {                                                                                                     
            singular_maxwell_hypersingular_ce(kappa, order, trial->RWG->sbf[i]->T2->v, test->RWG->sbf[j]->T2->v, trial->RWG->sbf[i]->p2, test->RWG->sbf[j]->p2, cres1,NULL, trial->RWG->sbf[i]->len, test->RWG->sbf[j]->len, trial->RWG->sbf[i]->T2->a, test->RWG->sbf[j]->T2->a, trial->RWG->sbf[i]->T2->edges, test->RWG->sbf[j]->T2->edges, 1, 1);
            M1[ i*test->RWG->nedge + j ]+= cres1[0];                                                           
            M2[ i*test->RWG->nedge + j ]+= cres1[1]; 
            singular_maxwell_weakly_ce(kappa, order, trial->RWG->sbf[i]->T2->v, test->RWG->sbf[j]->T2->v, trial->RWG->sbf[i]->p2, test->RWG->sbf[j]->p2, cres2,NULL, trial->RWG->sbf[i]->len, test->RWG->sbf[j]->len, trial->RWG->sbf[i]->T2->a, test->RWG->sbf[j]->T2->a, trial->RWG->sbf[i]->T2->edges, test->RWG->sbf[j]->T2->edges, 1, 1);
            M3[ i*test->RWG->nedge + j ]+= cres2[0];                                                           
            M4[ i*test->RWG->nedge + j ]+= cres2[1];                                                          
        }                                                                                                     
        else if(trial->RWG->sing4[m]==1)                                                                       
        {                                                                                                     
            singular_maxwell_hypersingular_cv(kappa, order, trial->RWG->sbf[i]->T2->v, test->RWG->sbf[j]->T2->v, trial->RWG->sbf[i]->p2, test->RWG->sbf[j]->p2, cres1,NULL, trial->RWG->sbf[i]->len, test->RWG->sbf[j]->len, trial->RWG->sbf[i]->T2->a, test->RWG->sbf[j]->T2->a, trial->RWG->sbf[i]->T2->dofs, test->RWG->sbf[j]->T2->dofs, 1, 1);
            M1[ i*test->RWG->nedge + j ]+= cres1[0];                                                           
            M2[ i*test->RWG->nedge + j ]+= cres1[1];  
            singular_maxwell_weakly_cv(kappa, order, trial->RWG->sbf[i]->T2->v, test->RWG->sbf[j]->T2->v, trial->RWG->sbf[i]->p2, test->RWG->sbf[j]->p2, cres2,NULL, trial->RWG->sbf[i]->len, test->RWG->sbf[j]->len, trial->RWG->sbf[i]->T2->a, test->RWG->sbf[j]->T2->a, trial->RWG->sbf[i]->T2->dofs, test->RWG->sbf[j]->T2->dofs, 1, 1);
            M3[ i*test->RWG->nedge + j ]+= cres2[0];                                                           
            M4[ i*test->RWG->nedge + j ]+= cres2[1];                                                         
        }                                                                                                     
        else if(trial->RWG->sing4[m]==0)                                                                       
        {                                                                                                     
            integration_elements_PRIMAL( sz[0],sz[1], eta, chi, trial->RWG->sbf[i]->T2->v, IE_x, P0_x );       
            divergence(trial->RWG->sbf[i]->T2->a, trial->RWG->sbf[i]->len, trial->RWG->sbf[i]->T2->v,"RWG",1, divx, length_arr);
            vec_prod( length_arr, divx, weights, IE_x, di );
            ev_basis_function_maxwell(P0_x, phii, trial->RWG->sbf[i]->p2, trial->RWG->sbf[i]->T2->a, trial->RWG->sbf[i]->len, trial->RWG->sbf[i]->T2->n, length_arr, 1, "RWG");
            vec_prod_maxwell( length_arr, weights, IE_x, phii, di2 );                                                    
            integration_elements_PRIMAL( sz[0],sz[1], eta, chi , test->RWG->sbf[j]->T2->v, IE_y, P0_y );      
            divergence(test->RWG->sbf[j]->T2->a,test->RWG->sbf[j]->len, test->RWG->sbf[j]->T2->v,"RWG",1, divy, length_arr);
            vec_prod( length_arr, divy, weights, IE_y, dj );      
            ev_basis_function_maxwell(P0_y, phij, test->RWG->sbf[j]->p2, test->RWG->sbf[j]->T2->a, test->RWG->sbf[j]->len, test->RWG->sbf[j]->T2->n, length_arr, 1, "RWG");
            vec_prod_maxwell( length_arr, weights, IE_y, phij, dj2 );                                             
            evaluate_kernel_maxwell_T( P0_x,P0_y, kappa, length_arr, length_arr, di, dj, di2, dj2, cres1, cres2 );                        
            M1[ i*test->RWG->nedge + j ]+= cres1[0];                                                           
            M2[ i*test->RWG->nedge + j ]+= cres1[1];
            M3[ i*test->RWG->nedge + j ]+= cres2[0];                                                          
            M4[ i*test->RWG->nedge + j ]+= cres2[1];                                                            
        }                                                                                                     
                                                                                                              
    }
}

void assemble_maxwell_hypersingular(// Input:                                                                                    
                 // Space data                                                                                
                 pSpace trial,
                 pSpace test,
                 // wavenumber                                                                                
                 double kappa,                                                                                
                 // Order of quadrature                                                                       
                 int *order,                                                                                  
                 // Output: complex and real part                                                             
                 double *M1,                                                                                  
                 double *M2 )
{                                                                                                             
    double weights[4], weights2[4],eta[4], eta2[4], chi2[4], chi[4];                                                      
    int sz[2], sz2[2], i, j, m, length_arr;                                                              
    quad_rules_PRIMAL(eta, chi, sz, order[0]);                                                                
    quad_rules_PRIMAL(eta2, chi2, sz2, order[1]); 
    weights_PRIMAL(weights2, order[1], sz2[0]);
    weights_PRIMAL(weights, order[0], sz[0]);                                                                 
    length_arr = sz[0];                                                                                       
    double P0_x[3*length_arr], P0_y[3*length_arr], IE_x[length_arr], IE_y[length_arr];                                                                                                         
    double di[length_arr], dj[length_arr];                                                                
    double cres[2];                                                                                           
    double divx[length_arr], divy[length_arr]; 
    for( m = 0; m < trial->RWG->sz_idx; m++ )
    {


        i = trial->RWG->idx1[m]/trial->RWG->nedge;                                                                               
        j = test->RWG->idx1[m]%test->RWG->nedge;
        integration_elements_PRIMAL( sz[0],sz[1], eta, chi, trial->RWG->sbf[i]->T1->v, IE_x, P0_x );           
        divergence(trial->RWG->sbf[i]->T1->a,trial->RWG->sbf[i]->len, trial->RWG->sbf[i]->T1->v,"RWG",0, divx, length_arr);
        vec_prod( length_arr, divx, weights, IE_x, di );
        integration_elements_PRIMAL( sz[0],sz[1], eta, chi , test->RWG->sbf[j]->T1->v, IE_y, P0_y );      
        divergence(test->RWG->sbf[j]->T1->a,test->RWG->sbf[j]->len, test->RWG->sbf[j]->T1->v,"RWG",0, divy, length_arr);
        vec_prod( length_arr, divy, weights, IE_y, dj );                                                  
        evaluate_kernel( P0_x,P0_y, kappa, length_arr, length_arr, di, dj, cres );
        M1[ i*test->RWG->nedge + j ]+= cres[0];
        M2[ i*test->RWG->nedge + j ]+= cres[1];
        integration_elements_PRIMAL( sz[0],sz[1], eta, chi , test->RWG->sbf[j]->T2->v, IE_y, P0_y );      
        divergence(test->RWG->sbf[j]->T2->a,test->RWG->sbf[j]->len, test->RWG->sbf[j]->T2->v,"RWG",1, divy, length_arr);
        vec_prod( length_arr, divy, weights, IE_y, dj );                                                  
        evaluate_kernel( P0_x,P0_y, kappa, length_arr, length_arr, di, dj, cres );
        M1[ i*test->RWG->nedge + j ]+= cres[0];                                                               
        M2[ i*test->RWG->nedge + j ]+= cres[1];  
        integration_elements_PRIMAL( sz[0],sz[1], eta, chi, trial->RWG->sbf[i]->T2->v, IE_x, P0_x );           
        divergence(trial->RWG->sbf[i]->T2->a,trial->RWG->sbf[i]->len, trial->RWG->sbf[i]->T2->v,"RWG",1, divx, length_arr);
        vec_prod( length_arr, divx, weights, IE_x, di );
        integration_elements_PRIMAL( sz[0],sz[1], eta, chi ,test->RWG->sbf[j]->T1->v, IE_y, P0_y );      
        divergence(test->RWG->sbf[j]->T1->a, test->RWG->sbf[j]->len, test->RWG->sbf[j]->T1->v,"RWG",0, divy, length_arr);
        vec_prod( length_arr, divy, weights, IE_y, dj );                                                  
        evaluate_kernel( P0_x,P0_y, kappa, length_arr, length_arr, di, dj, cres );
        M1[ i*test->RWG->nedge + j ]+= cres[0];                                                               
        M2[ i*test->RWG->nedge + j ]+= cres[1];                          
        integration_elements_PRIMAL( sz[0],sz[1], eta, chi , test->RWG->sbf[j]->T2->v, IE_y, P0_y );      
        divergence(test->RWG->sbf[j]->T2->a, test->RWG->sbf[j]->len, test->RWG->sbf[j]->T2->v,"RWG",1, divy, length_arr);
        vec_prod( length_arr, divy, weights, IE_y, dj );                                                  
        evaluate_kernel( P0_x,P0_y, kappa, length_arr, length_arr, di, dj, cres ); 
        M1[ i*test->RWG->nedge + j ]+= cres[0];                                                               
        M2[ i*test->RWG->nedge + j ]+= cres[1]; 

    }
    for( m=0; m< trial->RWG->sz_sing; m++)
    {
        i = trial->RWG->idx_sing[m]/trial->RWG->nedge;                                                              
        j = test->RWG->idx_sing[m]%test->RWG->nedge;


        if(trial->RWG->sing1[m]==3)                                                                                 
        {    
            singular_maxwell_hypersingular2(kappa, order, trial->RWG->sbf[i]->T1->v, trial->RWG->sbf[i]->p1, test->RWG->sbf[j]->p1, cres, trial->RWG->sbf[i]->T1->n, trial->RWG->sbf[i]->T1->a, trial->RWG->sbf[i]->len, test->RWG->sbf[j]->len, 0, 0); 
             M1[ i*test->RWG->nedge + j ]+= cres[0];                                                          
             M2[ i*test->RWG->nedge + j ]+= cres[1]; 

        }   
        else if(trial->RWG->sing1[m]==2)                                                                            
        {                                                                                                     
            singular_maxwell_hypersingular_ce(kappa, order, trial->RWG->sbf[i]->T1->v, test->RWG->sbf[j]->T1->v, trial->RWG->sbf[i]->p1, test->RWG->sbf[j]->p1, cres,NULL, trial->RWG->sbf[i]->len, test->RWG->sbf[j]->len, trial->RWG->sbf[i]->T1->a, test->RWG->sbf[j]->T1->a, trial->RWG->sbf[i]->T1->edges, test->RWG->sbf[j]->T1->edges, 0, 0);            
            M1[ i*test->RWG->nedge + j ]+= cres[0];                                                           
            M2[ i*test->RWG->nedge + j ]+= cres[1];

        }
        else if(trial->RWG->sing1[m]==1)
        {
            singular_maxwell_hypersingular_cv(kappa, order, trial->RWG->sbf[i]->T1->v, test->RWG->sbf[j]->T1->v, trial->RWG->sbf[i]->p1, test->RWG->sbf[j]->p1, cres,NULL, trial->RWG->sbf[i]->len, test->RWG->sbf[j]->len, trial->RWG->sbf[i]->T1->a, test->RWG->sbf[j]->T1->a, trial->RWG->sbf[i]->T1->dofs, test->RWG->sbf[j]->T1->dofs, 0, 0);
            M1[ i*test->RWG->nedge + j ]+= cres[0];                                                           
            M2[ i*test->RWG->nedge + j ]+= cres[1];

        }
        else if(trial->RWG->sing1[m]==0) 
        {
            integration_elements_PRIMAL( sz[0],sz[1], eta, chi, trial->RWG->sbf[i]->T1->v, IE_x, P0_x );           
            divergence(trial->RWG->sbf[i]->T1->a, trial->RWG->sbf[i]->len, trial->RWG->sbf[i]->T1->v,"RWG",0, divx, length_arr);
            vec_prod( length_arr, divx, weights, IE_x, di );
            integration_elements_PRIMAL( sz[0],sz[1], eta, chi , test->RWG->sbf[j]->T1->v, IE_y, P0_y );          
            divergence(test->RWG->sbf[j]->T1->a,test->RWG->sbf[j]->len, test->RWG->sbf[j]->T1->v,"RWG",0, divy, length_arr);
            vec_prod( length_arr, divy, weights, IE_y, dj );                                                      
            evaluate_kernel( P0_x,P0_y, kappa, length_arr, length_arr, di, dj, cres );                            
            M1[ i*test->RWG->nedge + j ]+= cres[0];                                                               
            M2[ i*test->RWG->nedge + j ]+= cres[1]; 

        }
        if(trial->RWG->sing2[m]==3)                                                                                 
        {                                                                                                     
            singular_maxwell_hypersingular2(kappa, order, trial->RWG->sbf[i]->T1->v, trial->RWG->sbf[i]->p1, test->RWG->sbf[j]->p2, cres, trial->RWG->sbf[i]->T1->n, trial->RWG->sbf[i]->T1->a, trial->RWG->sbf[i]->len, test->RWG->sbf[j]->len,0, 1);

             M1[ i*test->RWG->nedge + j ]+= cres[0];                                                          
             M2[ i*test->RWG->nedge + j ]+= cres[1];
        }                                                                                                     
        else if(trial->RWG->sing2[m]==2)                                                                            
        {                                                                                                     
            singular_maxwell_hypersingular_ce(kappa, order, trial->RWG->sbf[i]->T1->v, test->RWG->sbf[j]->T2->v, trial->RWG->sbf[i]->p1, test->RWG->sbf[j]->p2, cres,NULL, trial->RWG->sbf[i]->len, test->RWG->sbf[j]->len, trial->RWG->sbf[i]->T1->a, test->RWG->sbf[j]->T2->a, trial->RWG->sbf[i]->T1->edges, test->RWG->sbf[j]->T2->edges, 0, 1);              
            M1[ i*test->RWG->nedge + j ]+= cres[0];                                                         
            M2[ i*test->RWG->nedge + j ]+= cres[1];  
        } 
        else if(trial->RWG->sing2[m]==1)                                                                       
        {                                                                                                     
            singular_maxwell_hypersingular_cv(kappa, order, trial->RWG->sbf[i]->T1->v, test->RWG->sbf[j]->T2->v, trial->RWG->sbf[i]->p1, test->RWG->sbf[j]->p2, cres,NULL, trial->RWG->sbf[i]->len, test->RWG->sbf[j]->len, trial->RWG->sbf[i]->T1->a, test->RWG->sbf[j]->T2->a, trial->RWG->sbf[i]->T1->dofs, test->RWG->sbf[j]->T2->dofs, 0, 1);
            M1[ i*test->RWG->nedge + j ]+= cres[0];                                                           
            M2[ i*test->RWG->nedge + j ]+= cres[1]; 
        }
        else if(trial->RWG->sing2[m]==0) 
        {
            integration_elements_PRIMAL( sz[0],sz[1], eta, chi, trial->RWG->sbf[i]->T1->v, IE_x, P0_x );       
            divergence(trial->RWG->sbf[i]->T1->a,trial->RWG->sbf[i]->len, trial->RWG->sbf[i]->T1->v,"RWG",0, divx, length_arr);
            vec_prod( length_arr, divx, weights, IE_x, di );
            integration_elements_PRIMAL( sz[0],sz[1], eta, chi , test->RWG->sbf[j]->T2->v, IE_y, P0_y );          
            divergence(test->RWG->sbf[j]->T2->a,test->RWG->sbf[j]->len, test->RWG->sbf[j]->T2->v,"RWG",1, divy, length_arr);
            vec_prod( length_arr, divy, weights, IE_y, dj );                                                      
            evaluate_kernel( P0_x,P0_y, kappa, length_arr, length_arr, di, dj, cres );                            
            M1[ i*test->RWG->nedge + j ]+= cres[0];                                                               
            M2[ i*test->RWG->nedge + j ]+= cres[1];
        }
        if(trial->RWG->sing3[m]==3)                                                                                 
        {                                                                       
            singular_maxwell_hypersingular2(kappa, order, trial->RWG->sbf[i]->T2->v, trial->RWG->sbf[i]->p2, test->RWG->sbf[j]->p1, cres, trial->RWG->sbf[i]->T1->n, trial->RWG->sbf[i]->T2->a, trial->RWG->sbf[i]->len, test->RWG->sbf[j]->len, 1, 0);
            M1[ i*test->RWG->nedge + j ]+= cres[0];                                                           
            M2[ i*test->RWG->nedge + j ]+= cres[1];  
        }                                                                                                     
        else if(trial->RWG->sing3[m]==2)                                                                            
        {                                                                                                     
            singular_maxwell_hypersingular_ce(kappa, order, trial->RWG->sbf[i]->T2->v, test->RWG->sbf[j]->T1->v, trial->RWG->sbf[i]->p2, test->RWG->sbf[j]->p1, cres,NULL, trial->RWG->sbf[i]->len, test->RWG->sbf[j]->len, trial->RWG->sbf[i]->T2->a, test->RWG->sbf[j]->T1->a, trial->RWG->sbf[i]->T2->edges, test->RWG->sbf[j]->T1->edges, 1, 0);              
            M1[ i*test->RWG->nedge + j ]+= cres[0];                                                         
            M2[ i*test->RWG->nedge + j ]+= cres[1];
        }
        else if(trial->RWG->sing3[m]==1)                                                                       
        {                                                                                                     
            singular_maxwell_hypersingular_cv(kappa, order, trial->RWG->sbf[i]->T2->v, test->RWG->sbf[j]->T1->v, trial->RWG->sbf[i]->p2, test->RWG->sbf[j]->p1, cres,NULL, trial->RWG->sbf[i]->len, test->RWG->sbf[j]->len, trial->RWG->sbf[i]->T2->a, test->RWG->sbf[j]->T1->a, trial->RWG->sbf[i]->T2->dofs, test->RWG->sbf[j]->T1->dofs, 1, 0);
            M1[ i*test->RWG->nedge + j ]+= cres[0];                                                           
            M2[ i*test->RWG->nedge + j ]+= cres[1]; 
        }
        else if(trial->RWG->sing3[m]==0) 
        {
            integration_elements_PRIMAL( sz[0],sz[1], eta, chi, trial->RWG->sbf[i]->T2->v, IE_x, P0_x );           
            divergence(trial->RWG->sbf[i]->T2->a, trial->RWG->sbf[i]->len, trial->RWG->sbf[i]->T2->v,"RWG",1, divx, length_arr);
            vec_prod( length_arr, divx, weights, IE_x, di );                                                      
            integration_elements_PRIMAL( sz[0],sz[1], eta, chi , test->RWG->sbf[j]->T1->v, IE_y, P0_y );          
            divergence(test->RWG->sbf[j]->T1->a, test->RWG->sbf[j]->len, test->RWG->sbf[j]->T1->v,"RWG",0, divy, length_arr);
            vec_prod( length_arr, divy, weights, IE_y, dj );                                                      
            evaluate_kernel( P0_x,P0_y, kappa, length_arr, length_arr, di, dj, cres );                            
            M1[ i*test->RWG->nedge + j ]+= cres[0];                                                               
            M2[ i*test->RWG->nedge + j ]+= cres[1];
        }
        if(trial->RWG->sing4[m]==3)                                                                                 
        {                                                                                                     
            singular_maxwell_hypersingular2(kappa, order, trial->RWG->sbf[i]->T2->v, trial->RWG->sbf[i]->p2, test->RWG->sbf[j]->p2, cres, trial->RWG->sbf[i]->T1->n, trial->RWG->sbf[i]->T2->a, trial->RWG->sbf[i]->len, test->RWG->sbf[j]->len, 1, 1);
            M1[ i*test->RWG->nedge + j ]+= cres[0];
            M2[ i*test->RWG->nedge + j ]+= cres[1]; 
        }                                                                                                     
        else if(trial->RWG->sing4[m]==2)                                                                            
        {                                                                                                     
            singular_maxwell_hypersingular_ce(kappa, order, trial->RWG->sbf[i]->T2->v, test->RWG->sbf[j]->T2->v, trial->RWG->sbf[i]->p2, test->RWG->sbf[j]->p2, cres,NULL, trial->RWG->sbf[i]->len, test->RWG->sbf[j]->len, trial->RWG->sbf[i]->T2->a, test->RWG->sbf[j]->T2->a, trial->RWG->sbf[i]->T2->edges, test->RWG->sbf[j]->T2->edges, 1, 1);
            M1[ i*test->RWG->nedge + j ]+= cres[0];                                                         
            M2[ i*test->RWG->nedge + j ]+= cres[1];
        } 
        else if(trial->RWG->sing4[m]==1)                                                                       
        {                                                                                                     
            singular_maxwell_hypersingular_cv(kappa, order, trial->RWG->sbf[i]->T2->v, test->RWG->sbf[j]->T2->v, trial->RWG->sbf[i]->p2, test->RWG->sbf[j]->p2, cres,NULL, trial->RWG->sbf[i]->len, test->RWG->sbf[j]->len, trial->RWG->sbf[i]->T2->a, test->RWG->sbf[j]->T2->a, trial->RWG->sbf[i]->T2->dofs, test->RWG->sbf[j]->T2->dofs, 1, 1);
            M1[ i*test->RWG->nedge + j ]+= cres[0];                                                           
            M2[ i*test->RWG->nedge + j ]+= cres[1]; 
        } 
        else if(trial->RWG->sing4[m]==0) 
        {
            integration_elements_PRIMAL( sz[0],sz[1], eta, chi, trial->RWG->sbf[i]->T2->v, IE_x, P0_x );       
            divergence(trial->RWG->sbf[i]->T2->a, trial->RWG->sbf[i]->len, trial->RWG->sbf[i]->T2->v,"RWG",1, divx, length_arr);
            vec_prod( length_arr, divx, weights, IE_x, di ); 
            integration_elements_PRIMAL( sz[0],sz[1], eta, chi , test->RWG->sbf[j]->T2->v, IE_y, P0_y );          
            divergence(test->RWG->sbf[j]->T2->a, test->RWG->sbf[j]->len, test->RWG->sbf[j]->T2->v,"RWG",1, divy, length_arr);
            vec_prod( length_arr, divy, weights, IE_y, dj );                                                      
            evaluate_kernel( P0_x,P0_y, kappa, length_arr, length_arr, di, dj, cres );                            
            M1[ i*test->RWG->nedge + j ]+= cres[0];                                                               
            M2[ i*test->RWG->nedge + j ]+= cres[1]; 
        }

    }

} 

void assemble_maxwell_identity(// Input:                                                                                    
                 // Space data                                                                                
                 pMeshes Mesh,                                                                                
                 // wavenumber                                                                                
                 double kappa,                                                                                
                 // Order of quadrature                                                                       
                 int *order,                                                                                  
                 // Output: complex and real part                                                             
                 double *M1,                                                                                                                                                                  
                 char *trialbf,                                                                               
                 char *testbf )                                                                               
{                                                                                                             
    double weights[4], eta[4], eta2[4], chi2[4], chi[4];                                                      
    int sz[2], sz2[2], i, j, k, m, length_arr;                                                              
    quad_rules_PRIMAL(eta, chi, sz, order[0]);                                                                
    quad_rules_PRIMAL(eta2, chi2, sz2, order[1]);                                                             
    weights_PRIMAL(weights, order[0], sz[0]);                                                                 
    length_arr = sz[0];                                                                                       
    double P0_x[3*length_arr], IE_x[length_arr];                        
    double *phii = (double*)malloc(3*length_arr*sizeof(double));                                              
    double *phij = (double*)malloc(3*length_arr*sizeof(double));                                                                                                                                      
    double dj[3*length_arr];                                                                


    if(strcmp(trialbf,"RWG")==0 && strcmp(testbf,"RWG")==0)
    {
    for( m = 0; m < Mesh->RWG->nedge*Mesh->RWG->nedge; m++ )                                                  
    {                                                                                                         
        i = m/Mesh->RWG->nedge;                                                                               
        j = m%Mesh->RWG->nedge;                                                                               
        integration_elements_PRIMAL( sz[0],sz[1], eta, chi, Mesh->RWG->sbf[i]->T1->v, IE_x, P0_x );           
        ev_basis_function_maxwell(P0_x, phii, Mesh->RWG->sbf[i]->p1, Mesh->RWG->sbf[i]->T1->a, Mesh->RWG->sbf[i]->len, Mesh->RWG->sbf[i]->T1->n, length_arr, 0, "RWG");

        if( Mesh->RWG->sbf[i]->num[0] == Mesh->RWG->sbf[j]->num[0] )                                          
        {                                                                                                        
            ev_basis_function_maxwell(P0_x, phij, Mesh->RWG->sbf[j]->p1, Mesh->RWG->sbf[j]->T1->a, Mesh->RWG->sbf[j]->len, Mesh->RWG->sbf[j]->T1->n, length_arr, 0, "RWG");
            vec_prod_maxwell( length_arr, weights, IE_x, phij, dj );                                          
            M1[ i*Mesh->RWG->nedge + j ]+= dot_prod2(length_arr, dj, phii);                                                                                                                   
        }
        else if( Mesh->RWG->sbf[i]->num[0] == Mesh->RWG->sbf[j]->num[1] )                                          
        {                                                                                                        
            ev_basis_function_maxwell(P0_x, phij, Mesh->RWG->sbf[j]->p2, Mesh->RWG->sbf[j]->T2->a, Mesh->RWG->sbf[j]->len, Mesh->RWG->sbf[j]->T2->n, length_arr, 1, "RWG");
            vec_prod_maxwell( length_arr, weights, IE_x, phij, dj );                                          
            M1[ i*Mesh->RWG->nedge + j ]+= dot_prod2(length_arr, dj, phii);                                                          
        }

        integration_elements_PRIMAL( sz[0],sz[1], eta, chi, Mesh->RWG->sbf[i]->T2->v, IE_x, P0_x );           
        ev_basis_function_maxwell(P0_x, phii, Mesh->RWG->sbf[i]->p2, Mesh->RWG->sbf[i]->T2->a, Mesh->RWG->sbf[i]->len, Mesh->RWG->sbf[i]->T2->n, length_arr, 1, "RWG");
        if( Mesh->RWG->sbf[i]->num[1] == Mesh->RWG->sbf[j]->num[0] )                                          
        {                                                                                                     
            ev_basis_function_maxwell(P0_x, phij, Mesh->RWG->sbf[j]->p1, Mesh->RWG->sbf[j]->T1->a, Mesh->RWG->sbf[j]->len, Mesh->RWG->sbf[j]->T1->n, length_arr, 0, "RWG");
            vec_prod_maxwell( length_arr, weights, IE_x, phij, dj );                                          
            M1[ i*Mesh->RWG->nedge + j ]+= dot_prod2(length_arr,dj, phii);                                                                                                      
        }                                                                                                     
                                                                                                              
        else if( Mesh->RWG->sbf[i]->num[1] == Mesh->RWG->sbf[j]->num[1] )                                          
        {                                                                                                     
            ev_basis_function_maxwell(P0_x, phij, Mesh->RWG->sbf[j]->p2, Mesh->RWG->sbf[j]->T2->a, Mesh->RWG->sbf[j]->len, Mesh->RWG->sbf[j]->T2->n, length_arr, 1, "RWG");
            vec_prod_maxwell( length_arr, weights, IE_x, phij, dj );                                          
            M1[ i*Mesh->RWG->nedge + j ]+= dot_prod2(length_arr, dj, phii);                                                                                                                  
        } 

    } 
    }
    else if(strcmp(trialbf,"RWG")==0 && strcmp(testbf,"P0")==0) 
    {
        for(i=0; i< Mesh->RWG->nedge; i++)
        {
         
            //if(i==0)
            //{
            integration_elements_PRIMAL( sz[0],sz[1], eta, chi, Mesh->RWG->sbf[i]->T1->v, IE_x, P0_x ); 
            for(k=0;k<length_arr;k++)
            {
                M1[ i*Mesh->RWG->nedge + Mesh->RWG->sbf[i]->num[0] ] += weights[k]*(P0_x[3*k]-Mesh->RWG->sbf[i]->p1[0]+P0_x[3*k+1]-Mesh->RWG->sbf[i]->p1[1]+ P0_x[3*k+2]-Mesh->RWG->sbf[i]->p1[2])*IE_x[k]*Mesh->RWG->sbf[i]->len/(2* Mesh->RWG->sbf[i]->T1->a);
                //printf("%f %f %f %f %f %f\n", P0_x[3*k]-Mesh->RWG->sbf[i]->p1[0],P0_x[3*k+1]-Mesh->RWG->sbf[i]->p1[1], P0_x[3*k+2]-Mesh->RWG->sbf[i]->p1[2],eta[k],chi[k], IE_x[k]);
                
            }
            integration_elements_PRIMAL( sz[0],sz[1], eta, chi, Mesh->RWG->sbf[i]->T2->v, IE_x, P0_x );
            for(k=0;k<length_arr;k++)                                                                      
            { 
                //printf("%f %f %f %f %f %f\n",P0_x[3*k]-Mesh->RWG->sbf[i]->p2[0],P0_x[3*k+1]-Mesh->RWG->sbf[i]->p2[1], P0_x[3*k+2]-Mesh->RWG->sbf[i]->p2[2], eta[k], chi[k], IE_x[k]);
                M1[ i*Mesh->RWG->nedge + Mesh->RWG->sbf[i]->num[1] ] += weights[k]*(P0_x[3*k]-Mesh->RWG->sbf[i]->p2[0]+P0_x[3*k+1]-Mesh->RWG->sbf[i]->p2[1]+ P0_x[3*k+2]-Mesh->RWG->sbf[i]->p2[2])*IE_x[k]*Mesh->RWG->sbf[i]->len/(2* Mesh->RWG->sbf[i]->T2->a);
            }                                                                                         
            //}
        }     
    
    }
}

void assemble_maxwell_laplace_beltrami(// Input:                                                                                    
                 // Space data                                                                                
                 pMeshes Mesh,                                                                                
                 // wavenumber                                                                                
                 double kappa,                                                                                
                 // Order of quadrature                                                                       
                 int *order,                                                                                  
                 // Output: complex and real part                                                             
                 double *M1,                                                                                  
                 char *trialbf,                                                                               
                 char *testbf ) 
{
    double weights[4], eta[4], eta2[4], chi2[4], chi[4];                                                      
    int sz[2], sz2[2], i, j, m, length_arr;                                                              
    quad_rules_PRIMAL(eta, chi, sz, order[0]);                                                                
    quad_rules_PRIMAL(eta2, chi2, sz2, order[1]);                                                             
    weights_PRIMAL(weights, order[0], sz[0]);                                                                 
    length_arr = sz[0];                                                                                       
    double P0_x[3*length_arr], IE_x[length_arr];                                                                                                             
    double dj[length_arr];                                                                                                                                                            
    double divx[length_arr], divy[length_arr];                                                                
    for( m = 0; m < Mesh->RWG->nedge*Mesh->RWG->nedge; m++ )                                                  
    {                                                                                                         
        i = m/Mesh->RWG->nedge;                                                                               
        j = m%Mesh->RWG->nedge;                                                                               
        integration_elements_PRIMAL( sz[0],sz[1], eta, chi, Mesh->RWG->sbf[i]->T1->v, IE_x, P0_x );           
        divergence(Mesh->RWG->sbf[i]->T1->a,Mesh->RWG->sbf[i]->len, Mesh->RWG->sbf[i]->T1->v,"RWG",0, divx, length_arr);
        if( Mesh->RWG->sbf[i]->num[0] == Mesh->RWG->sbf[j]->num[0] )                                          
        {                                                                                                     
            divergence(Mesh->RWG->sbf[j]->T1->a,Mesh->RWG->sbf[j]->len, Mesh->RWG->sbf[j]->T1->v,"RWG",0, divy, length_arr);
            vec_prod( length_arr, divy, weights, IE_x, dj );                                                  
            M1[ i*Mesh->RWG->nedge + j ]+= dot_prod(length_arr, 0, 0, dj, divx);                                                                                                     
        }                                                                                                                                       
        else if( Mesh->RWG->sbf[i]->num[0] == Mesh->RWG->sbf[j]->num[1] )                                          
        {                                                                                                       
            divergence(Mesh->RWG->sbf[j]->T2->a,Mesh->RWG->sbf[j]->len, Mesh->RWG->sbf[j]->T2->v,"RWG",1, divy, length_arr);
            vec_prod( length_arr, divy, weights, IE_x, dj );                                                                   
            M1[ i*Mesh->RWG->nedge + j ]+= dot_prod(length_arr, 0, 0, dj, divx);                                                                                          
        }                                                                                                     
        integration_elements_PRIMAL( sz[0],sz[1], eta, chi, Mesh->RWG->sbf[i]->T2->v, IE_x, P0_x );           
        divergence(Mesh->RWG->sbf[i]->T2->a,Mesh->RWG->sbf[i]->len, Mesh->RWG->sbf[i]->T2->v,"RWG",1, divx, length_arr);
        if( Mesh->RWG->sbf[i]->num[1] == Mesh->RWG->sbf[j]->num[0] )                                          
        {                                                                                                        
            divergence(Mesh->RWG->sbf[j]->T1->a,Mesh->RWG->sbf[j]->len, Mesh->RWG->sbf[j]->T1->v,"RWG",0, divy, length_arr);
            vec_prod( length_arr, divy, weights, IE_x, dj );                                                                        
            M1[ i*Mesh->RWG->nedge + j ]+= dot_prod(length_arr, 0, 0, dj, divx);                              
        } 
         else if( Mesh->RWG->sbf[i]->num[1] == Mesh->RWG->sbf[j]->num[1] )                                          
        {                                                                                                        
            divergence(Mesh->RWG->sbf[j]->T2->a,Mesh->RWG->sbf[j]->len, Mesh->RWG->sbf[j]->T2->v,"RWG",1, divy, length_arr);
            vec_prod( length_arr, divy, weights, IE_x, dj );                                                                         
            M1[ i*Mesh->RWG->nedge + j ]+= dot_prod(length_arr, 0, 0, dj, divx);                              
        }                                                                                                     
                                                                                                              
    }
}

gsl_complex assemble_maxwell_dual_aux2( psRWG trial, psRWG test, double kappa, int *order)
{
    double weights[4], weights2[4],eta[4], eta2[4], chi2[4], chi[4];                                          
    int sz[2], sz2[2], i, j, m, length_arr;                                                                   
    quad_rules_PRIMAL(eta, chi, sz, order[0]);                                                                
    quad_rules_PRIMAL(eta2, chi2, sz2, order[1]);                                                             
    weights_PRIMAL(weights2, order[1], sz2[0]);                                                               
    weights_PRIMAL(weights, order[0], sz[0]);                                                                 
    length_arr = sz[0];                                                                                       
    double P0_x[3*length_arr], P0_y[3*length_arr], IE_x[length_arr], IE_y[length_arr];                        
    double di[length_arr], dj[length_arr];                                                                    
    double di2[3*length_arr], dj2[3*length_arr];                                                              
    double cres1[2], cres2[2];                                                                                
    double divx[length_arr], divy[length_arr];                                                                
    double phii[length_arr], phij[length_arr];  

        integration_elements_PRIMAL( sz[0],sz[1], eta, chi, trial->T1->v, IE_x, P0_x );          
        divergence(trial->RWG->sbf[i]->T1->a,trial->RWG->sbf[i]->len, trial->T1->v,"RWG",0, divx, length_arr);
        vec_prod( length_arr, divx, weights, IE_x, di );                                                      
        ev_basis_function_maxwell(P0_x, phii, trial->p1, trial->T1->a, trial->len, trial->T1->n, length_arr, 0, "RWG");
        vec_prod_maxwell( length_arr, weights, IE_x, phii, di2 );                                             
        integration_elements_PRIMAL( sz[0],sz[1], eta, chi , test->T1->v, IE_y, P0_y );          
        divergence(test->T1->a, test->len, test->T1->v,"RWG",0, divy, length_arr);
        vec_prod( length_arr, divy, weights, IE_y, dj );                                                      
        ev_basis_function_maxwell(P0_y, phij, test->p1, test->T1->a, test->len, test->T1->n, length_arr, 0, "RWG");
        vec_prod_maxwell( length_arr, weights, IE_y, phij, dj2 );                                             
        evaluate_kernel_maxwell_T( P0_x,P0_y, kappa, length_arr, length_arr, di, dj, di2, dj2, cres1, cres2 );
                                                                                                              
        integration_elements_PRIMAL( sz[0],sz[1], eta, chi , test->T2->v, IE_y, P0_y );          
        divergence(test->T2->a, test->len, test->T2->v,"RWG",1, divy, length_arr);
        vec_prod( length_arr, divy, weights, IE_y, dj );                                                      
        ev_basis_function_maxwell(P0_y, phij, test->p2, test->T2->a, test->len, test->T2->n, length_arr, 1, "RWG");
        vec_prod_maxwell( length_arr, weights, IE_y, phij, dj2 );                                             
        evaluate_kernel_maxwell_T( P0_x,P0_y, kappa, length_arr, length_arr, di, dj, di2, dj2, cres1, cres2 );
                                                                                                              
        integration_elements_PRIMAL( sz[0],sz[1], eta, chi, trial->T2->v, IE_x, P0_x );          
        divergence(trial->T2->a, trial->len, trial->T2->v,"RWG",1, divx, length_arr);
        vec_prod( length_arr, divx, weights, IE_x, di );                                                      
        ev_basis_function_maxwell(P0_x, phii, trial->p2, trial->T2->a, trial->len, trial->T2->n, length_arr, 1, "RWG");
        vec_prod_maxwell( length_arr, weights, IE_x, phii, di2 );                                             
        integration_elements_PRIMAL( sz[0],sz[1], eta, chi , test->T1->v, IE_y, P0_y );          
        divergence(test->T1->a, test->len, test->T1->v,"RWG",0, divy, length_arr);
        vec_prod( length_arr, divy, weights, IE_y, dj );                                                      
        ev_basis_function_maxwell(P0_y, phij, test->p1, test->T1->a, test->len, test->T1->n, length_arr, 0, "RWG");
        vec_prod_maxwell( length_arr, weights, IE_y, phij, dj2 );                                             
        evaluate_kernel_maxwell_T( P0_x,P0_y, kappa, length_arr, length_arr, di, dj, di2, dj2, cres1, cres2 );

        integration_elements_PRIMAL( sz[0],sz[1], eta, chi , test->T2->v, IE_y, P0_y );          
        divergence(test->T2->a, test->len, test->T2->v,"RWG",1, divy, length_arr);
        vec_prod( length_arr, divy, weights, IE_y, dj );                                                      
        ev_basis_function_maxwell(P0_y, phij, test->p2, test->T2->a, test->len, test->T2->n, length_arr, 1, "RWG");
        vec_prod_maxwell( length_arr, weights, IE_y, phij, dj2 );                                             
        evaluate_kernel_maxwell_T( P0_x,P0_y, kappa, length_arr, length_arr, di, dj, di2, dj2, cres1, cres2 );
}

gsl_complex assemble_maxwell_dual_aux(psRWG trial, psBC test, double kappa, int *order)
{

    int i;
    gsl_vector_complex *res1 = gsl_vector_complex_alloc(15);                                                  
    gsl_vector_complex *res2 = gsl_vector_complex_alloc(15);                                                  
    gsl_vector_complex *res3 = gsl_vector_complex_alloc(15);                                                  
    gsl_vector_complex *res4 = gsl_vector_complex_alloc(15);                                                  
    gsl_vector_complex *res5 = gsl_vector_complex_alloc(2); 
    
    for( i = 0; i< test->len1; i++)
    {
        assemble_maxwell_dual_aux2( trial, test->sbf1[i], kappa, order);
    }
    for( i = 0; i< test->len2; i++)                                                                           
    {                                          
        assemble_maxwell_dual_aux2( trial, test->sbf2[i], kappa, order);                                                                                                      
    } 
    for( i = 0; i< test->len3; i++)                                                                           
    {                                                                                                         
        assemble_maxwell_dual_aux2( trial, test->sbf3[i], kappa, order);                                                          
    } 
    for( i = 0; i< test->len4; i++)                                                                           
    {                                                                                                         
        assemble_maxwell_dual_aux2( trial, test->sbf4[i], kappa, order);                                                                             
    } 
    for( i = 0; i< 2; i++)                                                                           
    {                                                                                                         
        assemble_maxwell_dual_aux2( trial, test->sbf5[i], kappa, order);                                                
    } 
}

gsl_complex sing_maxwell_dual(psBC bf, gsl_matrix_complex *Aux)
{

    int i;
    for(i = 0; i< bc->len1; i++)
    {
        singular_maxwell_hypersingular2(kappa, order, bf->sbf1[i]->T1->v, bf->sbf1[i]->p1, bf->sbf1[i]->p1, cres1, bf->sbf1[i]->T1->n, bf->sbf1[i]->T1->a, bf->sbf1[i]->len, bf->sbf1[i]->len, 0, 0);
        singular_maxwell_weakly2(kappa, order, bf->sbf1[i]->T1->v, bf->sbf1[i]->p1, bf->sbf1[i]->p1, cres2, bf->sbf1[i]->T1->n, bf->sbf1[i]->T1->a, bf->sbf1[i]->len, bf->sbf1[i]->len, 0, 0);
        singular_maxwell_hypersingular_ce(kappa, order, bf->sbf1[i]->T1->v, bf->sbf1[i]->T2->v, bf->sbf1[i]->p1, bf->sbf1[i]->p2, cres1,NULL, bf->sbf1[i]->len, bf->sbf1[i]->len, bf->sbf1[i]->T1->a, bf->sbf1[i]->T2->a, bf->sbf1[i]->T1->edges, bf->sbf[i]->T2->edges, 0, 0);
        singular_maxwell_weakly_ce(kappa, order, bf->sbf1[i]->T1->v, bf->sbf1[i]->T2->v, bf->sbf1[i]->p1, bf->sbf1[i]->p2, cres1,NULL, bf->sbf1[i]->len, bf->sbf1[i]->len, bf->sbf1[i]->T1->a, bf->sbf1[i]->T2->a, bf->sbf1[i]->T1->edges, bf->sbf[i]->T2->edges, 0, 0);
    }
    

}

void assemble_maxwell_dual(pSpace trial,                                                                                
                 pSpace test,                                                                                 
                 // wavenumber                                                                                
                 double kappa,                                                                                
                 // Order of quadrature                                                                       
                 int *order,                                                                                  
                 // Output: complex and real part                                                             
                 double *M1,                                                                                  
                 double *M2 )
{

    int i,j, k;
    gsl_vector_complex *res1 = gsl_vector_complex_alloc(15);
    gsl_vector_complex *res2 = gsl_vector_complex_alloc(15);
    gsl_vector_complex *res3 = gsl_vector_complex_alloc(15);
    gsl_vector_complex *res4 = gsl_vector_complex_alloc(15);
    gsl_vector_complex *res5 = gsl_vector_complex_alloc(2);
    gsl_matrix_complex *Aux = gsl_matrix_complex_calloc(trial->nvp, trial->nvp);

    for(i=0; i< trial->BC->nedge; i++)                                                                        
    {

    }
    for(i=0; i< trial->BC->nedge; i++)
    {
        for(j=0; j<test->BC->nedge; j++)
        {
            for(k=0; k<trial->BC->sbf[i]->len1; k++)
            {
                assemble_maxwell_dual_aux(trial->BC->sbf[i]->sbf1[k], test->BC->sbf[j], kappa, order);
            }
            for(k=0; k<trial->BC->sbf[i]->len2; k++)                                                                     
            {                                                                                                 
                assemble_maxwell_dual_aux(trial->BC->sbf[i]->sbf2[k], test->BC->sbf[j], kappa, order);                                          
            } 
            for(k=0; k<trial->BC->sbf[i]->len2; k++)                                                                     
            {                                                                                                 
                assemble_maxwell_dual_aux(trial->BC->sbf[i]->sbf3[k], test->BC->sbf[j], kappa, order);                                             
            } 
            for(k=0; k<trial->BC->sbf[i]->len4; k++)                                                                     
            {                                                                                                 
                assemble_maxwell_dual_aux(trial->BC->sbf[i]->sbf4[k], test->BC->sbf[j], kappa, order);                                             
            } 
            for(k=0; k<2; k++)                                                                     
            {                                                                                                 
                assemble_maxwell_dual_aux(trial->BC->sbf[i]->sbf5[k], test->BC->sbf[j], kappa, order);                                             
            } 

        }
    }
}
void project_mat(                                                                                             
                 // Input:                                                                                    
                 // Space data                                                                                
                 pMeshes Mesh,                                                                                                                                                                
                 // Output: complex and real part                                                             
                 double *Matr,
                 double *Mati,
                 double *M1,                                                                                  
                 double *M2,                                                                                  
                 char *trialbf,                                                                               
                 char *testbf                                                                                 
                 ){
    
    int i,j,k,l;
    if( strcmp(testbf,"P1d") == 0 && strcmp(trialbf,"P1d") == 0)                                         
    {                                                                                                                                          
        gsl_complex aux1, aux2;                                                                               
        aux1 = gsl_complex_rect(0.0,0.0);                                                                     
        aux2 = gsl_complex_rect(0.0,0.0);                                                                     
        int idx1[7], idx2[7];                                                                                 
        gsl_vector_complex *f1 = gsl_vector_complex_alloc(7);                                                 
        gsl_vector_complex *f2 = gsl_vector_complex_alloc(7);                                                 
        gsl_vector_complex *f3 = gsl_vector_complex_alloc(7);                                                 
        gsl_vector_complex *f4 = gsl_vector_complex_alloc(7);                                                 
        for( i = 0; i < Mesh->PD1->nt; i++ )                                                                  
        {                                                                                                     
            for( j = 0; j < Mesh->PD1->nt; j++ )                                                              
            {                                                                                                 
                idx2[0] = Mesh->P1P->globaldofs[Mesh->PD1->T[j]->dofs[0]];                                    
                idx2[1] = Mesh->P1P->globaldofs[Mesh->PD1->T[j]->dofs[1]];                                    
                idx2[2] = Mesh->P1P->globaldofs[Mesh->PD1->T[j]->dofs[2]];                                    
                idx2[3] = Mesh->P1P->globaldofs[Mesh->PD1->T[j]->T[0]->dofs[2]];                              
                idx2[4] = Mesh->P1P->globaldofs[Mesh->PD1->T[j]->T[0]->dofs[1]];                              
                idx2[5] = Mesh->P1P->globaldofs[Mesh->PD1->T[j]->T[5]->dofs[2]];                              
                idx2[6] = Mesh->P1P->globaldofs[Mesh->PD1->T[j]->T[3]->dofs[2]];                              
                idx1[0] = Mesh->P1P->globaldofs[Mesh->PD1->T[i]->dofs[0]];                                    
                idx1[1] = Mesh->P1P->globaldofs[Mesh->PD1->T[i]->dofs[1]];                                    
                idx1[2] = Mesh->P1P->globaldofs[Mesh->PD1->T[i]->dofs[2]];                                    
                idx1[3] = Mesh->P1P->globaldofs[Mesh->PD1->T[i]->T[0]->dofs[2]];                              
                idx1[4] = Mesh->P1P->globaldofs[Mesh->PD1->T[i]->T[0]->dofs[1]];                              
                idx1[5] = Mesh->P1P->globaldofs[Mesh->PD1->T[i]->T[5]->dofs[2]];                              
                idx1[6] = Mesh->P1P->globaldofs[Mesh->PD1->T[i]->T[3]->dofs[2]];                              
                gsl_vector_complex_set(f1, 0, gsl_complex_rect(1.0/(Mesh->PD1->lensbf[i][0]),0));             
                gsl_vector_complex_set(f1, 1, gsl_complex_rect(1.0/(Mesh->PD1->lensbf[i][1]),0));             
                gsl_vector_complex_set(f1, 2, gsl_complex_rect(1.0/(Mesh->PD1->lensbf[i][2]),0));             
                gsl_vector_complex_set(f1, 3, gsl_complex_rect(1.0,0));                                       
                gsl_vector_complex_set(f1, 4, gsl_complex_rect(0.5,0));                                       
                gsl_vector_complex_set(f1, 5, gsl_complex_rect(0.5,0));                                       
                gsl_vector_complex_set(f1, 6, gsl_complex_rect(0.5,0));                                       
                gsl_vector_complex_set(f2, 0, gsl_complex_rect(1.0/(Mesh->PD1->lensbf[j][0]),0));             
                gsl_vector_complex_set(f2, 1, gsl_complex_rect(1.0/(Mesh->PD1->lensbf[j][1]),0));             
                gsl_vector_complex_set(f2, 2, gsl_complex_rect(1.0/(Mesh->PD1->lensbf[j][2]),0));             
                gsl_vector_complex_set(f2, 3, gsl_complex_rect(1.0,0));                                       
                gsl_vector_complex_set(f2, 4, gsl_complex_rect(0.5,0));                                       
                gsl_vector_complex_set(f2, 5, gsl_complex_rect(0.5,0));                                       
                gsl_vector_complex_set(f2, 6, gsl_complex_rect(0.5,0));                                       
                                                                                                              
                for(k=0;k<7;k++)                                                                              
                {                                                                                             
                    for(l=0;l<7;l++)                                                                          
                    {                                                                                                                                                                          
                       gsl_vector_complex_set(f3, l,gsl_complex_rect(Matr[idx1[k]*Mesh->P1P->ndof +idx2[l]],Mati[idx1[k]*Mesh->P1P->ndof +idx2[l]]));
                                                                                                             
                    }                                                                                         
                    gsl_blas_zdotu(f3,f2,&aux1);                                                              
                    gsl_vector_complex_set(f4, k, aux1);                                                      
                }                                                                                             
                                                                                                              
                gsl_blas_zdotu(f4,f1,&aux2);                                                                  
                M1[i*Mesh->PD1->nt+j] = GSL_REAL(aux2);                                                       
                M2[i*Mesh->PD1->nt+j] = GSL_IMAG(aux2);                                                       
                                                                                                              
            }                                                                                                 
        }                                                                                                                           
    } 
}
gsl_complex assemble_P1_HP_entries(int *order, double kappa, int i, int j, double *sing, pMeshes Mesh)
{

    double weights[4], weights2[4], eta[4], eta2[4], chi2[4], chi[4];                                                      
    int sz[2], sz2[2], k, l, length_arr;                                                              
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
    double cres[2];  
    gsl_complex aux = gsl_complex_rect(0.0,0.0);
    gsl_complex aux2 = gsl_complex_rect(0.0,0.0);
   for( k = 0; k <Mesh->P1P->sbf[i]->len; k++ )                                                 
   {                                                                                             
        integration_elements_PRIMAL( sz[0],sz[1], eta, chi, Mesh->P1P->sbf[i]->T[k]->v, IE_x, P0_x ); 
        //n_x = normal(P0_x, Mesh->P1P->sbf[i]->T[k]->v);                                                                         
        //curl_x = curl_bf(Mesh->P1P->sbf[i]->T[k]->v, 0, n_x);  
        for(l = 0; l < Mesh->P1P->sbf[j]->len; l++)                                              
        {  
           if( Mesh->P1P->sbf[i]->num[k] != Mesh->P1P->sbf[j]->num[l] )                        
           {     
                integration_elements_PRIMAL( sz[0],sz[1], eta, chi , Mesh->P1P->sbf[j]->T[l]->v, IE_y, P0_y );
                res = HP_kern(kappa, weights, weights, phii, phij, P0_x, P0_y, IE_x, IE_y, Mesh->P1P->sbf[i]->T[k]->n, Mesh->P1P->sbf[j]->T[l]->n, Mesh->P1P->sbf[i]->T[k]->c, Mesh->P1P->sbf[j]->T[l]->c, length_arr, length_arr );    
                //n_y = normal(P0_y, Mesh->P1P->sbf[j]->T[l]->v);                                       
                //curl_y = curl_bf(Mesh->P1P->sbf[j]->T[l]->v, 0, n_y);                                 
                //res = HP_kern(kappa, weights, weights, phii, phij, P0_x, P0_y, IE_x, IE_y, n_x, n_y, curl_x, curl_y, length_arr, length_arr );
                aux = gsl_complex_add(aux, res);
           }                                                                                     
           else                                                                                  
           {  
                double resc, resn;                                                                    
                integration_elements_PRIMAL( sz2[0],sz2[1], eta2, chi2, Mesh->P1P->sbf[i]->T[k]->v, IE_x, P0_x );                                
                integration_elements_PRIMAL( sz[0],sz[1], eta, chi , Mesh->P1P->sbf[j]->T[l]->v, IE_y, P0_y );
                singular(P0_x, Mesh->P1P->sbf[j]->T[l]->v, weights2, cres, IE_x,sz2, order[2], 0, kappa, "P1", phij2[0]);  
                gsl_blas_ddot(Mesh->P1P->sbf[i]->T[k]->n, Mesh->P1P->sbf[j]->T[l]->n, &resn);   
                //n_y = normal(P0_y, Mesh->P1P->sbf[j]->T[l]->v);                                       
                //curl_y = curl_bf(Mesh->P1P->sbf[j]->T[l]->v, 0, n_y);                                 
                //gsl_blas_ddot(n_x, n_y, &resn);  
                res =  gsl_complex_rect(-resn*kappa*kappa*cres[0],-resn*kappa*kappa*cres[1]);
                aux2 = gsl_complex_add(aux2, res);                                                          
                singular(P0_x, Mesh->P1P->sbf[j]->T[l]->v, weights2, cres, IE_x, sz2, order[2], 1, kappa, "P0", 1.0);
                gsl_blas_ddot(Mesh->P1P->sbf[i]->T[k]->c, Mesh->P1P->sbf[j]->T[l]->c, &resc);
                //gsl_blas_ddot(curl_x, curl_y, &resc);
                res = gsl_complex_rect(cres[0]*resc,cres[1]*resc);
                aux2 = gsl_complex_add(aux2, res);  
           
           }                                                                                     
        }                                                                                         
     } 
     
     sing[0] = GSL_REAL(aux2);
     sing[1] = GSL_IMAG(aux2);   
     return aux;
}

gsl_complex s_part(double kappa, int *order, pMeshes Mesh, int i, int k, int j, int l, int opt)
{
    
    double weights[4], weights2[4], eta[4], eta2[4], chi2[4], chi[4];
    int sz[2], sz2[2], length_arr;
    quad_rules_PRIMAL(eta, chi, sz, order[0]);
    quad_rules_PRIMAL(eta2, chi2, sz2, order[1]);
    weights_PRIMAL(weights, order[0], sz[0]);
    weights_PRIMAL(weights2, order[1], sz2[0]); 
    length_arr = sz[0];
    double P0_x[3*length_arr], P0_y[3*length_arr], IE_x[length_arr], IE_y[length_arr];
    double phii[16] = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
    double phij[16] = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
    gsl_complex res;
    gsl_vector *n_x;
    gsl_vector *n_y;
    gsl_vector *curl_x;
    gsl_vector *curl_y;
    double cres[2];
    gsl_complex aux = gsl_complex_rect(0.0,0.0);
    double resc, resn;
    
    if (opt==0)
    {
        integration_elements_PRIMAL( sz[0],sz[1], eta, chi, Mesh->P1P->sbf[i]->T[k]->v, IE_x, P0_x );
        n_x = normal(P0_x, Mesh->P1P->sbf[i]->T[k]->v);
        ev_basis_function( length_arr, 0, phii, eta, chi, "P1" );
        curl_x = curl_bf(Mesh->P1P->sbf[i]->T[k]->v, 0, n_x);
        integration_elements_PRIMAL( sz[0],sz[1], eta, chi , Mesh->P1P->sbf[i]->T[k+1]->v, IE_y, P0_y );
        ev_basis_function( length_arr, 0, phij, eta, chi, "P1" );
        n_y = normal(P0_y, Mesh->P1P->sbf[i]->T[k+1]->v);
        curl_y = curl_bf(Mesh->P1P->sbf[i]->T[k+1]->v, 0, n_y);
        res = HP_kern(kappa, weights, weights, phii, phij, P0_x, P0_y, IE_x, IE_y, n_x, n_y, curl_x, curl_y, length_arr, length_arr );
        aux = gsl_complex_add(aux, res);
                                                                             
        integration_elements_PRIMAL( sz[0],sz[1], eta, chi, Mesh->P1P->sbf[i]->T[k+1]->v, IE_x, P0_x );
        n_x = normal(P0_x, Mesh->P1P->sbf[i]->T[k+1]->v);
        ev_basis_function( length_arr, 0, phii, eta, chi, "P1" );
        curl_x = curl_bf(Mesh->P1P->sbf[i]->T[k+1]->v, 0, n_x);
        integration_elements_PRIMAL( sz[0],sz[1], eta, chi , Mesh->P1P->sbf[i]->T[k]->v, IE_y, P0_y );
        ev_basis_function( length_arr, 0, phij, eta, chi, "P1" );
        n_y = normal(P0_y, Mesh->P1P->sbf[i]->T[k]->v);
        curl_y = curl_bf(Mesh->P1P->sbf[i]->T[k]->v, 0, n_y);
        res = HP_kern(kappa, weights, weights, phii, phij, P0_x, P0_y, IE_x, IE_y, n_x, n_y, curl_x, curl_y, length_arr, length_arr );
        aux = gsl_complex_add(aux, res);

        integration_elements_PRIMAL( sz[0],sz[1], eta, chi, Mesh->P1P->sbf[i]->T[k]->v, IE_x, P0_x );         
        n_x = normal(P0_x, Mesh->P1P->sbf[i]->T[k]->v);                                                       
        ev_basis_function( length_arr, 0, phii, eta, chi, "P1" );                                             
        curl_x = curl_bf(Mesh->P1P->sbf[i]->T[k]->v, 0, n_x); 
        integration_elements_PRIMAL( sz2[0],sz2[1], eta2, chi2, Mesh->P1P->sbf[i]->T[k]->v, IE_x, P0_x );
        ev_basis_function( 1, 0, phij, eta2, chi2, "P1" );                                            
        integration_elements_PRIMAL( sz[0],sz[1], eta, chi , Mesh->P1P->sbf[j]->T[l]->v, IE_y, P0_y );
        singular(P0_x, Mesh->P1P->sbf[j]->T[l]->v, weights2, cres, IE_x,sz2, order[2], 0, kappa, "P1", phij[0]);
        n_y = normal(P0_y, Mesh->P1P->sbf[j]->T[l]->v);                                               
        curl_y = curl_bf(Mesh->P1P->sbf[j]->T[l]->v, 0, n_y);                                         
        gsl_blas_ddot(n_x, n_y, &resn);                                                               
        aux = gsl_complex_add(aux, gsl_complex_rect(-resn*kappa*kappa*cres[0],-resn*kappa*kappa*cres[1]));
        singular(P0_x, Mesh->P1P->sbf[j]->T[l]->v, weights2, cres, IE_x, sz2, order[2], 1, kappa, "P0", 1.0);   
        gsl_blas_ddot(curl_x, curl_y, &resc);                                                         
        aux = gsl_complex_add(aux, gsl_complex_rect(cres[0]*resc,cres[1]*resc)); 
        
        integration_elements_PRIMAL( sz[0],sz[1], eta, chi, Mesh->P1P->sbf[i]->T[k+1]->v, IE_x, P0_x );         
        n_x = normal(P0_x, Mesh->P1P->sbf[i]->T[k+1]->v);                                                       
        ev_basis_function( length_arr, 0, phii, eta, chi, "P1" );                                             
        curl_x = curl_bf(Mesh->P1P->sbf[i]->T[k+1]->v, 0, n_x);                                                 
        integration_elements_PRIMAL( sz2[0],sz2[1], eta2, chi2, Mesh->P1P->sbf[i]->T[k+1]->v, IE_x, P0_x );     
        ev_basis_function( 1, 0, phij, eta2, chi2, "P1" );                                                    
        integration_elements_PRIMAL( sz[0],sz[1], eta, chi , Mesh->P1P->sbf[j]->T[l+1]->v, IE_y, P0_y );        
        singular(P0_x, Mesh->P1P->sbf[j]->T[l+1]->v, weights2, cres, IE_x,sz2, order[2], 0, kappa, "P1", phij[0]);        
        n_y = normal(P0_y, Mesh->P1P->sbf[j]->T[l+1]->v);                                                       
        curl_y = curl_bf(Mesh->P1P->sbf[j]->T[l+1]->v, 0, n_y);                                                 
        gsl_blas_ddot(n_x, n_y, &resn);                                                                       
        aux = gsl_complex_add(aux, gsl_complex_rect(-resn*kappa*kappa*cres[0],-resn*kappa*kappa*cres[1]));    
        singular(P0_x, Mesh->P1P->sbf[j]->T[l+1]->v, weights2, cres, IE_x, sz2, order[2], 1, kappa, "P0", 1.0);           
        gsl_blas_ddot(curl_x, curl_y, &resc);                                                                 
        aux = gsl_complex_add(aux, gsl_complex_rect(cres[0]*resc,cres[1]*resc)); 
    }
    else if(opt==1)
    {
        integration_elements_PRIMAL( sz[0],sz[1], eta, chi, Mesh->P1P->sbf[i]->T[k]->v, IE_x, P0_x );
        n_x = normal(P0_x, Mesh->P1P->sbf[i]->T[k]->v);
        ev_basis_function( length_arr, 0, phii, eta, chi, "P1" );
        curl_x = curl_bf(Mesh->P1P->sbf[i]->T[k]->v, 0, n_x);
        integration_elements_PRIMAL( sz[0],sz[1], eta, chi , Mesh->P1P->sbf[j]->T[l]->v, IE_y, P0_y );
        ev_basis_function( length_arr, 0, phij, eta, chi, "P1" );
        n_y = normal(P0_y, Mesh->P1P->sbf[j]->T[l]->v);
        curl_y = curl_bf(Mesh->P1P->sbf[j]->T[l]->v, 0, n_y);
        res = HP_kern(kappa, weights, weights, phii, phij, P0_x, P0_y, IE_x, IE_y, n_x, n_y, curl_x, curl_y, length_arr, length_arr );
        aux = gsl_complex_add(aux, res);
        
        ev_basis_function( 1, 0, phij, eta2, chi2, "P1" );
        integration_elements_PRIMAL( sz[0],sz[1], eta, chi , Mesh->P1P->sbf[i]->T[k]->v, IE_y, P0_y );
        singular(P0_x, Mesh->P1P->sbf[i]->T[k]->v, weights2, cres, IE_x,sz2, order[2], 0, kappa, "P1", phij[0]);
        n_y = normal(P0_y, Mesh->P1P->sbf[i]->T[k]->v);
        curl_y = curl_bf(Mesh->P1P->sbf[i]->T[k]->v, 0, n_y);
        gsl_blas_ddot(n_x, n_y, &resn);
        aux = gsl_complex_add(aux, gsl_complex_rect(-resn*kappa*kappa*cres[0],-resn*kappa*kappa*cres[1]));
        singular(P0_x, Mesh->P1P->sbf[i]->T[k]->v, weights2, cres, IE_x, sz2, order[2], 1, kappa, "P0", 1.0);
        gsl_blas_ddot(curl_x, curl_y, &resc);
        aux = gsl_complex_add(aux, gsl_complex_rect(cres[0]*resc,cres[1]*resc));
        
        integration_elements_PRIMAL( sz[0],sz[1], eta, chi, Mesh->P1P->sbf[i]->T[k+1]->v, IE_x, P0_x );
        n_x = normal(P0_x, Mesh->P1P->sbf[i]->T[k+1]->v);
        ev_basis_function( length_arr, 0, phii, eta, chi, "P1" );
        curl_x = curl_bf(Mesh->P1P->sbf[i]->T[k+1]->v, 0, n_x);
        integration_elements_PRIMAL( sz[0],sz[1], eta, chi , Mesh->P1P->sbf[j]->T[l]->v, IE_y, P0_y );
        ev_basis_function( length_arr, 0, phij, eta, chi, "P1" );
        n_y = normal(P0_y, Mesh->P1P->sbf[j]->T[l]->v);
        curl_y = curl_bf(Mesh->P1P->sbf[j]->T[l]->v, 0, n_y);
        res = HP_kern(kappa, weights, weights, phii, phij, P0_x, P0_y, IE_x, IE_y, n_x, n_y, curl_x, curl_y, length_arr, length_arr );
        aux = gsl_complex_add(aux, res);
        
        integration_elements_PRIMAL( sz[0],sz[1], eta, chi, Mesh->P1P->sbf[i]->T[k+1]->v, IE_x, P0_x );
        n_x = normal(P0_x, Mesh->P1P->sbf[i]->T[k+1]->v);
        ev_basis_function( length_arr, 0, phii, eta, chi, "P1" );
        curl_x = curl_bf(Mesh->P1P->sbf[i]->T[k+1]->v, 0, n_x);
        integration_elements_PRIMAL( sz[0],sz[1], eta, chi , Mesh->P1P->sbf[j]->T[l+1]->v, IE_y, P0_y );
        ev_basis_function( length_arr, 0, phij, eta, chi, "P1" );
        n_y = normal(P0_y, Mesh->P1P->sbf[j]->T[l+1]->v);
        curl_y = curl_bf(Mesh->P1P->sbf[j]->T[l+1]->v, 0, n_y);
        res = HP_kern(kappa, weights, weights, phii, phij, P0_x, P0_y, IE_x, IE_y, n_x, n_y, curl_x, curl_y, length_arr, length_arr );
        aux = gsl_complex_add(aux, res);
    }
    else if(opt==2)
    {
        integration_elements_PRIMAL( sz[0],sz[1], eta, chi, Mesh->P1P->sbf[i]->T[k]->v, IE_x, P0_x );
        n_x = normal(P0_x, Mesh->P1P->sbf[i]->T[k]->v);
        ev_basis_function( length_arr, 0, phii, eta, chi, "P1" );
        curl_x = curl_bf(Mesh->P1P->sbf[i]->T[k]->v, 0, n_x);
        integration_elements_PRIMAL( sz[0],sz[1], eta, chi , Mesh->P1P->sbf[j]->T[l]->v, IE_y, P0_y );
        ev_basis_function( length_arr, 0, phij, eta, chi, "P1" );
        n_y = normal(P0_y, Mesh->P1P->sbf[j]->T[l]->v);
        curl_y = curl_bf(Mesh->P1P->sbf[j]->T[l]->v, 0, n_y);
        res = HP_kern(kappa, weights, weights, phii, phij, P0_x, P0_y, IE_x, IE_y, n_x, n_y, curl_x, curl_y, length_arr, length_arr );
        aux = gsl_complex_add(aux, res);
        
        integration_elements_PRIMAL( sz[0],sz[1], eta, chi, Mesh->P1P->sbf[i]->T[k]->v, IE_x, P0_x );
        n_x = normal(P0_x, Mesh->P1P->sbf[i]->T[k]->v);
        ev_basis_function( length_arr, 0, phii, eta, chi, "P1" );
        curl_x = curl_bf(Mesh->P1P->sbf[i]->T[k]->v, 0, n_x);
        integration_elements_PRIMAL( sz[0],sz[1], eta, chi , Mesh->P1P->sbf[j]->T[l+1]->v, IE_y, P0_y );
        ev_basis_function( length_arr, 0, phij, eta, chi, "P1" );
        n_y = normal(P0_y, Mesh->P1P->sbf[j]->T[l+1]->v);
        curl_y = curl_bf(Mesh->P1P->sbf[j]->T[l+1]->v, 0, n_y);
        res = HP_kern(kappa, weights, weights, phii, phij, P0_x, P0_y, IE_x, IE_y, n_x, n_y, curl_x, curl_y, length_arr, length_arr );
        aux = gsl_complex_add(aux, res);
        
        ev_basis_function( 1, 0, phij, eta2, chi2, "P1" );
        integration_elements_PRIMAL( sz[0],sz[1], eta, chi , Mesh->P1P->sbf[i]->T[k+1]->v, IE_y, P0_y );
        singular(P0_x, Mesh->P1P->sbf[i]->T[k+1]->v, weights2, cres, IE_x,sz2, order[2], 0, kappa, "P1", phij[0]);
        n_y = normal(P0_y, Mesh->P1P->sbf[i]->T[k+1]->v);
        curl_y = curl_bf(Mesh->P1P->sbf[i]->T[k+1]->v, 0, n_y);
        gsl_blas_ddot(n_x, n_y, &resn);
        aux = gsl_complex_add(aux, gsl_complex_rect(-resn*kappa*kappa*cres[0],-resn*kappa*kappa*cres[1]));
        singular(P0_x, Mesh->P1P->sbf[i]->T[k+1]->v, weights2, cres, IE_x, sz2, order[2], 1, kappa, "P0", 1.0);
        gsl_blas_ddot(curl_x, curl_y, &resc);
        aux = gsl_complex_add(aux, gsl_complex_rect(cres[0]*resc,cres[1]*resc));
        
        integration_elements_PRIMAL( sz[0],sz[1], eta, chi, Mesh->P1P->sbf[i]->T[k+1]->v, IE_x, P0_x );
        n_x = normal(P0_x, Mesh->P1P->sbf[i]->T[k+1]->v);
        ev_basis_function( length_arr, 0, phii, eta, chi, "P1" );
        curl_x = curl_bf(Mesh->P1P->sbf[i]->T[k+1]->v, 0, n_x);
        integration_elements_PRIMAL( sz[0],sz[1], eta, chi , Mesh->P1P->sbf[j]->T[l+1]->v, IE_y, P0_y );
        ev_basis_function( length_arr, 0, phij, eta, chi, "P1" );
        n_y = normal(P0_y, Mesh->P1P->sbf[j]->T[l+1]->v);
        curl_y = curl_bf(Mesh->P1P->sbf[j]->T[l+1]->v, 0, n_y);
        res = HP_kern(kappa, weights, weights, phii, phij, P0_x, P0_y, IE_x, IE_y, n_x, n_y, curl_x, curl_y, length_arr, length_arr );
        aux = gsl_complex_add(aux, res);
    }
    else if (opt==3)                                                                                               
    {  
        double vert1[9];
        double vert2[9];
        vert1[0] = Mesh->PD1->T[i]->T[k]->v[6];
        vert1[1] = Mesh->PD1->T[i]->T[k]->v[7]; 
        vert1[2] = Mesh->PD1->T[i]->T[k]->v[8]; 
        vert1[3] = Mesh->PD1->T[i]->T[k]->v[0]; 
        vert1[4] = Mesh->PD1->T[i]->T[k]->v[1]; 
        vert1[5] = Mesh->PD1->T[i]->T[k]->v[2]; 
        vert1[6] = Mesh->PD1->T[i]->T[k]->v[3]; 
        vert1[7] = Mesh->PD1->T[i]->T[k]->v[4]; 
        vert1[8] = Mesh->PD1->T[i]->T[k]->v[5]; 
        vert2[0] = Mesh->PD1->T[i]->T[k+1]->v[3];                                                             
        vert2[1] = Mesh->PD1->T[i]->T[k+1]->v[4];                                                             
        vert2[2] = Mesh->PD1->T[i]->T[k+1]->v[5];                                                             
        vert2[3] = Mesh->PD1->T[i]->T[k+1]->v[6];                                                             
        vert2[4] = Mesh->PD1->T[i]->T[k+1]->v[7];                                                             
        vert2[5] = Mesh->PD1->T[i]->T[k+1]->v[8];                                                             
        vert2[6] = Mesh->PD1->T[i]->T[k+1]->v[0];                                                             
        vert2[7] = Mesh->PD1->T[i]->T[k+1]->v[1];                                                             
        vert2[8] = Mesh->PD1->T[i]->T[k+1]->v[2];  

        integration_elements_PRIMAL( sz[0],sz[1], eta, chi, vert2, IE_x, P0_x );         
        n_x = normal(P0_x, vert2);                                                       
        ev_basis_function( length_arr, 0, phii, eta, chi, "P1" );                                             
        curl_x = curl_bf(vert2, 0, n_x);                                                 
        integration_elements_PRIMAL( sz[0],sz[1], eta, chi , Mesh->P1P->sbf[j]->T[l]->v, IE_y, P0_y );      
        ev_basis_function( length_arr, 0, phij, eta, chi, "P1" );                                             
        n_y = normal(P0_y, Mesh->P1P->sbf[j]->T[l]->v);                                                     
        curl_y = curl_bf(Mesh->P1P->sbf[j]->T[l]->v, 0, n_y);                                               
        res = HP_kern(kappa, weights, weights, phii, phij, P0_x, P0_y, IE_x, IE_y, n_x, n_y, curl_x, curl_y, length_arr, length_arr );
        aux = gsl_complex_add(aux, res);                                                                    


        integration_elements_PRIMAL( sz[0],sz[1], eta, chi, vert2, IE_x, P0_x );         
        n_x = normal(P0_x, vert2);                                                       
        ev_basis_function( length_arr, 0, phii, eta, chi, "P1" );                                             
        curl_x = curl_bf(vert2, 0, n_x);                                                 
        integration_elements_PRIMAL( sz[0],sz[1], eta, chi , Mesh->P1P->sbf[j]->T[l+1]->v, IE_y, P0_y );        
        ev_basis_function( length_arr, 0, phij, eta, chi, "P1" );                                             
        n_y = normal(P0_y, Mesh->P1P->sbf[j]->T[l+1]->v);                                                       
        curl_y = curl_bf(Mesh->P1P->sbf[j]->T[l+1]->v, 0, n_y);                                                 
        res = HP_kern(kappa, weights, weights, phii, phij, P0_x, P0_y, IE_x, IE_y, n_x, n_y, curl_x, curl_y, length_arr, length_arr );
        aux = gsl_complex_add(aux, res);  


        integration_elements_PRIMAL( sz[0],sz[1], eta, chi, vert1, IE_x, P0_x );         
        n_x = normal(P0_x, vert1);                                                       
        ev_basis_function( length_arr, 0, phii, eta, chi, "P1" );                                             
        curl_x = curl_bf(vert1, 0, n_x);                                                 
        integration_elements_PRIMAL( sz[0],sz[1], eta, chi , Mesh->P1P->sbf[j]->T[l]->v, IE_y, P0_y );        
        ev_basis_function( length_arr, 0, phij, eta, chi, "P1" );                                             
        n_y = normal(P0_y, Mesh->P1P->sbf[j]->T[l]->v);                                                       
        curl_y = curl_bf(Mesh->P1P->sbf[j]->T[l]->v, 0, n_y);                                                 
        res = HP_kern(kappa, weights, weights, phii, phij, P0_x, P0_y, IE_x, IE_y, n_x, n_y, curl_x, curl_y, length_arr, length_arr );
        aux = gsl_complex_add(aux, res);  

                                                                                        
                                                                            
    }
    else if(opt==4)
    {
        double vert1[9];                                                                                      
        double vert2[9];                                                                                      
        vert1[0] = Mesh->PD1->T[i]->T[k]->v[6];                                                               
        vert1[1] = Mesh->PD1->T[i]->T[k]->v[7];                                                               
        vert1[2] = Mesh->PD1->T[i]->T[k]->v[8];                                                               
        vert1[3] = Mesh->PD1->T[i]->T[k]->v[0];                                                               
        vert1[4] = Mesh->PD1->T[i]->T[k]->v[1];                                                               
        vert1[5] = Mesh->PD1->T[i]->T[k]->v[2];                                                               
        vert1[6] = Mesh->PD1->T[i]->T[k]->v[3];                                                               
        vert1[7] = Mesh->PD1->T[i]->T[k]->v[4];                                                               
        vert1[8] = Mesh->PD1->T[i]->T[k]->v[5];                                                               
        vert2[0] = Mesh->PD1->T[i]->T[k+1]->v[3];                                                             
        vert2[1] = Mesh->PD1->T[i]->T[k+1]->v[4];                                                             
        vert2[2] = Mesh->PD1->T[i]->T[k+1]->v[5];                                                             
        vert2[3] = Mesh->PD1->T[i]->T[k+1]->v[6];                                                             
        vert2[4] = Mesh->PD1->T[i]->T[k+1]->v[7];                                                             
        vert2[5] = Mesh->PD1->T[i]->T[k+1]->v[8];                                                             
        vert2[6] = Mesh->PD1->T[i]->T[k+1]->v[0];                                                             
        vert2[7] = Mesh->PD1->T[i]->T[k+1]->v[1];                                                             
        vert2[8] = Mesh->PD1->T[i]->T[k+1]->v[2];

        integration_elements_PRIMAL( sz[0],sz[1], eta, chi, vert1, IE_x, P0_x );                              
        n_x = normal(P0_x, vert1);                                                                            
        ev_basis_function( length_arr, 0, phii, eta, chi, "P1" );                                             
        curl_x = curl_bf(vert1, 0, n_x);                                                                      
        integration_elements_PRIMAL( sz[0],sz[1], eta, chi , Mesh->P1P->sbf[j]->T[l+1]->v, IE_y, P0_y );        
        ev_basis_function( length_arr, 0, phij, eta, chi, "P1" );                                             
        n_y = normal(P0_y, Mesh->P1P->sbf[j]->T[l+1]->v);                                                       
        curl_y = curl_bf(Mesh->P1P->sbf[j]->T[l+1]->v, 0, n_y);                                                 
        res = HP_kern(kappa, weights, weights, phii, phij, P0_x, P0_y, IE_x, IE_y, n_x, n_y, curl_x, curl_y, length_arr, length_arr );
        aux = gsl_complex_add(aux, res);                                                                      

        integration_elements_PRIMAL( sz[0],sz[1], eta, chi, vert1, IE_x, P0_x );                              
        n_x = normal(P0_x, vert1);                                                                            
        ev_basis_function( length_arr, 0, phii, eta, chi, "P1" );                                             
        curl_x = curl_bf(vert1, 0, n_x);                                                                      
        integration_elements_PRIMAL( sz[0],sz[1], eta, chi , Mesh->P1P->sbf[j]->T[l]->v, IE_y, P0_y );      
        ev_basis_function( length_arr, 0, phij, eta, chi, "P1" );                                             
        n_y = normal(P0_y, Mesh->P1P->sbf[j]->T[l]->v);                                                     
        curl_y = curl_bf(Mesh->P1P->sbf[j]->T[l]->v, 0, n_y);                                               
        res = HP_kern(kappa, weights, weights, phii, phij, P0_x, P0_y, IE_x, IE_y, n_x, n_y, curl_x, curl_y, length_arr, length_arr );
        aux = gsl_complex_add(aux, res);                                                                      


        integration_elements_PRIMAL( sz[0],sz[1], eta, chi, vert2, IE_x, P0_x );                              
        n_x = normal(P0_x, vert2);                                                                            
        ev_basis_function( length_arr, 0, phii, eta, chi, "P1" );                                             
        curl_x = curl_bf(vert2, 0, n_x);                                                                      
        integration_elements_PRIMAL( sz[0],sz[1], eta, chi , Mesh->P1P->sbf[j]->T[l+1]->v, IE_y, P0_y );        
        ev_basis_function( length_arr, 0, phij, eta, chi, "P1" );                                             
        n_y = normal(P0_y, Mesh->P1P->sbf[j]->T[l+1]->v);                                                       
        curl_y = curl_bf(Mesh->P1P->sbf[j]->T[l+1]->v, 0, n_y);                                                 
        res = HP_kern(kappa, weights, weights, phii, phij, P0_x, P0_y, IE_x, IE_y, n_x, n_y, curl_x, curl_y, length_arr, length_arr );
        aux = gsl_complex_add(aux, res); 

        
    }
    else if(opt==5)                                                                                           
    {                                                                                                         
        double vert1[9];                                                                                      
        double vert2[9];                                                                                      
        vert1[0] = Mesh->PD1->T[i]->T[k]->v[6];                                                               
        vert1[1] = Mesh->PD1->T[i]->T[k]->v[7];                                                               
        vert1[2] = Mesh->PD1->T[i]->T[k]->v[8];                                                               
        vert1[3] = Mesh->PD1->T[i]->T[k]->v[0];                                                               
        vert1[4] = Mesh->PD1->T[i]->T[k]->v[1];                                                               
        vert1[5] = Mesh->PD1->T[i]->T[k]->v[2];                                                               
        vert1[6] = Mesh->PD1->T[i]->T[k]->v[3];                                                               
        vert1[7] = Mesh->PD1->T[i]->T[k]->v[4];                                                               
        vert1[8] = Mesh->PD1->T[i]->T[k]->v[5];                                                               
        vert2[0] = Mesh->PD1->T[i]->T[k+1]->v[3];                                                             
        vert2[1] = Mesh->PD1->T[i]->T[k+1]->v[4];                                                             
        vert2[2] = Mesh->PD1->T[i]->T[k+1]->v[5];                                                             
        vert2[3] = Mesh->PD1->T[i]->T[k+1]->v[6];                                                             
        vert2[4] = Mesh->PD1->T[i]->T[k+1]->v[7];                                                             
        vert2[5] = Mesh->PD1->T[i]->T[k+1]->v[8];                                                             
        vert2[6] = Mesh->PD1->T[i]->T[k+1]->v[0];                                                             
        vert2[7] = Mesh->PD1->T[i]->T[k+1]->v[1];                                                             
        vert2[8] = Mesh->PD1->T[i]->T[k+1]->v[2];                                                              
        
        integration_elements_PRIMAL( sz[0],sz[1], eta, chi, vert1, IE_x, P0_x );                              
        n_x = normal(P0_x, vert1);                                                                            
        ev_basis_function( length_arr, 0, phii, eta, chi, "P1" );                                             
        curl_x = curl_bf(vert1, 0, n_x);                                                                      
        integration_elements_PRIMAL( sz[0],sz[1], eta, chi , Mesh->P1P->sbf[j]->T[l+1]->v, IE_y, P0_y );      
        ev_basis_function( length_arr, 0, phij, eta, chi, "P1" );                                             
        n_y = normal(P0_y, Mesh->P1P->sbf[j]->T[l+1]->v);                                                     
        curl_y = curl_bf(Mesh->P1P->sbf[j]->T[l+1]->v, 0, n_y);                                               
        res = HP_kern(kappa, weights, weights, phii, phij, P0_x, P0_y, IE_x, IE_y, n_x, n_y, curl_x, curl_y, length_arr, length_arr );
        aux = gsl_complex_add(aux, res);                                                                      

        integration_elements_PRIMAL( sz[0],sz[1], eta, chi, vert1, IE_x, P0_x );                              
        n_x = normal(P0_x, vert1);                                                                            
        ev_basis_function( length_arr, 0, phii, eta, chi, "P1" );                                             
        curl_x = curl_bf(vert1, 0, n_x);                                                                      
        integration_elements_PRIMAL( sz[0],sz[1], eta, chi , Mesh->P1P->sbf[j]->T[l]->v, IE_y, P0_y );        
        ev_basis_function( length_arr, 0, phij, eta, chi, "P1" );                                             
        n_y = normal(P0_y, Mesh->P1P->sbf[j]->T[l]->v);                                                       
        curl_y = curl_bf(Mesh->P1P->sbf[j]->T[l]->v, 0, n_y);                                                 
        res = HP_kern(kappa, weights, weights, phii, phij, P0_x, P0_y, IE_x, IE_y, n_x, n_y, curl_x, curl_y, length_arr, length_arr );
        aux = gsl_complex_add(aux, res);                                                                      

        integration_elements_PRIMAL( sz[0],sz[1], eta, chi, vert2, IE_x, P0_x );                              
        n_x = normal(P0_x, vert2);                                                                            
        ev_basis_function( length_arr, 0, phii, eta, chi, "P1" );                                             
        curl_x = curl_bf(vert2, 0, n_x);                                                                      
        integration_elements_PRIMAL( sz[0],sz[1], eta, chi , Mesh->P1P->sbf[j]->T[l+1]->v, IE_y, P0_y );      
        ev_basis_function( length_arr, 0, phij, eta, chi, "P1" );                                             
        n_y = normal(P0_y, Mesh->P1P->sbf[j]->T[l+1]->v);                                                     
        curl_y = curl_bf(Mesh->P1P->sbf[j]->T[l+1]->v, 0, n_y);                                               
        res = HP_kern(kappa, weights, weights, phii, phij, P0_x, P0_y, IE_x, IE_y, n_x, n_y, curl_x, curl_y, length_arr, length_arr );
        aux = gsl_complex_add(aux, res); 


        integration_elements_PRIMAL( sz[0],sz[1], eta, chi, vert2, IE_x, P0_x );                              
        n_x = normal(P0_x, vert2);                                                                            
        ev_basis_function( length_arr, 0, phii, eta, chi, "P1" );                                             
        curl_x = curl_bf(vert2, 0, n_x);                                                                      
        integration_elements_PRIMAL( sz[0],sz[1], eta, chi , Mesh->P1P->sbf[j]->T[l]->v, IE_y, P0_y );      
        ev_basis_function( length_arr, 0, phij, eta, chi, "P1" );                                             
        n_y = normal(P0_y, Mesh->P1P->sbf[j]->T[l]->v);                                                     
        curl_y = curl_bf(Mesh->P1P->sbf[j]->T[l]->v, 0, n_y);                                               
        res = HP_kern(kappa, weights, weights, phii, phij, P0_x, P0_y, IE_x, IE_y, n_x, n_y, curl_x, curl_y, length_arr, length_arr );
        aux = gsl_complex_add(aux, res);  

    }
    else if(opt==6)
    {
        double vert1[9];                                                                                                                                                                         
        double vert2[9];
        vert1[0] = Mesh->PD1->T[i]->T[k]->v[6];                                                               
        vert1[1] = Mesh->PD1->T[i]->T[k]->v[7];                                                               
        vert1[2] = Mesh->PD1->T[i]->T[k]->v[8];                                                               
        vert1[3] = Mesh->PD1->T[i]->T[k]->v[0];                                                               
        vert1[4] = Mesh->PD1->T[i]->T[k]->v[1];                                                               
        vert1[5] = Mesh->PD1->T[i]->T[k]->v[2];                                                               
        vert1[6] = Mesh->PD1->T[i]->T[k]->v[3];                                                               
        vert1[7] = Mesh->PD1->T[i]->T[k]->v[4];                                                               
        vert1[8] = Mesh->PD1->T[i]->T[k]->v[5];                                                               
        vert2[0] = Mesh->PD1->T[i]->T[k+1]->v[3];                                                             
        vert2[1] = Mesh->PD1->T[i]->T[k+1]->v[4];                                                             
        vert2[2] = Mesh->PD1->T[i]->T[k+1]->v[5];                                                             
        vert2[3] = Mesh->PD1->T[i]->T[k+1]->v[6];                                                             
        vert2[4] = Mesh->PD1->T[i]->T[k+1]->v[7];                                                             
        vert2[5] = Mesh->PD1->T[i]->T[k+1]->v[8];                                                             
        vert2[6] = Mesh->PD1->T[i]->T[k+1]->v[0];                                                             
        vert2[7] = Mesh->PD1->T[i]->T[k+1]->v[1];                                                             
        vert2[8] = Mesh->PD1->T[i]->T[k+1]->v[2];                                                             
                                                                                                              
        integration_elements_PRIMAL( sz[0],sz[1], eta, chi, vert1, IE_x, P0_x );                              
        n_x = normal(P0_x, vert1);                                                                            
        ev_basis_function( length_arr, 0, phii, eta, chi, "P1" );                                             
        curl_x = curl_bf(vert1, 0, n_x);                                                                      
        integration_elements_PRIMAL( sz2[0],sz2[1], eta2, chi2, vert1, IE_x, P0_x );                          
        ev_basis_function( 1, 0, phij, eta2, chi2, "P1" );                                                    
        integration_elements_PRIMAL( sz[0],sz[1], eta, chi , Mesh->P1P->sbf[j]->T[l+1]->v, IE_y, P0_y );      
        singular(P0_x, Mesh->P1P->sbf[j]->T[l+1]->v, weights2, cres, IE_x,sz2, order[2], 0, kappa, "P1", phij[0]);      
        n_y = normal(P0_y, Mesh->P1P->sbf[j]->T[l+1]->v);                                                     
        curl_y = curl_bf(Mesh->P1P->sbf[j]->T[l+1]->v, 0, n_y);                                               
        gsl_blas_ddot(n_x, n_y, &resn);                                                                       
        aux = gsl_complex_add(aux, gsl_complex_rect(-resn*kappa*kappa*cres[0],-resn*kappa*kappa*cres[1]));    
        singular(P0_x, Mesh->P1P->sbf[j]->T[l+1]->v, weights2, cres, IE_x, sz2, order[2], 1, kappa, "P0", 1.0);         
        gsl_blas_ddot(curl_x, curl_y, &resc);                                                                 
        aux = gsl_complex_add(aux, gsl_complex_rect(cres[0]*resc,cres[1]*resc)); 
    }
    else if(opt==7)
    {

        double vert1[9];
        double vert2[9];                                                                                      
        vert1[0] = Mesh->PD1->T[i]->T[k]->v[6];                                                               
        vert1[1] = Mesh->PD1->T[i]->T[k]->v[7];                                                               
        vert1[2] = Mesh->PD1->T[i]->T[k]->v[8];                                                               
        vert1[3] = Mesh->PD1->T[i]->T[k]->v[0];                                                               
        vert1[4] = Mesh->PD1->T[i]->T[k]->v[1];                                                               
        vert1[5] = Mesh->PD1->T[i]->T[k]->v[2];                                                               
        vert1[6] = Mesh->PD1->T[i]->T[k]->v[3];                                                               
        vert1[7] = Mesh->PD1->T[i]->T[k]->v[4];                                                               
        vert1[8] = Mesh->PD1->T[i]->T[k]->v[5];                                                               
        vert2[0] = Mesh->PD1->T[i]->T[k+1]->v[3];                                                             
        vert2[1] = Mesh->PD1->T[i]->T[k+1]->v[4];                                                             
        vert2[2] = Mesh->PD1->T[i]->T[k+1]->v[5];                                                             
        vert2[3] = Mesh->PD1->T[i]->T[k+1]->v[6];                                                             
        vert2[4] = Mesh->PD1->T[i]->T[k+1]->v[7];                                                             
        vert2[5] = Mesh->PD1->T[i]->T[k+1]->v[8];                                                             
        vert2[6] = Mesh->PD1->T[i]->T[k+1]->v[0];                                                             
        vert2[7] = Mesh->PD1->T[i]->T[k+1]->v[1];                                                             
        vert2[8] = Mesh->PD1->T[i]->T[k+1]->v[2];                                                             
                                                                                                              
        integration_elements_PRIMAL( sz[0],sz[1], eta, chi, vert2, IE_x, P0_x );                              
        n_x = normal(P0_x, vert2);                                                                            
        ev_basis_function( length_arr, 0, phii, eta, chi, "P1" );                                             
        curl_x = curl_bf(vert2, 0, n_x);                                                                      
        integration_elements_PRIMAL( sz2[0],sz2[1], eta2, chi2, vert2, IE_x, P0_x );                          
        ev_basis_function( 1, 0, phij, eta2, chi2, "P1" );                                                    
        integration_elements_PRIMAL( sz[0],sz[1], eta, chi , Mesh->P1P->sbf[j]->T[l]->v, IE_y, P0_y );        
        singular(P0_x, Mesh->P1P->sbf[j]->T[l]->v, weights2, cres, IE_x,sz2, order[2], 0, kappa, "P1", phij[0]);        
        n_y = normal(P0_y, Mesh->P1P->sbf[j]->T[l]->v);                                                       
        curl_y = curl_bf(Mesh->P1P->sbf[j]->T[l]->v, 0, n_y);                                                 
        gsl_blas_ddot(n_x, n_y, &resn);                                                                       
        aux = gsl_complex_add(aux, gsl_complex_rect(-resn*kappa*kappa*cres[0],-resn*kappa*kappa*cres[1]));    
        singular(P0_x, Mesh->P1P->sbf[j]->T[l]->v, weights2, cres, IE_x, sz2, order[2], 1, kappa, "P0", 1.0);           
        gsl_blas_ddot(curl_x, curl_y, &resc);                                                                 
        aux = gsl_complex_add(aux, gsl_complex_rect(cres[0]*resc,cres[1]*resc));  
    }
/*    else if(opt==8)
    {

        integration_elements_PRIMAL( sz[0],sz[1], eta, chi, Mesh->P1P->sbf[i]->T[k+1]->v, IE_x, P0_x );                              
        n_x = normal(P0_x, Mesh->P1P->sbf[i]->T[k+1]->v);                                                                            
        ev_basis_function( length_arr, 0, phii, eta, chi, "P1" );                                             
        curl_x = curl_bf(Mesh->P1P->sbf[i]->T[k+1]->v, 0, n_x);                                                                      
        integration_elements_PRIMAL( sz[0],sz[1], eta, chi , Mesh->P1P->sbf[j]->T[l]->v, IE_y, P0_y );      
        ev_basis_function( length_arr, 0, phij, eta, chi, "P1" );                                             
        n_y = normal(P0_y, Mesh->P1P->sbf[j]->T[l]->v);                                                     
        curl_y = curl_bf(Mesh->P1P->sbf[j]->T[l]->v, 0, n_y);                                               
        res = HP_kern(kappa, weights, weights, phii, phij, P0_x, P0_y, IE_x, IE_y, n_x, n_y, curl_x, curl_y, length_arr, length_arr );
        aux = gsl_complex_add(aux, res); 
        
        integration_elements_PRIMAL( sz[0],sz[1], eta, chi, Mesh->P1P->sbf[i]->T[k]->v, IE_x, P0_x );                              
        n_x = normal(P0_x,Mesh->P1P->sbf[i]->T[k]->v);                                                                            
        ev_basis_function( length_arr, 0, phii, eta, chi, "P1" );                                             
        curl_x = curl_bf(Mesh->P1P->sbf[i]->T[k]->v, 0, n_x);                                                                      
        integration_elements_PRIMAL( sz2[0],sz2[1], eta2, chi2, Mesh->P1P->sbf[i]->T[k]->v, IE_x, P0_x );                          
        ev_basis_function( 1, 0, phij, eta2, chi2, "P1" );                                                    
        integration_elements_PRIMAL( sz[0],sz[1], eta, chi , Mesh->P1P->sbf[j]->T[l]->v, IE_y, P0_y );      
        singular(P0_x, Mesh->P1P->sbf[j]->T[l]->v, cres, IE_x,sz2, order[2], 0, kappa, "P1", phij[0]);      
        n_y = normal(P0_y, Mesh->P1P->sbf[j]->T[l]->v);                                                     
        curl_y = curl_bf(Mesh->P1P->sbf[j]->T[l]->v, 0, n_y);                                               
        gsl_blas_ddot(n_x, n_y, &resn);                                                                       
        aux = gsl_complex_add(aux, gsl_complex_rect(-resn*kappa*kappa*cres[0],-resn*kappa*kappa*cres[1]));    
        singular(P0_x, Mesh->P1P->sbf[j]->T[l]->v, cres, IE_x, sz2, order[2], 1, kappa, "P0", 1.0);         
        gsl_blas_ddot(curl_x, curl_y, &resc);                                                                 
        aux = gsl_complex_add(aux, gsl_complex_rect(cres[0]*resc,cres[1]*resc));   
    }
    else if(opt==9)                                                                                           
    {                                                                                                         
                                                                                                              
        integration_elements_PRIMAL( sz[0],sz[1], eta, chi, Mesh->P1P->sbf[i]->T[k]->v, IE_x, P0_x );       
        n_x = normal(P0_x, Mesh->P1P->sbf[i]->T[k]->v);                                                     
        ev_basis_function( length_arr, 0, phii, eta, chi, "P1" );                                             
        curl_x = curl_bf(Mesh->P1P->sbf[i]->T[k]->v, 0, n_x);                                               
        integration_elements_PRIMAL( sz[0],sz[1], eta, chi , Mesh->P1P->sbf[j]->T[l]->v, IE_y, P0_y );        
        ev_basis_function( length_arr, 0, phij, eta, chi, "P1" );                                             
        n_y = normal(P0_y, Mesh->P1P->sbf[j]->T[l]->v);                                                       
        curl_y = curl_bf(Mesh->P1P->sbf[j]->T[l]->v, 0, n_y);                                                 
        res = HP_kern(kappa, weights, weights, phii, phij, P0_x, P0_y, IE_x, IE_y, n_x, n_y, curl_x, curl_y, length_arr, length_arr );
        aux = gsl_complex_add(aux, res);                                                                      
                                                                                                              
        integration_elements_PRIMAL( sz[0],sz[1], eta, chi, Mesh->P1P->sbf[i]->T[k+1]->v, IE_x, P0_x );         
        n_x = normal(P0_x,Mesh->P1P->sbf[i]->T[k+1]->v);                                                        
        ev_basis_function( length_arr, 0, phii, eta, chi, "P1" );                                             
        curl_x = curl_bf(Mesh->P1P->sbf[i]->T[k+1]->v, 0, n_x);                                                 
        integration_elements_PRIMAL( sz2[0],sz2[1], eta2, chi2, Mesh->P1P->sbf[i]->T[k+1]->v, IE_x, P0_x );                          
        ev_basis_function( 1, 0, phij, eta2, chi2, "P1" );                                                    
        integration_elements_PRIMAL( sz[0],sz[1], eta, chi , Mesh->P1P->sbf[j]->T[l]->v, IE_y, P0_y );        
        singular(P0_x, Mesh->P1P->sbf[j]->T[l]->v, cres, IE_x,sz2, order[2], 0, kappa, "P1", phij[0]);        
        n_y = normal(P0_y, Mesh->P1P->sbf[j]->T[l]->v);                                                       
        curl_y = curl_bf(Mesh->P1P->sbf[j]->T[l]->v, 0, n_y);                                                 
        gsl_blas_ddot(n_x, n_y, &resn);                                                                       
        aux = gsl_complex_add(aux, gsl_complex_rect(-resn*kappa*kappa*cres[0],-resn*kappa*kappa*cres[1]));    
        singular(P0_x, Mesh->P1P->sbf[j]->T[l]->v, cres, IE_x, sz2, order[2], 1, kappa, "P0", 1.0);           
        gsl_blas_ddot(curl_x, curl_y, &resc);                                                                 
        aux = gsl_complex_add(aux, gsl_complex_rect(cres[0]*resc,cres[1]*resc));                              
    }*/
    return aux;
 
    
}

void assembly_P1d1( pMeshes Mesh, int i, int dof, int k, int *sz, double *P0, double *IE, double *eta, double *chi)
{                                                                                  
    integration_elements_DUAL( sz[0],sz[1], eta, chi,Mesh->PD1->sbfb[i]->sbf[dof]->Q[k]->v, IE, P0);
}



void assembly_P1d2(pMeshes Mesh, int i, int k, int *sz, double *P0, double *IE, double *eta, double *chi)
{                                                                               
    integration_elements_PRIMAL( sz[0],sz[1], eta, chi,Mesh->P1P->sbf[i]->T[k]->v, IE, P0);                                         
}


void assembly_P1d3(pMeshes Mesh, int i, int dof, int k, int *sz, double *P0, double *IE, double *eta, double *chi)
{                                                                                                                                                                                               
    integration_elements_PRIMAL( sz[0],sz[1], eta, chi,Mesh->PD1->sbfb[i]->sbf[dof]->T[k/2]->v, IE, P0);                                           
}

gsl_complex assemble_P1_HP_entriesd(int *order, double kappa, int i, int j, int idx, int idy, double *sing, pMeshes Mesh, int dofx, int dofy)
{
    
    double weights_x[16], weights_y[16], etax[4], etay[4], chix[4], chiy[4];
    int sz_x[2], sz_y[2], k, l;
    double phii[16] = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
    double phij[16] = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
    gsl_complex res;
    gsl_complex aux = gsl_complex_rect(0.0,0.0);
    int length_arr = order[0]*order[0];
    double P0_x[3*length_arr], P0_y[3*length_arr], IE_x[length_arr], IE_y[length_arr];

    if((dofx==0 || dofx==1 || dofx==2 ) &&(dofy==0 || dofy==1 || dofy==2 ))
    {
        int stop =0;
        quad_rules_DUAL(etax, chix, sz_x, order[0]); 
        quad_rules_DUAL(etay, chiy, sz_y, order[0]);
        weights_DUAL(weights_x, order[0], sz_x[0]);
        weights_DUAL(weights_y, order[0], sz_y[0]);
        ev_basis_function( 2*sz_x[0], 0, phii, etax, chix, "P1di" ); 
        ev_basis_function( 2*sz_y[0], 0, phij, etay, chiy, "P1di" );
        
        for( k = 0; k <Mesh->P1P->sbf[idx]->len && stop==0; k=k+2 )
        {
            assembly_P1d1(Mesh, i, dofx, k/2, sz_x, P0_x, IE_x, etax, chix);
            for(l = 0; l < Mesh->P1P->sbf[idy]->len &&stop ==0; l=l+2)
            {
                if( Mesh->P1P->sbf[idx]->num[k] != Mesh->P1P->sbf[idy]->num[l] )
                {
                    assembly_P1d1(Mesh, j, dofy, l/2, sz_y, P0_y, IE_y, etay, chiy);
                    res = HP_kern2(kappa, weights_x, weights_y, phii, phij, P0_x, P0_y, IE_x, IE_y, Mesh->PD1->sbfb[i]->sbf[dofx]->Q[k/2]->n, Mesh->PD1->sbfb[j]->sbf[dofy]->Q[l/2]->n,  Mesh->PD1->sbfb[i]->sbf[dofx]->Q[k/2]->c[order[0]-1],  Mesh->PD1->sbfb[j]->sbf[dofy]->Q[l/2]->c[order[0]-1], length_arr, length_arr );
                    aux = gsl_complex_add(aux, res);
                }
                else
                {
                    aux = assemble_P1_HP_entries(order, kappa, idx, idy, sing, Mesh);                         
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
            assembly_P1d3(Mesh, i, dofx, k, sz_x, P0_x, IE_x,etax, chix);                
            for(l = 0; l < Mesh->P1P->sbf[idy]->len && stop==0; l=l+2)                                                   
            {                                                                                             
                if( Mesh->PD1->T[i]->num[k] != Mesh->P1P->sbf[idy]->num[l+1] && Mesh->PD1->T[i]->num[k+1] != Mesh->P1P->sbf[idy]->num[l] )
                {                                                                                             
                    assembly_P1d1(Mesh, j, dofy, l/2, sz_y, P0_y, IE_y, etay, chiy);      
                    res = HP_kern3(kappa, weights_x, weights_y, phii, phij, P0_x, P0_y, IE_x, IE_y, Mesh->PD1->sbfb[i]->sbf[dofx]->T[k/2]->n, Mesh->PD1->sbfb[j]->sbf[dofy]->Q[l/2]->n, Mesh->PD1->sbfb[i]->sbf[dofx]->T[k/2]->c, Mesh->PD1->sbfb[j]->sbf[dofy]->Q[l/2]->c[order[0]-1], sz_x[0], length_arr);
                    aux = gsl_complex_add(aux, res);                                                          
                }
                else
                {                                                                                             
                    aux = assemble_P1_HP_entries(order, kappa, idx, idy, sing, Mesh);                         
                    stop = 1;                                                                                 
                } 
            }
        }

    }
    else if((dofx==4 || dofx==5 || dofx==6)&&(dofy==0 || dofy==1 || dofy==2))
    {

        int stop = 0;
        quad_rules_PRIMAL(etax, chix, sz_x, order[0]); 
        quad_rules_DUAL(etay, chiy, sz_y, order[0]); 
        weights_PRIMAL(weights_x, order[0], sz_x[0]);
        weights_DUAL(weights_y, order[0], sz_y[0]); 
        ev_basis_function( sz_x[0], 0, phii, etax, chix, "P1" );  
        ev_basis_function( 2*sz_y[0], 0, phij, etay, chiy, "P1di" ); 

        for( k = 0; k <Mesh->P1P->sbf[idx]->len && stop==0; k=k+1 )                                           
        {                                                                                                     
            assembly_P1d2(Mesh, idx, k, sz_x, P0_x, IE_x, etax, chix);                                                                                                      
            for(l = 0; l < Mesh->P1P->sbf[idy]->len && stop==0; l=l+2)                                        
            {                                                                                                 
                if( Mesh->P1P->sbf[idx]->num[k] != Mesh->P1P->sbf[idy]->num[l] && Mesh->P1P->sbf[idx]->num[k] != Mesh->P1P->sbf[idy]->num[l+1] )
                {                                                                                             
                    assembly_P1d1(Mesh, j, dofy, l/2, sz_y, P0_y, IE_y, etay, chiy);      
                    res = HP_kern3(kappa, weights_x, weights_y, phii, phij, P0_x, P0_y, IE_x, IE_y, Mesh->P1P->sbf[idx]->T[k]->n, Mesh->PD1->sbfb[j]->sbf[dofy]->Q[l/2]->n, Mesh->P1P->sbf[idx]->T[k]->c, Mesh->PD1->sbfb[j]->sbf[dofy]->Q[l/2]->c[order[0]-1], sz_x[0], length_arr );
                    aux = gsl_complex_add(aux, res);                                                          
                }                                                                                             
                else if(Mesh->P1P->sbf[idx]->num[k] == Mesh->P1P->sbf[idy]->num[l] || Mesh->P1P->sbf[idx]->num[k] == Mesh->P1P->sbf[idy]->num[l+1] )
                {                                                                                             
                    aux = assemble_P1_HP_entries(order, kappa, idx, idy, sing, Mesh);                         
                    stop = 1;                                                                                 
                }    
            }                                                                                                
        }

    }
    else if((dofx==0 || dofx==1 || dofx==2) &&dofy==3)                                     
    {                                                                                                         
        int stop = 0;
        quad_rules_DUAL(etax, chix, sz_x, order[0]);                                                            
        quad_rules_PRIMAL(etay, chiy, sz_y, order[0]);
        weights_DUAL(weights_x, order[0], sz_x[0]); 
        weights_PRIMAL(weights_y, order[0], sz_y[0]);
        ev_basis_function( 2*sz_x[0], 0, phii, etax, chix, "P1di" ); 
        ev_basis_function( sz_y[0], 0, phij, etay, chiy, "P1" ); 
        for( k = 0; k <Mesh->P1P->sbf[idx]->len && stop==0; k=k+2 )                                                                         
        {   
            
            assembly_P1d1(Mesh, i, dofx, k/2, sz_x, P0_x, IE_x, etax, chix);  
            for(l = 0; l < 6 && stop==0; l=l+2)                                               
            {                                                                                             
                if( Mesh->PD1->T[j]->num[l] != Mesh->P1P->sbf[idx]->num[k+1] && Mesh->PD1->T[j]->num[l+1] != Mesh->P1P->sbf[idx]->num[k] )
                {                                                                                         
                    assembly_P1d3(Mesh, j, dofy, l, sz_y, P0_y, IE_y, etay, chiy);  
                    res = HP_kern4(kappa, weights_x, weights_y, phii, phij, P0_x, P0_y, IE_x, IE_y, Mesh->PD1->sbfb[i]->sbf[dofx]->Q[k/2]->n, Mesh->PD1->sbfb[j]->sbf[dofy]->T[l/2]->n, Mesh->PD1->sbfb[i]->sbf[dofx]->Q[k/2]->c[order[0]-1], Mesh->PD1->sbfb[j]->sbf[dofy]->T[l/2]->c, length_arr, sz_y[0]);
                    aux = gsl_complex_add(aux, res);                                                      
                }
                 else 
                {                                                                                             
                    aux = assemble_P1_HP_entries(order, kappa, idx, idy, sing, Mesh);                         
                    stop = 1;                                                                                 
                                                                                                              
                } 
            }
        }

    }
    else if(dofx==3 &&dofy==3)                                                        
    {                                      
        quad_rules_PRIMAL(etax, chix, sz_x,order[0]);
        quad_rules_PRIMAL(etay, chiy, sz_y,order[0]); 
        weights_PRIMAL(weights_x, order[0], sz_x[0]);
        weights_PRIMAL(weights_y, order[0], sz_y[0]); 
        ev_basis_function( sz_x[0], 0, phii, etax, chix, "P1" );  
        ev_basis_function( sz_y[0], 0, phij, etay, chiy, "P1" ); 
        if(i!=j)
        {
            for( k = 0; k <6; k=k+2 )                                                        
            {    
                assembly_P1d3(Mesh, i, dofx, k, sz_x,P0_x, IE_x, etax, chix);                
                for(l = 0; l < 6; l=l+2)                                                     
                {                                                                                                 
                    assembly_P1d3(Mesh, j, dofy, l, sz_y, P0_y, IE_y, etay, chiy);        
                    res = HP_kern(kappa, weights_x, weights_y, phii, phij, P0_x, P0_y, IE_x, IE_y, Mesh->PD1->sbfb[i]->sbf[dofx]->T[k/2]->n, Mesh->PD1->sbfb[j]->sbf[dofy]->T[l/2]->n, Mesh->PD1->sbfb[i]->sbf[dofx]->T[k/2]->c, Mesh->PD1->sbfb[j]->sbf[dofy]->T[l/2]->c, sz_x[0], sz_y[0] );
                    aux = gsl_complex_add(aux, res);                                                       
                }                                                                                                 
            }
        }
        else
        {
            res = assemble_P1_HP_entries(order, kappa, idx, idy, sing, Mesh); 
            aux = gsl_complex_add(aux, res);
        }
    }
    else if((dofx==4 || dofx==5 || dofx==6)&&dofy==3)                                 
    {
        int stop = 0;
        quad_rules_PRIMAL(etax, chix, sz_x,order[0]);
        quad_rules_PRIMAL(etay, chiy, sz_y,order[0]); 
        weights_PRIMAL(weights_x, order[0], sz_x[0]);                                                         
        weights_PRIMAL(weights_y, order[0], sz_y[0]);
        ev_basis_function( sz_x[0], 0, phii, etax, chix, "P1" );
        ev_basis_function( sz_y[0], 0, phij, etay, chiy, "P1" );  
        for( k = 0; k<Mesh->P1P->sbf[idx]->len && stop==0; k=k+1 )                                                                             
        {        
            assembly_P1d2(Mesh, idx, k, sz_x, P0_x, IE_x, etax, chix);                                                                                              
            
            for(l = 0; l < 6 && stop==0; l=l+2)                                                   
            { 
                if( Mesh->P1P->sbf[idx]->num[k] != Mesh->PD1->T[j]->num[l+1] && Mesh->P1P->sbf[idx]->num[k] != Mesh->PD1->T[j]->num[l] )
                {
                    assembly_P1d3(Mesh, j, dofy, l, sz_y, P0_y, IE_y, etay, chiy); 
                    res = HP_kern(kappa, weights_x, weights_y, phii, phij, P0_x, P0_y, IE_x, IE_y, Mesh->P1P->sbf[idx]->T[k]->n, Mesh->PD1->sbfb[j]->sbf[dofy]->T[l/2]->n, Mesh->P1P->sbf[idx]->T[k]->c, Mesh->PD1->sbfb[j]->sbf[dofy]->T[l/2]->c, sz_x[0], sz_y[0] );
                    aux = gsl_complex_add(aux, res); 
                }
                else if( Mesh->P1P->sbf[idx]->num[k] == Mesh->PD1->T[j]->num[l] || Mesh->P1P->sbf[idx]->num[k] == Mesh->PD1->T[j]->num[l+1])                                   
                {   
                    aux = assemble_P1_HP_entries(order, kappa, idx, idy, sing, Mesh);                         
                    stop = 1;
    
                } 
            }
        }
    }
    else if((dofx==0 || dofx==1 || dofx==2) &&(dofy==4 || dofy==5 || dofy==6))
    {
        int stop = 0;
        quad_rules_DUAL(etax, chix, sz_x, order[0]);
        quad_rules_PRIMAL(etay,chiy,sz_y,order[0]); 
        weights_DUAL(weights_x, order[0], sz_x[0]); 
        weights_PRIMAL(weights_y, order[0], sz_y[0]); 
        ev_basis_function( 2*sz_x[0], 0, phii, etax, chix, "P1di" ); 
        ev_basis_function( sz_y[0], 0, phij, etay, chiy, "P1" ); 
        for( k = 0; k <Mesh->P1P->sbf[idx]->len && stop==0; k=k+2 )
        {
        
            assembly_P1d1(Mesh, i, dofx, k/2, sz_x, P0_x, IE_x, etax, chix);
            for(l = 0; l < Mesh->P1P->sbf[idy]->len && stop==0; l=l+1)
            {
                if( Mesh->P1P->sbf[idx]->num[k] != Mesh->P1P->sbf[idy]->num[l] && Mesh->P1P->sbf[idx]->num[k+1] != Mesh->P1P->sbf[idy]->num[l] )
                {
                    assembly_P1d2(Mesh, idy, l, sz_y, P0_y, IE_y, etay, chiy);
                    res = HP_kern4(kappa, weights_x, weights_y, phii, phij, P0_x, P0_y, IE_x, IE_y, Mesh->PD1->sbfb[i]->sbf[dofx]->Q[k/2]->n, Mesh->P1P->sbf[idy]->T[l]->n, Mesh->PD1->sbfb[i]->sbf[dofx]->Q[k/2]->c[order[0]-1], Mesh->P1P->sbf[idy]->T[l]->c, length_arr, sz_y[0] );
                    aux = gsl_complex_add(aux, res);
                }
                else if(Mesh->P1P->sbf[idx]->num[k] == Mesh->P1P->sbf[idy]->num[l] || Mesh->P1P->sbf[idx]->num[k+1] == Mesh->P1P->sbf[idy]->num[l] )
                {
                    aux = assemble_P1_HP_entries(order, kappa, idx, idy, sing, Mesh);                         
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
            assembly_P1d3(Mesh, i, dofx, k, sz_x, P0_x, IE_x, etax, chix);              
                                                                                                              
            for(l = 0; l < Mesh->P1P->sbf[idy]->len && stop==0; l=l+1)                                                                          
            {                                                                                                 
                if( Mesh->P1P->sbf[idy]->num[l] != Mesh->PD1->T[i]->num[k+1] && Mesh->P1P->sbf[idy]->num[l] != Mesh->PD1->T[i]->num[k] )
                {                                                                                             
                    assembly_P1d2(Mesh, idy, l, sz_y, P0_y, IE_y, etay, chiy);        
                    res = HP_kern(kappa, weights_x, weights_y, phii, phij, P0_x, P0_y, IE_x, IE_y, Mesh->PD1->sbfb[i]->sbf[dofx]->T[k/2]->n, Mesh->P1P->sbf[idy]->T[l]->n, Mesh->PD1->sbfb[i]->sbf[dofx]->T[k/2]->c, Mesh->P1P->sbf[idy]->T[l]->c, sz_x[0], sz_y[0] );
                    aux = gsl_complex_add(aux, res);                                                          
                }                                                                                             
                else if(Mesh->P1P->sbf[idy]->num[l] == Mesh->PD1->T[i]->num[k+1] || Mesh->P1P->sbf[idy]->num[l] == Mesh->PD1->T[i]->num[k])                              
                {       
                    aux = assemble_P1_HP_entries(order, kappa, idx, idy, sing, Mesh);                                 
                    stop = 1; 
                }                                                                                             
            }                                                                                                 
        }                                                                                                     
    }
    else if((dofx==4 || dofx==5 || dofx==6)&&(dofy==4 || dofy==5 || dofy==6))
    { 
        res = assemble_P1_HP_entries(order, kappa, idx, idy, sing, Mesh);
        aux = gsl_complex_add(aux, res);
    }
    
    return aux;
}


void assemble_IDP(
                  // Input:
                  // V: vertices of each element
                  pMeshes Mesh,
                  // wavenumber
                  // Order of quadrature
                  int order,
                  // Output: complex and real part
                  double* M,
                  char *trialbf,
                  char *testbf
                  ){
    
    double phii[16] = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
    double phij[16] = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
    double weights[16], eta[4], chi[4],  di[16];
    double P0[3*order*order], IE[order*order];
    int i, j, k, length_arr, sz[2];
    
    quad_rules_PRIMAL(eta, chi, sz, order);
    weights_PRIMAL(weights, order, sz[0]);
    length_arr = sz[0];
    
    if ( (strcmp(trialbf,"P1") == 0||strcmp(trialbf,"P1b") == 0) && (strcmp(testbf,"P0") == 0||strcmp(testbf,"P0b") == 0) )
    {   
        int nd;
        if(strcmp(trialbf,"P1") == 0)
        {
            nd = Mesh->ntp;
        }
        else
        {
            nd = Mesh->ndb;
        }
        for( i = 0; i < Mesh->P1P->ndof; i++ )
        {
            for( k = 0; k < Mesh->P1P->sbf[i]->len; k++ )
            {
        
                integration_elements_PRIMAL( sz[0],sz[1], eta, chi , Mesh->P1P->sbf[i]->T[k]->v, IE, P0 );
                ev_basis_function( length_arr, 0, phii, eta, chi, trialbf );
                vec_prod( length_arr, phii, weights, IE, di );
                M[ i*nd + Mesh->P1P->sbf[i]->num[k] ]+= sum_v( di, length_arr );
                    
            }
        }
    }
    if ( strcmp(trialbf,"P1") == 0 && strcmp(testbf,"P1") == 0 )
    {

        for( i = 0; i < Mesh->P1P->ndof; i++ )
        {
            for( j = 0; j < Mesh->P1P->sbf[i]->len; j++ )
            {   
                integration_elements_PRIMAL( sz[0], sz[1], eta, chi, Mesh->P1P->sbf[i]->T[j]->v, IE, P0 );
                ev_basis_function( length_arr, 0, phii, eta, chi, testbf );
                vec_prod( length_arr, phii, weights, IE, di ); 
                for( k = 0; k < 3; k++ )
                { 
                    ev_basis_function( length_arr,k, phij, eta, chi, trialbf ); 
                    M[ i*Mesh->P1P->ndof + Mesh->P1P->sbf[i]->T[j]->dofs[k] ]+= dot_prod( length_arr, 0, 0, di, phij );
                }
            }
        }
    }

    if ( strcmp(trialbf,"P0") == 0 && strcmp(testbf,"P0") == 0 )
    { 
        for( i = 0; i < Mesh->P0P->nt; i++ )
        {

            integration_elements_PRIMAL( sz[0], sz[1], eta, chi, Mesh->P0P->sbf[i]->T->v, IE, P0 );
            vec_prod( length_arr, phii, weights, IE, di );
            M[ i*Mesh->P0P->nt + i ]+= sum_v( di, length_arr );

        }
    }
    
}

void assemble_IDD(
                  // Input:
                  pMeshes Mesh,
                  // Order of quadrature
                  int order,
                  // Output: complex and real part
                  double* M,
                  char *trialbf,
                  char *testbf
                  ){
    
    double weights[16], eta[4], chi[4], phii[16], phij[16], di[16];
    double P0[3*order*order], IE[order*order];
    int i, j, length_arr, sz[2];

    quad_rules_DUAL( eta, chi, sz, order );
    weights_DUAL( weights, order, sz[0] );
    length_arr = order*order;
    

    if( strcmp(testbf, "P0d" ) == 0 && strcmp(trialbf, "P0d") == 0 )
    {
        for( i = 0; i < Mesh->PD0->ndof; i++ )
        {
            for( j = 0; j < Mesh->PD0->sbf[i]->len; j++ )
            {
                integration_elements_DUAL( sz[0],sz[1], eta, chi, Mesh->PD0->sbf[i]->Q[j]->v, IE, P0 );
                ev_basis_function( length_arr, 0, phii, eta, chi, testbf );
                ev_basis_function( length_arr, 0, phij, eta, chi, trialbf );
                vec_prod( length_arr, phii, weights, IE, di );
                M[ i*Mesh->PD0->ndof + i]+= dot_prod(length_arr, 0, 0, di, phij);
            }
        }
    }
    
}


void assemble_IDM(
                  pMeshes Mesh,
                  // wavenumber
                  // Order of quadrature                                                                      
                  int order,                                                                                                                                                                            
                  // Output
                  double* M,
                  char *trialbf,
                  char *testbf
                  ){                                                                                          
                                                                                                              
    double weights[16], eta[16], chi[16], di[16], eta_p[16], chi_p[16], P0[3*order*order], IE[order*order];
    int i, j, k, length_arr, sz[2];

    double phii[16] = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
    double phij[16] = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
    
    quad_rules_DUAL( eta, chi, sz, order );
    weights_DUAL( weights,order, sz[0] );
    length_arr = order*order;                                                                                 
    
    if (strcmp(testbf,"P0")==0)
    {
        for( i = 0; i < Mesh->PD0->ndof; i++ )
        {
            for( j = 0; j < Mesh->P0P->nt; j++ )
            {
                for( k = 0; k < (Mesh->PD0->sbf[i])->len; k++ )
                {
                    if( (Mesh->PD0->sbf[i])->num[k] == j)
                    {
                        integration_elements_DUAL( sz[0],sz[1], eta, chi , ((Mesh->PD0->sbf[i])->Q[k])->v, IE, P0 );
                        ev_basis_function( length_arr, 0, phii, eta, chi, trialbf );
                        vec_prod( length_arr, phii, weights, IE, di );
                        M[ i*Mesh->P0P->nt + j ]+= sum_v( di, length_arr );
                    }
                    
                }
            }
        }
    }
    else if (strcmp(testbf,"P1")==0)
    {
        for( i = 0; i < Mesh->PD0->ndof; i++ )
        {
            for( j = 0; j < Mesh->PD0->sbf[i]->len; j++ )
            {
                integration_elements_DUAL_PRIMAL( sz[0],sz[1], Mesh->P1P->sbf[i]->T[j]->v, Mesh->PD0->sbf[i]->Q[j]->v, eta, chi ,eta_p, chi_p, IE, P0 );
                ev_basis_function( length_arr, 0, phii, eta, chi, trialbf );
                 for( k = 0; k < 3; k++ )                                                                      
                { 
                    ev_basis_function( length_arr, k, phij, eta_p, chi_p, testbf );
                    vec_prod( length_arr, phii, weights, IE, di );
                    M[ i*Mesh->P1P->ndof + Mesh->P1P->sbf[i]->T[j]->dofs[k] ]+= dot_prod(length_arr, 0, 0, di, phij);
                }
            }
        }
    }

}


void assemble_IDDP(pMeshes Mesh,                                                                               
                  // wavenumber                                                                               
                  // Order of quadrature                                                                      
                  int order,                                                                                  
                  // Output                                                                                   
                  double* M,                                                                                  
                  char *trialbf,                                                                              
                  char *testbf){

    double phii[16] = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};                                                                     
    double weights[16], eta[4], chi[4],  di[16];                                                              
    double P0[3*order*order], IE[order*order];                                                                
    int length_arr, sz[2]; 
    int num;
    int i,k, l;
    int idx[7];
    double coef[7];
    quad_rules_PRIMAL(eta, chi, sz, order);                                                                   
    weights_PRIMAL(weights, order, sz[0]);                                                                    
    length_arr = sz[0];
    
    if (strcmp(testbf,"P1d")==0) 
    {
        coef[3] = 1.0;
        coef[4] = 0.5;
        coef[5] = 0.5;
        coef[6] = 0.5;
        for( i = 0; i < Mesh->PD1->nt; i++ )                                                                  
        {
            coef[0] = 1.0/(Mesh->PD1->lensbf[i][0]);
            coef[1] = 1.0/(Mesh->PD1->lensbf[i][1]);
            coef[2] = 1.0/(Mesh->PD1->lensbf[i][2]);
            idx[0] = Mesh->P1P->globaldofs[Mesh->PD1->T[i]->dofs[0]];                                    
            idx[1] = Mesh->P1P->globaldofs[Mesh->PD1->T[i]->dofs[1]];                                    
            idx[2] = Mesh->P1P->globaldofs[Mesh->PD1->T[i]->dofs[2]];                                    
            idx[3] = Mesh->P1P->globaldofs[Mesh->PD1->T[i]->T[0]->dofs[2]];                              
            idx[4] = Mesh->P1P->globaldofs[Mesh->PD1->T[i]->T[0]->dofs[1]];                              
            idx[5] = Mesh->P1P->globaldofs[Mesh->PD1->T[i]->T[5]->dofs[2]];                              
            idx[6] = Mesh->P1P->globaldofs[Mesh->PD1->T[i]->T[3]->dofs[2]];
            for(k=0;k<7;k++)                                                                              
            {
                for( l = 0; l <Mesh->P1P->sbf[idx[k]]->len; l++ )                                                               
                {
                                                         
                    integration_elements_PRIMAL( sz[0],sz[1], eta, chi, Mesh->P1P->sbf[idx[k]]->T[l]->v, IE, P0 );
                    ev_basis_function( length_arr, 0, phii, eta, chi, "P1" );                              
                    vec_prod( length_arr, phii, weights, IE, di );
                    num = (int) (Mesh->P1P->sbf[idx[k]]->num[l]/6);
                    M[ i*Mesh->PD1->nt + num]+= coef[k]*sum_v( di, length_arr );
                }
            }
        } 
    }
}
void assemble_IDDB(
                  pMeshes Mesh,
                  // Order of quadrature
                  int order,
                  // Output
                  double* M,
                  char *trialbf,
                  char *testbf
                  ){
    
    double weights[16], eta[16], chi[16], di[16], P0[3*order*order], IE[order*order];
    int i, k, length_arr, sz[2];
    int ndofb = Mesh->ndofb;
    double phii[16] = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
    double phij[16] = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
    double aux[4];
    
    quad_rules_DUAL( eta, chi, sz, order );
    weights_DUAL( weights,order, sz[0] );
    length_arr = order*order;


    if (strcmp(testbf,"P1b")==0)
    {
        for( i = 0; i < Mesh->PD0->ndof; i++ )
        {
            for( k = 0; k < Mesh->PD0->sbf[i]->len; k++ )
            {
               integration_elements_DUAL( sz[0],sz[1], eta, chi, Mesh->PD0->sbf[i]->Q[k]->v, IE, P0 );
                ev_basis_function( length_arr, 0, phii, eta, chi, trialbf );
                ev_basis_function( length_arr, 0, phij, eta, chi, testbf );
                vec_prod( length_arr, phii, weights, IE, di );
                M[ i*(ndofb+1) + i]+= dot_prod(length_arr, 0, 0, di, phij);
            
                integration_elements_DUAL( sz[0],sz[1], eta, chi , Mesh->PD0->sbf[i]->Q[k]->v, IE, P0 );
                ev_basis_function( length_arr, 0, phij, eta, chi, "P1bi" );
                M[ i*ndofb + Mesh->PD0->sbf[i]->Q[k]->dofs[2] ]+= dot_prod(length_arr, 0, 0, di, phij);
                aux[0] = Mesh->PD0->sbf[i]->Q[k]->v[0];
                aux[1] = Mesh->PD0->sbf[i]->Q[k]->v[2];
                if(k < Mesh->PD0->sbf[i]->len-1)
                {
                    aux[2] = Mesh->PD0->sbf[i]->Q[k+1]->v[1];
                    aux[3] = Mesh->PD0->sbf[i]->Q[k+1]->v[2];
                }
                else
                {
                    aux[2] = Mesh->PD0->sbf[i]->Q[0]->v[1];
                    aux[3] = Mesh->PD0->sbf[i]->Q[0]->v[2];
                }
                integration_elements_DUAL( sz[0],sz[1], eta, chi , aux, IE, P0 );
                ev_basis_function( length_arr, 0, phii, eta, chi, trialbf );
                ev_basis_function( length_arr, 0, phij, eta, chi, "P1bi" );
                vec_prod( length_arr, phii, weights, IE, di );
                M[ i*ndofb + Mesh->PD0->sbf[i]->Q[k]->dofs[3] ]+= dot_prod(length_arr, 0, 0, di, phij);
                    
                
            }
        }
    }

    
}

void average_matrix_full(
                   pMeshes Mesh,
                   double* M,
                   int opt
                   ){
    

    int i,k;                                                                                        
    int ndofb;                                                                                                
    if(opt==0)                                                                                                
    {                                                                                                                                                                              
        ndofb = Mesh->P1P->ndof;  
        for( i = 0; i < Mesh->PD1->nt; i++ )                                                                  
        {   
          M[i*ndofb+Mesh->P1P->globaldofs[Mesh->PD1->T[i]->dofs[0]]] = 1.0/(Mesh->PD1->lensbf[i][0]);
          M[i*ndofb+Mesh->P1P->globaldofs[Mesh->PD1->T[i]->dofs[1]]] = 1.0/(Mesh->PD1->lensbf[i][1]);
          M[i*ndofb+Mesh->P1P->globaldofs[Mesh->PD1->T[i]->dofs[2]]] = 1.0/(Mesh->PD1->lensbf[i][2]);
          M[i*ndofb+Mesh->P1P->globaldofs[Mesh->PD1->T[i]->T[0]->dofs[1]]] = 1.0;
          M[i*ndofb+Mesh->P1P->globaldofs[Mesh->PD1->T[i]->T[0]->dofs[2]]] = 0.5;
          M[i*ndofb+Mesh->P1P->globaldofs[Mesh->PD1->T[i]->T[1]->dofs[1]]] = 0.5;
          M[i*ndofb+Mesh->P1P->globaldofs[Mesh->PD1->T[i]->T[3]->dofs[1]]] = 0.5;
        }                                                                                                                                  
    }
    else if(opt==1)
    {
        int num;
        for( i = 0; i < Mesh->P1P->ndof; i++ )                                                                
        {                                                                                                     
            for( k = 0; k <Mesh->P1P->sbf[i]->len; k++ )                                                      
            {

                num = (int)(Mesh->P1P->sbf[i]->num[k]/6);                                                     
                M[ num*Mesh->ndb + Mesh->P1P->sbf[i]->num[k]]=1;  
            }

        }

    }
}


void projection_matrix(
                   pMeshes Mesh,
                   // Order of quadrature
                   // Output
                   double* M,
                   int opt
                   ){
    
    int i, k, idx, idx2;
    int ndofb;
    if(opt==0)
    {
        ndofb = Mesh->nvb; 
        for( i = 0; i < Mesh->PD0->ndof; i++ )
        {
            
            for( k = 0; k < Mesh->PD0->sbf[i]->len; k++ )
            {
                idx = Mesh->P1P->sbf[i]->num[k];
                M[ idx*ndofb + Mesh->PD0->sbf[i]->Q[k]->dofs[0] ] = 1.0/Mesh->PD0->sbf[i]->len;
                idx2 = idx*ndofb + Mesh->PD0->sbf[i]->Q[k]->dofs[1];
                
                if(k == 0)
                {
                    if( Mesh->PD0->sbf[i]->Q[Mesh->PD0->sbf[i]->len-1]->dofs[3] == Mesh->PD0->sbf[i]->Q[k]->dofs[1])
                    {
                        M[idx2] = 0.5;
                    }
                    else
                    {
                        M[idx2] = 0.5;
                    }
                }
                else
                {
                    if(Mesh->PD0->sbf[i]->Q[k-1]->dofs[3] == Mesh->PD0->sbf[i]->Q[k]->dofs[1])
                    {
                        M[idx2] = 0.5;
                    }
                    else
                    {
                        M[idx2] = 0.5;
                    }
                }
                M[ idx*ndofb + Mesh->PD0->sbf[i]->Q[k]->dofs[2] ] = 1.0;
            }
        }     

    }

}
