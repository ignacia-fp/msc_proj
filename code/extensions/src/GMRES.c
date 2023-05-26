#include "GMRES.h"
/* gmres.c
 *
 * Copyright (C) 2014 Patrick Alken
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */


/*
 * The code in this module is based on the Householder GMRES
 * algorithm described in
 *
 * [1] H. F. Walker, Implementation of the GMRES method using
 *     Householder transformations, SIAM J. Sci. Stat. Comput.
 *     9(1), 1988.
 *
 * [2] Y. Saad, Iterative methods for sparse linear systems,
 *     2nd edition, SIAM, 2003.
 */

/*
 gmres_alloc()
 Allocate a GMRES workspace for solving an n-by-n system A x = b
 Inputs: n        - size of system
 krylov_m - size of Krylov subspace (ie: number of inner iterations)
 if this parameter is 0, the value GSL_MIN(n,10) is
 used
 Return: pointer to workspace
 */


gsl_complex gsl_linalg_complex_householder_transform2 (gsl_vector_complex * v)
{
    /* replace v[0:n-1] with a householder vector (v[0:n-1]) and
     coefficient tau that annihilate v[1:n-1] */
    
    const size_t n = v->size;
    
    if (n == 1)
    {
        gsl_complex tau;

        {
            GSL_REAL(tau) = 0.0;
            GSL_IMAG(tau) = 0.0;
        }
        
        return tau;
    }
    else
    {
        gsl_complex tau ;
        gsl_complex mu;
        gsl_complex aux;
        
        gsl_complex alpha = gsl_vector_complex_get (v, 0) ;
        double absa = gsl_complex_abs (alpha);;
        double vnorm = gsl_blas_dznrm2 (v);

        if(GSL_IMAG(alpha) ==0 && GSL_REAL(alpha) == 0)
        {
            GSL_REAL(mu) = 1;
            GSL_IMAG(mu) = 0;
        }
        else
        {
            mu = gsl_complex_div_real(alpha,gsl_complex_abs(alpha));
            
        }
        
        GSL_REAL(aux) = vnorm+absa;
        GSL_IMAG(aux) = 0.0;
        gsl_vector_complex_set(v,0,gsl_complex_mul(mu,aux));
        gsl_blas_zdotc(v,v,&tau);
        GSL_REAL(aux) = 2;
        GSL_IMAG(aux) = 0;
        tau = gsl_complex_mul(aux,gsl_complex_inverse(tau));

        return tau;
    }
    
}


int gsl_linalg_complex_householder_hv2 (gsl_complex tau, const gsl_vector_complex * v, gsl_vector_complex *  w)
{
 
    if (GSL_REAL(tau) == 0.0 && GSL_IMAG(tau) == 0.0)
        return GSL_SUCCESS;
    
    gsl_complex z, tz, ntz;
    gsl_blas_zdotc(v, w, &z);
    tz = gsl_complex_mul(tau, z);
    ntz = gsl_complex_negative(tz);
    gsl_blas_zaxpy(ntz, v, w);
    
    return GSL_SUCCESS;
}

void gsl_linalg_complex_givens (const gsl_complex a, const gsl_complex b, double complex *c, double complex *s)
{
    gsl_complex z;
    gsl_complex r;
    GSL_SET_COMPLEX(&z, 0, 0);
    if (gsl_complex_eq(b, z) == 1)
    {
        *c = 1;
        *s = 0;
    }
    else if (gsl_complex_eq(a, z) == 1)
    {
        *c = 0;
        z = gsl_complex_div_real(b,gsl_complex_abs(b));
        *s = GSL_REAL(z)+I*GSL_IMAG(z);
    }
    else
    {
        z = gsl_complex_conjugate(gsl_complex_mul(gsl_complex_inverse (a),b));
        GSL_SET_COMPLEX(&r,sqrt(1+GSL_REAL(gsl_complex_mul(z,gsl_complex_conjugate(z)))),0);
        r = gsl_complex_inverse(gsl_complex_mul(gsl_complex_div_real(a,gsl_complex_abs(a)),r)); 
        *c = GSL_REAL(r)+GSL_IMAG(r);
        *s = GSL_REAL(z)*(*c)+I*GSL_IMAG(z)*(*c);
    }
} /* gsl_linalg_givens() */

void gsl_linalg_complex_givens_gv (gsl_vector_complex * v, const size_t i, const size_t j, double complex a, double complex b)
{
    /* Apply rotation to vector v' = G^T v */
    gsl_complex c, s;
    GSL_SET_COMPLEX(&c, creal(a), cimag(a));
    GSL_SET_COMPLEX(&s, creal(b), cimag(b));
    gsl_complex vi = gsl_vector_complex_get (v, i);
    gsl_complex vj = gsl_vector_complex_get (v, j);
    gsl_vector_complex_set (v, i, gsl_complex_add(gsl_complex_mul(c, vi), gsl_complex_mul (s, vj)));
    gsl_vector_complex_set (v, j, gsl_complex_sub(gsl_complex_mul(gsl_complex_conjugate(c), vj), gsl_complex_mul (gsl_complex_conjugate(s), vi)));
}


void gsl_linalg_itersolve_alloc_c(const gsl_linalg_itersolve_type_c *T, gsl_linalg_itersolve_c *w, const size_t n, const size_t m)
{
    if (w == NULL)
    {
        //GSL_ERROR_NULL("failed to allocate space for itersolve struct", GSL_ENOMEM);
    }
    
    w->type = T;
    w->normr = 0.0;
    
    w->state = w->type->alloc(n, m);
    
    if (w->state == NULL)
    {
        gsl_linalg_itersolve_free_c(w);
        GSL_ERROR_NULL("failed to allocate space for itersolve state", GSL_ENOMEM);
    }

}

void gsl_linalg_itersolve_free_c(gsl_linalg_itersolve_c *w)
{
    RETURN_IF_NULL(w);
    
    if (w->state)
        w->type->free(w->state);
    
    free(w);
}

gmres_state_t_c *gmres_alloc_c(const size_t n, const size_t m)
{
    gmres_state_t_c *state;
    if (n == 0)
    {
        GSL_ERROR_NULL("matrix dimension n must be a positive integer",
                       GSL_EINVAL);
    }
    
    state = calloc(1, sizeof(gmres_state_t_c));
    if (!state)
    {
        GSL_ERROR_NULL("failed to allocate gmres state", GSL_ENOMEM);
    }
    
    state->n = n+1;
    
    /* compute size of Krylov subspace */
    if (m == 0)
        state->m = GSL_MIN(n, 10);
    else
        state->m = GSL_MIN(n, m);
    
    state->r = gsl_vector_complex_alloc(n);
    
    if (!state->r)
    {
        gmres_free_c(state);
        GSL_ERROR_NULL("failed to allocate r vector", GSL_ENOMEM);
    }
    
    state->H = gsl_matrix_complex_alloc(n+1, state->m + 1);
    if (!state->H)
    {
        gmres_free_c(state);
        GSL_ERROR_NULL("failed to allocate H matrix", GSL_ENOMEM);
    }
    
    state->tau = gsl_vector_complex_alloc(state->m + 1);
    if (!state->tau)
    {
        gmres_free_c(state);
        GSL_ERROR_NULL("failed to allocate tau vector", GSL_ENOMEM);
    }
    
    state->y = gsl_vector_complex_alloc(state->m + 1);
    if (!state->y)
    {
        gmres_free_c(state);
        GSL_ERROR_NULL("failed to allocate y vector", GSL_ENOMEM);
    }
    
    state->c = malloc(state->m * sizeof(double complex));
    state->s = malloc(state->m * sizeof(double complex));
    if (!state->c || !state->s)
    {
        gmres_free_c(state);
        GSL_ERROR_NULL("failed to allocate Givens vectors", GSL_ENOMEM);
    }
    
    state->normr = 0.0;
    
    return state;

} /* gmres_alloc() */

void gmres_free_c(gmres_state_t_c *state)
{
    
    if (state->r)
        gsl_vector_complex_free(state->r);
    
    if (state->H)
        gsl_matrix_complex_free(state->H);
    
    if (state->tau)
        gsl_vector_complex_free(state->tau);
    
    if (state->y)
        gsl_vector_complex_free(state->y);
    
    if (state->c)
        free(state->c);
    
    if (state->s)
        free(state->s);
    
    free(state);
} /* gmres_free() */

/*
 gmres_iterate()
 Solve A*x = b using GMRES algorithm
 Inputs: A    - sparse square matrix
 b    - right hand side vector
 tol  - stopping tolerance (see below)
 x    - (input/output) on input, initial estimate x_0;
 on output, solution vector
 work - workspace
 Return:
 GSL_SUCCESS if converged to solution (solution stored in x). In
 this case the following will be true:
 ||b - A*x|| <= tol * ||b||
 GSL_CONTINUE if not yet converged; in this case x contains the
 most recent solution vector and calling this function more times
 with the input x could result in convergence (ie: restarted GMRES)
 Notes:
 1) Based on algorithm 2.2 of (Walker, 1998 [1]) and algorithm 6.10 of
 (Saad, 2003 [2])
 2) On output, work->normr contains ||b - A*x||
 */

int gmres_iterate_c(const gsl_matrix_complex *A, const gsl_matrix_complex *P, int prec, const gsl_vector_complex *b, const double tol, gsl_vector_complex *x, gmres_state_t_c *state, double *res)
{
    const size_t N = A->size1+1;
    
    if (N != A->size2+1)
    {
        GSL_ERROR("matrix must be square", GSL_ENOTSQR);
    }
    else if (N != b->size+1)
    {
        GSL_ERROR("matrix does not match right hand side", GSL_EBADLEN);
    }
    else if (N != x->size+1)
    {
        GSL_ERROR("matrix does not match solution vector", GSL_EBADLEN);
    }
    else if (N != state->n)
    {
        GSL_ERROR("matrix does not match workspace ", GSL_EBADLEN);
    }
    else
    {
        gsl_complex  one;                                                                                     
        gsl_complex mone;                                                                                     
        gsl_complex zero;                                                                                     
        gsl_complex  tau;                             /* householder scalar */                                                                                                                  
        gsl_matrix_complex *H = state->H;               /* Hessenberg matrix */                               
        gsl_vector_complex *r = state->r;               /* residual vector */                                 
        gsl_vector_complex *w = state->y;               /* least squares RHS */                               
        gsl_matrix_complex_view Rm;                     /* R_m = H(1:m,2:m+1) */                              
        gsl_vector_complex_view ym;                     /* y(1:m) */                                          
        gsl_vector_complex_view h0 = gsl_matrix_complex_subcolumn(H, 0, 1, N - 1);
        GSL_SET_COMPLEX(&one,1,0);                                                                            
        GSL_SET_COMPLEX(&mone,-1,0);                                                                          
        GSL_SET_COMPLEX(&zero,0,0); 
        int status = 1;
        gsl_vector_complex *auxv2 = gsl_vector_complex_alloc(r->size);
        //gsl_vector_complex *raux = gsl_vector_complex_alloc(r->size); 
        const size_t maxit = state->m;
        double normb;
        if (prec == 1)                                                                                    
        {                                                                                                 
             gsl_blas_zgemv(CblasNoTrans, one, P, b, zero, auxv2);                                         
             normb = gsl_blas_dznrm2(auxv2);                                                        
        }
        else
        {
             normb = gsl_blas_dznrm2(b); /* ||b|| */
        }
        //const double reltol = tol * normb;      /* tol*||b|| */
        double normr;                           /* ||r|| */
        size_t m, k;
   

        /*
         * The Hessenberg matrix will have the following structure:
         *
         * H = [ ||r_0|| | v_1 v_2 ... v_m     ]
         *     [   u_1   | u_2 u_3 ... u_{m+1} ]
         *
         * where v_j are the orthonormal vectors spanning the Krylov
         * subpsace of length j + 1 and u_{j+1} are the householder
         * vectors of length n - j - 1.
         * In fact, u_{j+1} has length n - j since u_{j+1}[0] = 1,
         * but this 1 is not stored.
         */
        gsl_matrix_complex_set_zero(H);
        gsl_vector_complex_set_zero(auxv2);
        /* Step 1a: compute r = b - A*x_0 */
        gsl_vector_complex_memcpy(r, b);
        gsl_blas_zgemv(CblasNoTrans, mone, A, x, one, r);
 
        if (prec == 1)
        {
            gsl_blas_zgemv(CblasNoTrans, one, P, r, zero, auxv2);
            gsl_vector_complex_memcpy(r, auxv2); 
        }
        
        /* Step 1b */
        gsl_vector_complex_memcpy(&h0.vector, r);
        tau = gsl_linalg_complex_householder_transform2(&h0.vector);
        gsl_linalg_complex_householder_hv2(tau,&h0.vector,r);
        /* store tau_1 */
        gsl_vector_complex_set(state->tau, 0, tau);
        
        /* initialize w (stored in state->y) */
        gsl_vector_complex_set_zero(w);
        gsl_vector_complex_set(w, 0, gsl_vector_complex_get(r,0));
        
        for (m = 1; m <= maxit; ++m)
        {
            size_t j = m - 1; /* C indexing */
            double complex c, s;      /* Givens rotation */
            
            /* v_m */
            gsl_vector_complex_view vm = gsl_matrix_complex_subcolumn(H, m,0,N-1);
            /* v_m(m:end) */
            gsl_vector_complex_view vv = gsl_vector_complex_subvector(&vm.vector, j, N - (j+1));
            
            /* householder vector u_m for projection P_m */
            gsl_vector_complex_view um = gsl_matrix_complex_subcolumn(H, j, j+1, N - (j+1));
    
            /* Step 2a: form v_m = P_m e_m = e_m - tau_m w_m */
            gsl_vector_complex_set_zero(&vm.vector);
            gsl_vector_complex_memcpy(&vv.vector, &um.vector);
            tau = gsl_vector_complex_get(state->tau, j); /* tau_m */
            gsl_vector_complex_scale(&vv.vector, gsl_complex_mul(gsl_complex_mul(mone,tau),gsl_complex_conjugate(gsl_vector_complex_get(&vv.vector, 0))));
            gsl_vector_complex_set(&vv.vector, 0, gsl_complex_add(one,gsl_vector_complex_get(&vv.vector, 0)));
            /* Step 2a: v_m <- P_1 P_2 ... P_{m-1} v_m */
            
            for (k = j; k > 0 && k--; )
            {
                gsl_vector_complex_view uk = gsl_matrix_complex_subcolumn(H, k, k+1, N - (k+1));
                gsl_vector_complex_view vk = gsl_vector_complex_subvector(&vm.vector, k, N -(k+1));
                tau = gsl_vector_complex_get(state->tau, k);
                gsl_linalg_complex_householder_hv2(tau, &uk.vector, &vk.vector);
            }
            /* Step 2a: v_m <- A*v_m */
            gsl_blas_zgemv(CblasNoTrans, one, A, &vm.vector, zero, r);
            if (prec == 1)                                                                                        
            {                                                                                                     
                gsl_blas_zgemv(CblasNoTrans, one, P, r, zero, auxv2);                                             
                gsl_vector_complex_memcpy(r, auxv2);                                                              
            }  
            gsl_vector_complex_memcpy(&vm.vector, r);
            
            /* Step 2a: v_m <- P_m ... P_1 v_m */
            for (k = 0; k <= j; k++)
            {
                gsl_vector_complex_view uk = gsl_matrix_complex_subcolumn(H, k, k+1, N - (k+1));
                gsl_vector_complex_view vk = gsl_vector_complex_subvector(&vm.vector, k, N-(k+1));
                tau = gsl_vector_complex_get(state->tau, k);
                gsl_linalg_complex_householder_hv2(tau, &uk.vector, &vk.vector);
            }

            /* Steps 2c,2d: find P_{m+1} and set v_m <- P_{m+1} v_m */
            if (m < N-1)
            {
                /* householder vector u_{m+1} for projection P_{m+1} */
                gsl_vector_complex_view ump1 = gsl_matrix_complex_subcolumn(H, m, m, N - m -1);
                gsl_vector_complex *auxv = gsl_vector_complex_alloc(N-(m+1));                                 
                gsl_vector_complex_memcpy(auxv,&ump1.vector);
                tau = gsl_linalg_complex_householder_transform2(auxv);
                gsl_vector_complex_set(state->tau, j + 1, tau);
                gsl_linalg_complex_householder_hv2(tau,auxv,&ump1.vector);
                gsl_vector_complex_view auxv2 = gsl_matrix_complex_subcolumn(H, m,m+1,N-(m+1));
                gsl_vector_complex_memcpy(&auxv2.vector,auxv);
                gsl_vector_complex_free(auxv); 
            }
            
                
            /* Step 2e: v_m <- J_{m-1} ... J_1 v_m */
            for (k = 0; k < j; ++k)
            {
                gsl_linalg_complex_givens_gv(&vm.vector, k, k + 1, state->c[k], state->s[k]);
            }
            
            
            if (m < N-1)
            {
                /* Step 2g: find givens rotation J_m for v_m(m:m+1) */
                gsl_linalg_complex_givens(gsl_vector_complex_get(&vm.vector, j), gsl_vector_complex_get(&vm.vector, j + 1), &c, &s);
                /* store givens rotation for later use */
                state->c[j] = c;
                state->s[j] = s;
                /* Step 2h: v_m <- J_m v_m */
                 gsl_linalg_complex_givens_gv(&vm.vector, j, j + 1, c, s);
                /* Step 2h: w <- J_m w */
                gsl_linalg_complex_givens_gv(w, j, j + 1, c, s);
        
            }

            /*
             * Step 2i: R_m = [ R_{m-1}, v_m ] - already taken care
             * of due to our memory storage scheme
             */
            
            /* Step 2j: check residual w_{m+1} for convergence */
            normr = gsl_complex_abs(gsl_vector_complex_get(w, j + 1))/normb;                                                                 
            res[j] = normr;
            if (normr <= tol)
            {
                /*
                 * method has converged, break out of loop to compute
                 * update to solution vector x
                 */
                break;
            }

        }
        
        /*
         * At this point, we have either converged to a solution or
         * completed all maxit iterations. In either case, compute
         * an update to the solution vector x and test again for
         * convergence.
         */
        
        /* rewind m if we exceeded maxit iterations */
        if (m > maxit)
            m--;
        
        /* Step 3a: solve triangular system R_m y_m = w, in place */
        Rm = gsl_matrix_complex_submatrix(H, 0, 1, m, m);
        ym = gsl_vector_complex_subvector(w, 0, m);
        gsl_blas_ztrsv(CblasUpper,CblasNoTrans, CblasNonUnit, &Rm.matrix, &ym.vector);
        
        /*
         * Step 3b: update solution vector x; the loop below
         * uses a different but equivalent formulation from
         * Saad, algorithm 6.10, step 14; store Krylov projection
         * V_m y_m in 'r'
         */
        gsl_vector_complex_set_zero(r);
        for (k = m; k > 0 && k--; )
        {
            gsl_complex ymk = gsl_vector_complex_get(&ym.vector, k);
            gsl_vector_complex_view uk = gsl_matrix_complex_subcolumn(H, k, k+1, N - k-1);
            gsl_vector_complex_view rk = gsl_vector_complex_subvector(r, k, N - k-1);
            /* r <- n_k e_k + r */
            gsl_vector_complex_set(r, k, gsl_complex_add(gsl_vector_complex_get(r, k),ymk));
            
            /* r <- P_k r */
            tau = gsl_vector_complex_get(state->tau, k);
            gsl_linalg_complex_householder_hv2(tau, &uk.vector, &rk.vector);
        }
        /* x <- x + V_m y_m */
        gsl_vector_complex_add(x, r);
        /* compute new residual r = b - A*x */
        gsl_vector_complex_memcpy(r, b);
        gsl_blas_zgemv(CblasNoTrans, mone, A, x, one, r);
        if (prec == 1)                                                                                    
        {                                                                                                 
             gsl_blas_zgemv(CblasNoTrans, one, P, r, zero, auxv2);                                         
             gsl_vector_complex_memcpy(r, auxv2);                                                          
        }
        normr = gsl_blas_dznrm2(r)/normb;
        
        if (normr <= tol)
            status = 1;  /* converged */
        else
            status = 0; /* not yet converged */
        
        /* store residual norm */
        state->normr = normr;
        
        return status;
    }
} /* gmres_iterate() */

int gmres_iterate_prec(const gsl_matrix_complex *A, psupermataux sa, psupermatrix P, gsl_matrix_complex *MM, const gsl_vector_complex *b, const double tol, gsl_vector_complex *x, gmres_state_t_c *state, double *res, double *tm)
{
    const size_t N = A->size1+1;
    
    if (N != A->size2+1)
    {
        GSL_ERROR("matrix must be square", GSL_ENOTSQR);
    }
    else if (N != b->size+1)
    {
        GSL_ERROR("matrix does not match right hand side", GSL_EBADLEN);
    }
    else if (N != x->size+1)
    {
        GSL_ERROR("matrix does not match solution vector", GSL_EBADLEN);
    }
    else if (N != state->n)
    {
        GSL_ERROR("matrix does not match workspace ", GSL_EBADLEN);
    }
    else
    {   
        gsl_complex  one;
        gsl_complex mone;
        gsl_complex zero;
        gsl_complex  tau;                             /* householder scalar */ 
        gsl_matrix_complex *H = state->H;               /* Hessenberg matrix */
        gsl_vector_complex *r = state->r;               /* residual vector */
        gsl_vector_complex *w = state->y;               /* least squares RHS */
        gsl_matrix_complex_view Rm;                     /* R_m = H(1:m,2:m+1) */
        gsl_vector_complex_view ym;                     /* y(1:m) */
        gsl_vector_complex_view h0 = gsl_matrix_complex_subcolumn(H, 0, 1, N - 1);
        GSL_SET_COMPLEX(&one,1,0);
        GSL_SET_COMPLEX(&mone,-1,0);
        GSL_SET_COMPLEX(&zero,0,0);
        int status = 1; 
        gsl_vector_complex *auxv2 = gsl_vector_complex_alloc(r->size);
        const size_t maxit = state->m;
        double normb;        
        gsl_vector_complex * res2 = gsl_vector_complex_alloc(auxv2->size);                                    
        gsl_vector_complex_set_zero(res2);
        /*gettimeofday(&e, NULL);                                                                               
        tm[2] += (e.tv_sec - st.tv_sec) * 1000.0;      // sec to ms                                               
        tm[2] += (e.tv_usec - st.tv_usec) / 1000.0;
        mat_vec(sa,b,res2, tm);                    
        gettimeofday(&st, NULL);
        gsl_blas_zgemv(CblasNoTrans, one, MM, res2, zero, auxv2);
        gettimeofday(&e, NULL);
        tm[3] += (e.tv_sec - st.tv_sec) * 1000.0;      // sec to ms                                               
        tm[3] += (e.tv_usec - st.tv_usec) / 1000.0;
        normb = gsl_blas_dznrm2(auxv2);*/
        normb = gsl_blas_dznrm2(b);
  //      printf("Norm:%f\n",normb);
        //const double reltol = tol * normb;      /* tol*||b|| */
        double normr;                           /* ||r|| */
        size_t m, k;

                
        /*
         * The Hessenberg matrix will have the following structure:
         *
         * H = [ ||r_0|| | v_1 v_2 ... v_m     ]
         *     [   u_1   | u_2 u_3 ... u_{m+1} ]
         *
         * where v_j are the orthonormal vectors spanning the Krylov
         * subpsace of length j + 1 and u_{j+1} are the householder
         * vectors of length n - j - 1.
         * In fact, u_{j+1} has length n - j since u_{j+1}[0] = 1,
         * but this 1 is not stored.
         */
        gsl_matrix_complex_set_zero(H);
        gsl_vector_complex_set_zero(auxv2);
        /* Step 1a: compute r = b - A*x_0 */
        gsl_vector_complex_memcpy(r, b);
        gsl_blas_zgemv(CblasNoTrans, mone, A, x, one, r);
        /*gettimeofday(&st, NULL); 
        gsl_vector_complex_set_zero(res2);  
        gettimeofday(&e, NULL);                                                                               
        tm[2] += (e.tv_sec - st.tv_sec) * 1000.0;      // sec to ms                                              
        tm[2] += (e.tv_usec - st.tv_usec) / 1000.0;
        mat_vec(sa, r,res2, tm);
        gettimeofday(&st, NULL);
        gsl_blas_zgemv(CblasNoTrans, one, MM, res2, zero, auxv2);
        gettimeofday(&e, NULL);                                                                               
        tm[3] += (e.tv_sec - st.tv_sec) * 1000.0;      // sec to ms                                               
        tm[3] += (e.tv_usec - st.tv_usec) / 1000.0;
        gettimeofday(&st, NULL);
        gsl_vector_complex_memcpy(r, auxv2);
        gettimeofday(&e, NULL); 
        tm[2] += (e.tv_sec - st.tv_sec) * 1000.0;      // sec to ms                                            
        tm[2] += (e.tv_usec - st.tv_usec) / 1000.0;
        */
        /* Step 1b */
        gsl_vector_complex_memcpy(&h0.vector, r);
        tau = gsl_linalg_complex_householder_transform2(&h0.vector);
        gsl_linalg_complex_householder_hv2(tau,&h0.vector,r);
        /* store tau_1 */
        gsl_vector_complex_set(state->tau, 0, tau);
        
        /* initialize w (stored in state->y) */
        gsl_vector_complex_set_zero(w);
        gsl_vector_complex_set(w, 0, gsl_vector_complex_get(r,0));
        
        for (m = 1; m <= maxit; ++m)
        {
            size_t j = m - 1; /* C indexing */
            double complex c, s;      /* Givens rotation */
            
            /* v_m */
            gsl_vector_complex_view vm = gsl_matrix_complex_subcolumn(H, m,0,N-1);
            /* v_m(m:end) */
            gsl_vector_complex_view vv = gsl_vector_complex_subvector(&vm.vector, j, N - (j+1));
            
            /* householder vector u_m for projection P_m */
            gsl_vector_complex_view um = gsl_matrix_complex_subcolumn(H, j, j+1, N - (j+1));
            
            /* Step 2a: form v_m = P_m e_m = e_m - tau_m w_m */
            gsl_vector_complex_set_zero(&vm.vector);
            gsl_vector_complex_memcpy(&vv.vector, &um.vector);
            tau = gsl_vector_complex_get(state->tau, j); /* tau_m */
            gsl_vector_complex_scale(&vv.vector, gsl_complex_mul(gsl_complex_mul(mone,tau),gsl_complex_conjugate(gsl_vector_complex_get(&vv.vector, 0))));
            gsl_vector_complex_set(&vv.vector, 0, gsl_complex_add(one,gsl_vector_complex_get(&vv.vector, 0)));
            /* Step 2a: v_m <- P_1 P_2 ... P_{m-1} v_m */
            
            for (k = j; k > 0 && k--; )
            {
                gsl_vector_complex_view uk = gsl_matrix_complex_subcolumn(H, k, k+1, N - (k+1));
                gsl_vector_complex_view vk = gsl_vector_complex_subvector(&vm.vector, k, N -(k+1));
                tau = gsl_vector_complex_get(state->tau, k);
                gsl_linalg_complex_householder_hv2(tau, &uk.vector, &vk.vector);
            }
            /* Step 2a: v_m <- A*v_m */
            gsl_blas_zgemv(CblasNoTrans, one, A, &vm.vector, zero, r);
          /*  gettimeofday(&st, NULL);
            gsl_vector_complex_set_zero(res2);                                                              
            gettimeofday(&e, NULL);                                                                           
            tm[2] += (e.tv_sec - st.tv_sec) * 1000.0;      // sec to ms                                              
            tm[2] += (e.tv_usec - st.tv_usec) / 1000.0;
            mat_vec(sa, r, res2, tm); 
            gettimeofday(&st, NULL);
            gsl_blas_zgemv(CblasNoTrans, one, MM, res2, zero, auxv2);
            gettimeofday(&e, NULL);                                                                               
            tm[3] += (e.tv_sec - st.tv_sec) * 1000.0;      // sec to ms                                               
            tm[3] += (e.tv_usec - st.tv_usec) / 1000.0;
            gettimeofday(&st, NULL);
            gsl_vector_complex_memcpy(r, auxv2);
            gettimeofday(&e, NULL);                                                                               
            tm[2] += (e.tv_sec - st.tv_sec) * 1000.0;      // sec to ms                                              
            tm[2] += (e.tv_usec - st.tv_usec) / 1000.0; */
            gsl_vector_complex_memcpy(&vm.vector, r);
            
            /* Step 2a: v_m <- P_m ... P_1 v_m */
            for (k = 0; k <= j; k++)
            {
                gsl_vector_complex_view uk = gsl_matrix_complex_subcolumn(H, k, k+1, N - (k+1));
                gsl_vector_complex_view vk = gsl_vector_complex_subvector(&vm.vector, k, N-(k+1));
                tau = gsl_vector_complex_get(state->tau, k);
                gsl_linalg_complex_householder_hv2(tau, &uk.vector, &vk.vector);
            }
            
            /* Steps 2c,2d: find P_{m+1} and set v_m <- P_{m+1} v_m */
            if (m < N-1)
            {
                /* householder vector u_{m+1} for projection P_{m+1} */
                gsl_vector_complex_view ump1 = gsl_matrix_complex_subcolumn(H, m, m, N - m -1);
                gsl_vector_complex *auxv = gsl_vector_complex_alloc(N-(m+1));
                gsl_vector_complex_memcpy(auxv,&ump1.vector);
                tau = gsl_linalg_complex_householder_transform2(auxv);
                gsl_vector_complex_set(state->tau, j + 1, tau);
                gsl_linalg_complex_householder_hv2(tau,auxv,&ump1.vector);
                gsl_vector_complex_view auxv2 = gsl_matrix_complex_subcolumn(H, m,m+1,N-(m+1));
                gsl_vector_complex_memcpy(&auxv2.vector,auxv);
                gsl_vector_complex_free(auxv);
            }
            
            
            /* Step 2e: v_m <- J_{m-1} ... J_1 v_m */
            for (k = 0; k < j; ++k)
            {
                gsl_linalg_complex_givens_gv(&vm.vector, k, k + 1, state->c[k], state->s[k]);
            }
            
            
            if (m < N-1)
            {
                /* Step 2g: find givens rotation J_m for v_m(m:m+1) */
                gsl_linalg_complex_givens(gsl_vector_complex_get(&vm.vector, j), gsl_vector_complex_get(&vm.vector, j + 1), &c, &s);
                /* store givens rotation for later use */
               state->c[j] = c;
               state->s[j] = s;
                /* Step 2h: v_m <- J_m v_m */
                gsl_linalg_complex_givens_gv(&vm.vector, j, j + 1, c, s);
                /* Step 2h: w <- J_m w */
                gsl_linalg_complex_givens_gv(w, j, j + 1, c, s);
                
            }
            
            /*
             * Step 2i: R_m = [ R_{m-1}, v_m ] - already taken care
             * of due to our memory storage scheme
             */
            
            /* Step 2j: check residual w_{m+1} for convergence */
            normr = gsl_complex_abs(gsl_vector_complex_get(w, j + 1))/normb;
            res[j] = normr;
            if (normr <= tol)
            {
                /*
                 * method has converged, break out of loop to compute
                 * update to solution vector x
                 */
                break;
            }
            
        }
        
        /*
         * At this point, we have either converged to a solution or
         * completed all maxit iterations. In either case, compute
         * an update to the solution vector x and test again for
         * convergence.
         */
        
        /* rewind m if we exceeded maxit iterations */
        if (m > maxit)
            m--;
        
        /* Step 3a: solve triangular system R_m y_m = w, in place */
        Rm = gsl_matrix_complex_submatrix(H, 0, 1, m, m);
        ym = gsl_vector_complex_subvector(w, 0, m);
        gsl_blas_ztrsv(CblasUpper,CblasNoTrans, CblasNonUnit, &Rm.matrix, &ym.vector);
        
        /*
         * Step 3b: update solution vector x; the loop below
         * uses a different but equivalent formulation from
         * Saad, algorithm 6.10, step 14; store Krylov projection
         * V_m y_m in 'r'
         */
        gsl_vector_complex_set_zero(r);
        for (k = m; k > 0 && k--; )
        {
            gsl_complex ymk = gsl_vector_complex_get(&ym.vector, k);
            gsl_vector_complex_view uk = gsl_matrix_complex_subcolumn(H, k, k+1, N - k-1);
            gsl_vector_complex_view rk = gsl_vector_complex_subvector(r, k, N - k-1);
            /* r <- n_k e_k + r */
            gsl_vector_complex_set(r, k, gsl_complex_add(gsl_vector_complex_get(r, k),ymk));
            
            /* r <- P_k r */
            tau = gsl_vector_complex_get(state->tau, k);
            gsl_linalg_complex_householder_hv2(tau, &uk.vector, &rk.vector);
        }
        /* x <- x + V_m y_m */
        gsl_vector_complex_add(x, r);
        /* compute new residual r = b - A*x */
        gsl_vector_complex_memcpy(r, b);
        gsl_blas_zgemv(CblasNoTrans, mone, A, x, one, r);
        /*gettimeofday(&st, NULL); 
        gsl_vector_complex_set_zero(res2); 
        gettimeofday(&e, NULL);                                                                               
        tm[2] += (e.tv_sec - st.tv_sec) * 1000.0;      // sec to ms                                              
        tm[2] += (e.tv_usec - st.tv_usec) / 1000.0;                                                                                                                                 
        mat_vec(sa, r, res2, tm);
        gettimeofday(&st, NULL);
        gsl_blas_zgemv(CblasNoTrans, one, MM, res2, zero, auxv2);   
        gettimeofday(&e, NULL);                                                                               
        tm[3] += (e.tv_sec - st.tv_sec) * 1000.0;      // sec to ms                                              
        tm[3] += (e.tv_usec - st.tv_usec) / 1000.0; 
        gettimeofday(&st, NULL);                                      
        gsl_vector_complex_memcpy(r, auxv2);
        gettimeofday(&e, NULL);                                                                               
        tm[2] += (e.tv_sec - st.tv_sec) * 1000.0;      // sec to ms                                              
        tm[2] += (e.tv_usec - st.tv_usec) / 1000.0; */
        normr = gsl_blas_dznrm2(r)/normb; 
        if (normr <= tol)
            status = 1;  /* converged */
        else
            status = 0; /* not yet converged */
        
        /* store residual norm */
        state->normr = normr;
        
        return status;
    }
} /* gmres_iterate() */

double gmres_normr_c(const gmres_state_t_c *state)
{
    return state->normr;
} /* gmres_normr() */

