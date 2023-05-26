from ctypes import *
import numpy
from _bctlib import w_form
from _bctlib import assemble_matrix_dual
from _bctlib import assemble_matrix_primal
from _bctlib import assemble_matrix_mix
from _bctlib import gmres
from _bctlib import gmres_prec

def hmat( Vp, Vb, Dp, Db,  x, d, ls, eps, eta, trialbf, testbf, kappa, order ):
    
    xr = numpy.real(x).tolist()
    xi = numpy.imag(x).tolist()
    MR = numpy.zeros( (len( Vp ), len( Vp )))
    MI = numpy.zeros( (len( Vp ), len( Vp )))
    
    w_form( Vp, Vb, Dp, Db, MR, MI, xr, xi, d, ls, eps, eta, trialbf, testbf, kappa, order )
    
    return MR+(1j)*MI


def gmres_hlib(A, b, tol, maxit, P):
    xr= numpy.zeros(len(b))
    xi= numpy.zeros(len(b))
    r= numpy.zeros(maxit)
    br = numpy.real(b).tolist()
    bi = numpy.imag(b).tolist()
    AR = numpy.real(A)
    AI = numpy.imag(A)
    if P == 'None':
        gmres(AR, AI, br, bi, xr, xi, r, tol, maxit)
    else:
        PR = numpy.real(P)
        PI = numpy.imag(P)
        gmres_prec(AR, AI, PR, PI, br, bi, xr, xi, r, tol, maxit)
    return xr+1j*xi,r


def identity_primal(  V, T, D, order, k , trial, test ):
    '''
        Function to assemble matrix for helmholtz problems.
        
        B: list of arrays of vertex of shape (N elements primary mesh, 4, 3 )
        
        order: quadrature order.
        
        k: wavenumber
        '''
    if trial == "P0" and test == "P0":
        M1 = numpy.zeros( (len( D ), len( D )))
        assemble_matrix_primal( V, T, D, M1, M1, 0, trial, test, k, order)
    if trial == "P1" and test == "P1":
        M1 = numpy.zeros( (len( V ), len( V )))
        assemble_matrix_primal( V, T, D, M1, M1, 0, trial, test, k, order)
    if trial == "P1" and test == "P0":
        M1 = numpy.zeros( (len( V ), len( D )))
        assemble_matrix_primal( V, T, D, M1, M1, 0, trial, test, k, order)
    if trial == "P0" and test == "P1":
        M1 = numpy.zeros( (len( V ), len( D )))
        assemble_matrix_primal( V, T, D, M1, M1, 0, trial, test, k, order)
        M1 = M1.T
    # M is passed by reference

    return M1


def identity_dual(  Vp, Vb, Dp,Db, order, k , trial, test ):
    '''
        Function to assemble matrix for helmholtz problems.
        
        B: list of arrays of vertex of shape (N elements primary mesh, 4, 3 )
        
        order: quadrature order.
        
        k: wavenumber
        '''
    
    M1 = numpy.zeros( (len( Vp ), len( Vp )))
    
    # M is passed by reference
    assemble_matrix_dual( Vp, Vb, Dp, Db, M1, M1, 0, trial, test, k, order)
    
    return M1


def single_layer_primal( V,T, D, order, k , trial, test):
    '''
        Function to assemble matrix for helmholtz problems.
        
        B: list of arrays of vertex of shape (N elements primary mesh, 4, 3 )
        
        order: quadrature order.
        
        k: wavenumber
        '''
    
    # M is passed by reference
    if trial == "P0" and test == "P0":
        M1 = numpy.zeros( (len( D ), len( D )))
        M2 = numpy.zeros( (len( D ), len( D )))
        assemble_matrix_primal( V, T, D, M1, M2, 1, trial, test, k, order)
    if trial == "P1" and test == "P1":
        M1 = numpy.zeros( (len( V ), len( V )))
        M2 = numpy.zeros( (len( V ), len( V )))
        assemble_matrix_primal( V, T, D, M1, M2, 1, trial, test, k, order)
    M=M1+(1j)*M2

    return M

def singular_part_sld( V, T, D, order, k , trial, test):
    
    M1 = numpy.zeros( (len( D ), len( D )))
    M2 = numpy.zeros( (len( D ), len( D )))
    assemble_matrix_primal( V, T, D, M1, M2, 2, trial, test, k, order)
    M=M1+(1j)*M2
    
    return M


def single_layer_dual(Vp, Vb, Dp,Db, order, k , trial, test):
    '''
        Function to assemble matrix for helmholtz problems.
        
        B: list of arrays of vertex of shape (N elements primary mesh, 4, 3 )
        
        order: quadrature order.
        
        k: wavenumber
        '''
    if trial =="P0d" and test == "P0d":
        M1 = numpy.zeros( (len( Vp ), len( Vp )))
        M2 = numpy.zeros( (len( Vp ), len( Vp )))

    else:
        M1 = numpy.zeros( (len( Dp ), len( Dp )))
        M2 = numpy.zeros( (len( Dp ), len( Dp )))
    # M is passed by reference

    assemble_matrix_dual( Vp, Vb, Dp, Db, M1, M2, 1 , trial, test, k, order)


    M=M1+(1j)*M2

    return M


def identity_dual_primal(Vp, Vb, Dp, Db, T, order, k , trial, test ):
    '''
        Function to assemble matrix for helmholtz problems.
        
        B: list of arrays of vertex of shape (N elements primary mesh, 4, 3 )
        
        order: quadrature order.
        
        k: wavenumber   '''
    
    M1 = numpy.zeros( (len( Vp ), len( Vp )))
    
    assemble_matrix_mix( Vp, Vb, Dp, Db, T, M1, trial, test, k, order,0)
    
    return M1


def identity_dual_bary(  Vp, Vb, Dp, Db, order, k , trial, test ):
    '''
        Function to assemble matrix for helmholtz problems.
        
        B: list of arrays of vertex of shape (N elements primary mesh, 4, 3 )
        
        order: quadrature order.
        
        k: wavenumber
        '''
    
    M1 = numpy.zeros( (len( Vp ), len( Vb )))
    
    # M is passed by reference
    assemble_matrix_dual( Vp, Vb, Dp, Db, M1, M1, 2, trial, test, k, order)
    
    return M1

def projection_matrix(Vp, Vb, Dp, Db, T, order, k , trial, test ):
    '''
        Function to assemble matrix for helmholtz problems.
        
        B: list of arrays of vertex of shape (N elements primary mesh, 4, 3 )
        
        order: quadrature order.
        
        k: wavenumber
        '''
    
    # M is passed by reference
    if trial == "P1d" and test == "P1":
        M1 = numpy.zeros( (len( Dp ), len( Vb )))
        assemble_matrix_mix( Vp, Vb, Dp, Db, T, M1, "P0d", test, k, order,1)
    
    if trial == "P0d" and test == "P0":
        M1 = numpy.zeros( (len( Vp ), len( Db )))
        assemble_matrix_mix( Vp, Vb, Dp, Db, T, M1, trial, test, k, order,2)

    return M1

