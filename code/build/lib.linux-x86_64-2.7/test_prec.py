from bctlib.bctlib.core import identity_dual_bary
from bctlib.bctlib.core import identity_dual
from bctlib.bctlib.core import identity_dual_primal
from bctlib.bctlib.core import identity_primal
from bctlib.bctlib.core import single_layer_dual
from bctlib.bctlib.core import single_layer_primal
from bctlib.bctlib.core import projection_matrix
from bctlib.bctlib.core import singular_part_sld
from bctlib.bctlib.core import hmat
from bctlib.bctlib.core import gmres_hlib 
from test3 import new_space
from test3 import hierarchical_matrix
from test3 import GMRES
from test3 import mvec
from sys import getsizeof
from scipy.sparse import csc_matrix  
from scipy.sparse.linalg import inv
import matplotlib.pyplot as plt
import scipy.sparse.linalg 
import numpy as np
import bempp.api as bem                                                                                       
import numpy as np                                                                                            
from numpy.linalg import norm                                                                                 
from bempp.api.operators.boundary.sparse import identity                                                      
from bempp.api.assembly.discrete_boundary_operator import  SparseDiscreteBoundaryOperator                                                                            
from scipy.sparse import csc_matrix                                                                           
from bempp.api.assembly.discrete_boundary_operator import InverseSparseDiscreteBoundaryOperator                                                                     
from bempp.api.assembly import as_matrix                                                                      
from bempp.api.assembly import BoundaryOperator                                                               
import time
from scipy.sparse.linalg import LinearOperator

class it_counter(object):                                                                                     
                                                                                                              
    def __init__(self, store_residuals):                                                                      
        self._count = 0                                                                                       
        self._store_residuals = store_residuals                                                               
        self._residuals = []                                                                                  
                                                                                                              
    def __call__(self, x):                                                                                    
        self._count += 1                                                                                      
        if self._store_residuals:                                                                             
            self._residuals.append(np.linalg.norm(x))                                                         
                                                                                                              
                                                                                                              
    @property                                                                                                 
    def count(self):                                                                                          
        return self._count                                                                                    
                                                                                                              
    @property                                                                                                 
    def residuals(self):                                                                                      
        return self._residuals 

#def gmres(A, b, tol=1E-5, restart=None, maxiter=None, use_strong_form=False, return_residuals=False):  

#    callback = it_counter(return_residuals)                                                                   
#    print("Starting GMRES iteration")                                                                         
                                                                                                              
#    start_time = time.time()                                                                                  
#    x, info = scipy.sparse.linalg.gmres(A, b,                                                          
#                                    tol=tol, restart=restart, maxiter=maxiter,callback=callback)              
                                                                                                              
#    end_time = time.time()                                                                                    
    # bempp.api.LOGGER.info("GMRES finished in {0} iterations and took {1:.2E} sec.".format(callback.count, end_time - start_time))
#    print("GMRES finished in", callback.count,"iterations and took",  end_time - start_time,"sec.")

#    return x, info, callback.residuals 

def neumann_data(x, n, domain_index, result):                                                                 
    result[0] = 1j * kappa * n[0] * np.exp(1j *kappa * x[0])   

for num in range(2,3):
    eta = 0.5
    eps = 0.002
    kappa = 1.0 
    n = 20
    precision = 8.*num                                                                                            
    h = 2.0 * np.pi / (precision * kappa)                                                                     
    grid = bem.shapes.sphere(h=h)
    gridb = grid.barycentric_grid()                                                                           
    gridbb = gridb.barycentric_grid()                                                                         
    P0 = bem.function_space(grid, "DP", 0)                                                                    
    P1 = bem.function_space(grid, "P", 1)                                                                     
    P0d = bem.function_space(grid, "DUAL", 0)                                                                 
    P1b = bem.function_space(gridb, "P", 1)                                                                   
    P0b = bem.function_space(gridb, "DP", 0)                                                                  
    P0bd = bem.function_space(gridb, "DUAL", 0)                                                               
    P0bb = bem.function_space(gridbb, "DP", 0) 
    T2 =[]
    T = []   
    elements0 = list(grid.leaf_view.entity_iterator(0))                                                       
    elements1 = list(grid.leaf_view.entity_iterator(1))                                                       
    elements0b = list(gridb.leaf_view.entity_iterator(0))                                                     
    elements1b = list(gridb.leaf_view.entity_iterator(1))                                                     
    elements0bb = list(gridbb.leaf_view.entity_iterator(0))                                                   
    elements1bb = list(gridbb.leaf_view.entity_iterator(1))                                                   
    n0 = len(elements0)                                                                                       
    n1 = len(elements1)                                                                                       
    n0b = len(elements0b)                                                                                     
    n1b = len(elements1b)                                                                                     
    n0bb = len(elements0bb)                                                                                   
    n1bb = len(elements1bb)                                                                                   
    triangles = []                                                                                            
    trianglesb = []                                                                                           
    trianglesbb = []                                                                                          
    c1 = []                                                                                                   
    c1b = []                                                                                                  
    c1bb = []                                                                                                 
    v0 = []                                                                                                   
                                                                                                              
    edges = np.zeros((n1,2))                                                                                  
    edgesb = np.zeros((n1b,2))                                                                                
    edgesbb = np.zeros((n1bb,2))                                                                              
                                                                                                              
    for i in range(0,n0):                                                                                     
        triangles.append(elements0[i].geometry.corners.T[0].tolist()+elements0[i].geometry.corners.T[1].tolist()+elements0[i].geometry.corners.T[2].tolist())
        v0.append(elements0[i].geometry.volume)                                                               
    for i in range(0,n1):                                                                                     
        c1.append(elements1[i].geometry.corners.T.tolist())                                                   
                                                                                                              
    for i in range(0,n0b):                                                                                    
        trianglesb.append(elements0b[i].geometry.corners.T[0].tolist()+elements0b[i].geometry.corners.T[1].tolist()+elements0b[i].geometry.corners.T[2].tolist())
    for i in range(0,n1b):                                                                                    
        c1b.append(elements1b[i].geometry.corners.T.tolist())                                                 
                                                                                                              
    for i in range(0,n0bb):                                                                                   
        trianglesbb.append(elements0bb[i].geometry.corners.T[0].tolist()+elements0bb[i].geometry.corners.T[1].tolist()+elements0bb[i].geometry.corners.T[2].tolist())
    for i in range(0,n1bb):                                                                                   
        c1bb.append(elements1bb[i].geometry.corners.T.tolist())                                               
                                                                                                              
    dofs = grid.leaf_view.elements.T.tolist()                                                                 
    vertices = grid.leaf_view.vertices.T.tolist()                                                             
                                                                                                              
    dofsb = gridb.leaf_view.elements.T.tolist()                                                               
    verticesb = gridb.leaf_view.vertices.T.tolist()                                                           
                                                                                                              
    dofsbb = gridbb.leaf_view.elements.T.tolist()                                                             
    verticesbb = gridbb.leaf_view.vertices.T.tolist()                                                         
                                                                                                              
    for i in range(0,len(vertices)):                                                                          
        for j in range(0,n1):                                                                                 
            if (c1[j][0] == vertices[i]):                                                                     
                edges[j][0] =i;                                                                               
                if (c1[j][1] == vertices[i]):                                                                 
                    edges[j][1] =i;                                                                           
                                                                                                              
    for i in range(0,len(verticesb)):                                                                         
        for j in range(0,n1):                                                                                 
            if (c1b[j][0] == verticesb[i]):                                                                   
                edgesb[j][0] =i;                                                                              
                if (c1b[j][1] == verticesb[i]):                                                               
                    edgesb[j][1] =i;                                                                          
                                                                                                              
    for i in range(0,len(verticesbb)):                                                                        
        for j in range(0,n1b):                                                                                
            if (c1bb[j][0] == verticesbb[i]):                                                                 
                edgesbb[j][0] =i;                                                                             
                if (c1bb[j][1] == verticesbb[i]):                                                             
                    edgesbb[j][1] =i;                                                                         
    edges.tolist()                                                                                            
    edgesb.tolist()                                                                                           
    edgesbb.tolist()
    SP = new_space(vertices, verticesb, dofs, dofsb, triangles, trianglesb,"P0d","P0")
    H = hierarchical_matrix( SP, eta, n, eps, 1, kappa)
#    I_dp = bem.assembly.InverseSparseDiscreteBoundaryOperator(SparseDiscreteBoundaryOperator(csc_matrix(identity_dual_primal(vertices,verticesb, dofs,dofsb,triangles, 1, kappa , "P0d", "P1" ))))
    
    I_none = InverseSparseDiscreteBoundaryOperator(bem.operators.boundary.sparse.identity(P1, P1, P1).weak_form())
    I_none = as_matrix(I_none) 
    S_none = as_matrix(bem.operators.boundary.helmholtz.single_layer(P1, P1, P1,kappa).weak_form() )                                                                                           
    S_0d = as_matrix(bem.operators.boundary.helmholtz.single_layer(P0d, P0d, P0d, kappa).weak_form())
    HP_none =as_matrix( bem.operators.boundary.helmholtz.hypersingular(P1, P0, P1, kappa).weak_form())

    trace_fun = bem.GridFunction(P0, fun = neumann_data)                                                      
    ID_none = bem.operators.boundary.sparse.identity(P0, P0, P1).weak_form()                                  
    ADL_none = bem.operators.boundary.helmholtz.adjoint_double_layer(P0, P0, P1, kappa).weak_form()           
    rhs_none = (0.5*ID_none-ADL_none)*(trace_fun).projections()                                               
    P0op = bem.operators.boundary.sparse.identity(P0b,P0d,P0d).weak_form()*InverseSparseDiscreteBoundaryOperator(bem.operators.boundary.sparse.identity(P0b, P0b, P0b).weak_form())
    
#    start_time = time.time()                                                                                  
    S_0dp = SparseDiscreteBoundaryOperator(csc_matrix(single_layer_primal(verticesb, trianglesb, dofsb, 1, kappa , "P0", "P0")))
    S_0dp2 = as_matrix(S_0dp)
    S_0dp = as_matrix(P0op*S_0dp*P0op.T)
#    end_time = time.time()                                                                                    
#    T2.append(end_time-start_time) 

#    start_time = time.time()                                                                                  
  #  H=hmat(vertices, verticesb, dofs, dofsb, rhs_none, 3, n, eps, eta, "P0d", "P0d", kappa,1) +np.diag(np.diag(as_matrix(S_0dp)))                        
#    end_time = time.time()                                                                                    
#    plt.imshow(H.real)                                                                                        
#    plt.show()
#    print np.count_nonzero(H)
#    H = SparseDiscreteBoundaryOperator(csc_matrix(H))                                                         
#    T2.append(end_time-start_time) 
    
#    start_time = time.time()
#    S_0d2 = single_layer_dual(vertices, verticesb, dofs, dofsb, 1, kappa , "P0d", "P0d")
#    S_0d2 = np.diag(np.diag(as_matrix(S_0d)))+S_0d2
#    end_time = time.time()
#    T2.append(end_time-start_time)

#    start_time = time.time()
#    S_none2 = single_layer_primal( vertices,triangles, dofs, 2, kappa , "P1", "P1")
#    end_time = time.time()
#    T2.append(end_time-start_time)
#    x1,r1=gmres_hlib(HP_none,rhs_none,1e-5,50,np.dot(I_none,np.dot(S_0d,I_none)))  
    HM=hmat(vertices, verticesb, dofs, dofsb, rhs_none, 3, n, eps, eta, "P0d", "P0d", kappa,1)  
   # S_0d = S_0d-np.diag(np.diag(S_0d))
    x,r=gmres_hlib(HP_none,rhs_none,1e-5,50,np.dot(I_none,np.dot(S_0dp,I_none)))
    x2,r2=gmres_hlib(HP_none,rhs_none,1e-5,50,'None')
    x3,r3=gmres_hlib(np.dot(np.dot(I_none,np.dot(S_0d,I_none)),HP_none),np.dot(np.dot(I_none,np.dot(S_0d,I_none)),rhs_none),1e-5,50,'None')
#    x3,r3=gmres_hlib(np.dot(S_0d2,HP_none),np.dot(S_0d2,rhs_none),1e-4,50,'None')
    [r4,x4] = GMRES(HP_none,rhs_none, 1e-5, H, I_none)
    mvec(H,rhs_none)
#   P_none = I_none*S_none*I_none*HP_none
 #   P_none2 = I_none*S_none2*I_none*HP_none
#    P_0d = I_dp*S_0d*I_dp*HP_none
#    P_0dp = I_dp*S_0dp*I_dp*HP_none
#    P_0d2 = I_dp*S_0d2*I_dp*HP_none
#    P_H = I_dp*H*I_dp*HP_none

    
#    sz = np.array((getsizeof(S_none), getsizeof(S_none2), getsizeof(S_0d), getsizeof(S_0dp), getsizeof(S_0d2), getsizeof(H)))

#    rhs_none_p = I_none*S_none2*I_none*rhs_none
#    rhs_0d = I_none*S_none*I_none*rhs_none
#    rhs_0d2 = I_dp*S_0d2*I_dp*rhs_none
#    rhs_H = I_dp*H*I_dp*rhs_none

#    c_none = np.linalg.cond(as_matrix(P_none))
#    c_none2 = np.linalg.cond(as_matrix(P_none2))
#    c_0d = np.linalg.cond(as_matrix(P_0d))
#    c_0dp = np.linalg.cond(as_matrix(P_0dp))
#    c_0d2 = np.linalg.cond(as_matrix(P_0d2))
#    c_nn = np.linalg.cond(as_matrix(HP_none))
#    c_H = np.linalg.cond(as_matrix(P_H))

#    C = np.array((c_none, c_none2, c_0d, c_0dp, c_0d2, c_H))
    tolerance = 1e-5                                                                                           
#    restart = 1000                                                                                                
#    maxiter = 10000 

#    x_none, info_none, res_none = gmres(HP_none, rhs_none, tol=tolerance, restart=restart, return_residuals=True, maxiter=maxiter)
#    x_none_p, info_none_p, res_none_p = gmres(P_none, rhs_none_p, tol=tolerance, restart=restart, return_residuals=True, maxiter=maxiter)
#    x_0d, info_0d, res_0d = gmres(P_0dp, rhs_0d, tol=tolerance, restart=restart, return_residuals=True, maxiter=maxiter)
#    x_0d2, info_0d2, res_0d2 = gmres(P_0d2, rhs_0d2, tol=tolerance, restart=restart, return_residuals=True, maxiter=maxiter)
#    x_H, info_H, res_H = gmres(P_H, rhs_H, tol=tolerance, restart=restart, return_residuals=True, maxiter=maxiter)

#    print "c_none, c_none2, c_0d, c_0dp, c_0d2, c_H"
#    print c_none, c_none2, c_0d, c_0dp, c_0d2,c_H
#    print sz
#    print "S_0dp, S_0d2, S_0b2, S_none2, S_dd"
#    print T2
#    np.save('Sphere'+str(num)+'/x_none.npy',x_none); 
#    np.save('Sphere'+str(num)+'/x_none_p.npy',x_none_p);
#    np.save('Sphere'+str(num)+'/x_0d.npy',x_0d);
#    np.save('Sphere'+str(num)+'/x_0d2.npy',x_0d2);
#    np.save('Sphere'+str(num)+'/x_H.npy',x_H);
#    np.save('Sphere'+str(num)+'/res_none.npy',res_none);                                                          
#    np.save('Sphere'+str(num)+'/res_none_p.npy',res_none_p);
#    np.save('Sphere'+str(num)+'/res_0d.npy',res_0d);                                                              
#    np.save('Sphere'+str(num)+'/res_0d2.npy',res_0d2);
#    np.save('Sphere'+str(num)+'/res_H.npy',res_H);
#    np.save('Sphere'+str(num)+'/condition_number.npy',C);
#    np.save('Sphere'+str(num)+'/time2.npy',T2);
