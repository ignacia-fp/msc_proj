from ctypes import *
import numpy
import itertools
import time
import ctypes
import bempp.api  

class Triangle(Structure):
    _fields_=[("v",c_double * 9),("a",c_double),("dofs",c_int *3)]

class Trianglep(Structure):                                                                                    
    _fields_=[("T",POINTER(POINTER(Triangle))),("v",c_double * 9),("a",c_double),("dofs",c_int *3)] 

class Quadrilateral(Structure):
    _fields_=[("v",c_double * 12),("a",c_double),("dofs",c_int *4), ("dofsp",c_int *4)]

class sP0P(Structure):
    _fields_=[("T",POINTER(Triangle))]

class sP1P(Structure):
    _fields_=[("T",POINTER(POINTER(Triangle))), ("len",c_int), ("num",POINTER(c_int)), ("dv", c_int)]

class sPD0(Structure):
    _fields_=[("Q",POINTER(POINTER(Quadrilateral))), ("num",POINTER(c_int)), ("len",c_int), ("dv",c_int)]

class sPD1(Structure):
    _fields_=[("T",POINTER(POINTER(Trianglep))), ("num",POINTER(c_int)), ("len",c_int)]

class sRWG(Structure):                                                                                        
    _fields_=[("T1",POINTER(POINTER(Triangle))),("T2",POINTER(POINTER(Triangle))), ("num",c_int * 2),("len",c_double), ("p1",c_double * 3),("p2",c_double * 3), ("edge", c_int)] 

class sBC(Structure):
    _fields_=[("len1",c_int), ("len2",c_int), ("len3",c_int), ("len4",c_int), ("len5",c_int), ("edge",c_int), ("sbf1", POINTER(sRWG)),("sbf2", POINTER(sRWG)),("sbf3", POINTER(sRWG)),("sbf4", POINTER(sRWG)), ("sbf5", POINTER(sRWG)), ("coef1",POINTER(c_double)), ("coef2",POINTER(c_double)), ("coef3",POINTER(c_double)),("coef4",POINTER(c_double)),("coef5",POINTER(c_double)) ]
class P0Pmesh(Structure):
    _fields_=[("sbf",POINTER(POINTER(sP0P))), ("nt",c_int),("type", c_char_p), ("type2",c_int) ]

class P1Pmesh(Structure):
    _fields_=[("sbf",POINTER(POINTER(sP1P))), ("ndof",c_int),("type", c_char_p), ("type2",c_int),("globaldofs",POINTER(c_int)) ]

class PD0mesh(Structure):
    _fields_=[("sbf",POINTER(POINTER(sPD0))), ("ndof",c_int),("type", c_char_p), ("type2",c_int) ]

class PD1mesh(Structure):
    _fields_=[("sbf",POINTER(POINTER(sPD1))),("triangles",POINTER(c_int)), ("lensbf",POINTER(POINTER(c_int))),("dofsb",POINTER(c_int)),("ndofb",c_int),("ndof",c_int),("nt",c_int),("type", c_char_p), ("type2",c_int) ]

class RWGmesh(Structure):                                                                                     
    _fields_=[("sbf",POINTER(POINTER(sRWG))), ("nedge",c_int),("type", c_char_p), ("type2",c_int) ]

class BCmesh(Structure):                                                                                     
    _fields_=[("sbf",POINTER(POINTER(sBC))), ("nedge",c_int),("type", c_char_p) ]  

class Meshes(Structure):
    _fields_=[("P0P",POINTER(P0Pmesh)), ("P1P",POINTER(P1Pmesh)), ("PD0",POINTER(PD0mesh)), ("PD1",POINTER(PD1mesh)), ("RWG",POINTER(RWGmesh)), ("ntp",c_int),("ndp",c_int), ("ndb",c_int), ("nvp",c_int), ("nvb",c_int), ("ne",c_int), ("ndofb",c_int),("Vp",POINTER(POINTER(c_double))),("Vb",POINTER(POINTER(c_double))), ("trialbf", c_char_p),("testbf", c_char_p) ]

class Mesh(Structure):                                                                                      
    _fields_=[("ntp",c_int),("ntb",c_int),("ndp",c_int), ("ndb",c_int), ("nvp",c_int), ("nvb",c_int), ("nep",c_int), ("neb",c_int), ("ndofb",c_int),("Vp",POINTER(POINTER(c_double))),("Vb",POINTER(POINTER(c_double))), ("Tp",POINTER(POINTER(c_double))), ("Tb",POINTER(POINTER(c_double))), ("Dp",POINTER(POINTER(c_int))), ("Db",POINTER(POINTER(c_int))), ("Dq",POINTER(POINTER(c_int))), ("Ep",POINTER(POINTER(c_int))), ("Eb",POINTER(POINTER(c_int))),("Ep2",POINTER(POINTER(c_int))), ("Eb2",POINTER(POINTER(c_int)))]

class Space(Structure):                                                                                      
    _fields_=[("P0P",POINTER(P0Pmesh)), ("P1P",POINTER(P1Pmesh)), ("PD0",POINTER(PD0mesh)), ("PD1",POINTER(PD1mesh)), ("RWG",POINTER(RWGmesh)), ("BC",POINTER(BCmesh)),("M", POINTER(Mesh)), ("space", c_int)]

class cluster(Structure):
    pass

cluster._fields_=[("start",c_int), ("size",c_int), ("sons",c_int), ("d",c_int), ("bmin",POINTER(c_double)), ("bmax",POINTER(c_double)), ("M",POINTER(Meshes)), ("son",POINTER(POINTER(cluster))) ]

class clustertree(Structure):
    _fields_=[("ndof",c_int), ("nidx",c_int), ("dof2idx",POINTER(c_int)), ("idx2dof",POINTER(c_int)), ("root",POINTER(cluster)) ]

class gsl_block_complex(Structure):
    _fields_=[("size",c_int), ("data",POINTER(c_double))]

class gsl_matrix_complex(Structure):
    _fields_=[("size21",c_int), ("size2",c_int), ("tda",c_int), ("data",POINTER(c_double)), ("block",POINTER(gsl_block_complex)), ("owner",c_int) ]

class gsl_vector_complex(Structure):                                                                          
    _fields_=[("size",c_int), ("stride",c_int), ("data",POINTER(c_double)), ("block",POINTER(gsl_block_complex)), ("owner",c_int) ]

class fullmatrix(Structure):
    _fields_=[("e",c_int), ("cols",c_int), ("root",POINTER(gsl_matrix_complex)) ]

class blockcluster(Structure):
    pass                                                                                                      
                                                                                                              
blockcluster._fields_=[("row",POINTER(cluster)), ("col",POINTER(cluster)), ("type",c_int), ("son",POINTER(POINTER(blockcluster))), ("block_rows",c_int), ("block_cols",c_int)  ]

class rkmatrix(Structure):
    _fields_=[("k",c_int), ("kt",c_int),("rows",c_int), ("cols",c_int), ("a",POINTER(gsl_matrix_complex)), ("b",POINTER(gsl_matrix_complex)) ]

class supermatrix(Structure):
    pass

supermatrix._fields_=[("rows",POINTER(cluster)), ("cols",POINTER(cluster)), ("block_rows",c_int), ("block_cols",c_int), ("rk",POINTER(rkmatrix)), ("full",POINTER(fullmatrix)),("s",POINTER(POINTER(supermatrix))), ("in",POINTER(gsl_vector_complex)), ("out",POINTER(gsl_vector_complex)),("depth",c_int), ("nv",c_int),("index",POINTER(c_int)), ("hfi",POINTER(c_int)), ("hfo",POINTER(c_int))  ]

bem_lib=CDLL("/home/mfierrop/BEMLib2/build/lib.linux-x86_64-2.7/_beml.so")
meshes_pointer = POINTER(Meshes)              
mesh_pointer = POINTER(Mesh)                                                                  
space_pointer = POINTER(Space)
hmat_pointer = POINTER(supermatrix)                                                                           
gmres_pointer = POINTER(c_double)                                                                             
mat_pointer = POINTER(c_double)                                                                               
gsl_mat_pointer = POINTER(gsl_matrix_complex)                                                                 
gm = bem_lib.get_mat                                                                                          
gm.restype = mat_pointer                                                                                      
ge = bem_lib.get_element                                                                                      
ge.restype = c_double                                                                                         
ge_real = bem_lib.get_gsl_real                                                                                
ge_real.restype = c_double                                                                                    
ge_imag = bem_lib.get_gsl_imag                                                                                
ge_imag.restype = c_double                                                                                    
del_gsl_mat = bem_lib.del_gsl_mat                                                                             
del_gsl_mat.restype = None  

def new_mesh(grid):
    gridb = grid.barycentric_grid() 
    edges = list(grid.leaf_view.entity_iterator(1))                                                               
    edgesb = list(gridb.leaf_view.entity_iterator(1))
    vertices = list(grid.leaf_view.entity_iterator(2))                                                            
    verticesb = list(gridb.leaf_view.entity_iterator(2)) 
    elements = grid.leaf_view.elements.T                                                                          
    elementsb = gridb.leaf_view.elements.T 
    v = [x.geometry.corners.T[0].tolist() for x in vertices]                                                      
    vint = [[v.index(x.geometry.corners.T[0].tolist()),v.index(x.geometry.corners.T[1].tolist())] for x in edges] 
    edges = [[i for i,e in enumerate(elements) if set(e) >= set(x)] for x in vint]                                
    vb = [x.geometry.corners.T[0].tolist() for x in verticesb]                                                    
    vintb = [[vb.index(x.geometry.corners.T[0].tolist()),vb.index(x.geometry.corners.T[1].tolist())] for x in edgesb]
    edgesb = [[i for i,e in enumerate(elementsb) if set(e) >= set(x)] for x in vintb] 
    dofs = elements.tolist()
    dofsb = elementsb.tolist()
    vertices = grid.leaf_view.vertices.T.tolist() 
    verticesb = gridb.leaf_view.vertices.T.tolist()
    elements = list(grid.leaf_view.entity_iterator(0))
    elementsb = list(gridb.leaf_view.entity_iterator(0)) 
    triangles = [x.geometry.corners.T[0].tolist()+x.geometry.corners.T[1].tolist()+x.geometry.corners.T[2].tolist() for x in elements]
    trianglesb = [x.geometry.corners.T[0].tolist()+x.geometry.corners.T[1].tolist()+x.geometry.corners.T[2].tolist() for x in elementsb]
    dofsq = []
    idx = 0
    for i in range(0,len(dofsb)/2,3):                                                                             
        dofsq.append((dofs[idx][0], dofs[idx][1],dofs[idx][2], i))                                                
        dofsq.append((dofs[idx][1],dofs[idx][2],dofs[idx][0], i+1))                                               
        dofsq.append((dofs[idx][2],dofs[idx][0],dofs[idx][1], i+2))                                               
        idx=idx+1  
    dtype = [('a',int),('b',int),('c',int), ('d',int)]                                                            
    dofsq = numpy.array(dofsq, dtype = dtype)                                                                        
    dofsq = numpy.sort(dofsq, order = ['a','d']) 
    T = list(itertools.chain(*triangles))                                                               
    V = list(itertools.chain(*vertices))                                                                      
    D = list(itertools.chain(*dofs))                                                              
    Tb = list(itertools.chain(*trianglesb))
    Vb = list(itertools.chain(*verticesb))                                                        
    Db = list(itertools.chain(*dofsb))                                                                        
    Dq = list(itertools.chain(*dofsq))                                                                        
    Ep = list(itertools.chain(*edges)) 
    Eb = list(itertools.chain(*edgesb))
    Ep2 = list(itertools.chain(*vint))                                                                        
    Eb2 = list(itertools.chain(*vintb))
    T = (c_double * len(T))(*T)                                                                               
    V = (c_double * len(V))(*V)                                                                               
    D = (c_int * len(D))(*D)                                                                                  
    Tb = (c_double * len(Tb))(*Tb)
    Vb = (c_double * len(Vb))(*Vb)                                                                            
    Db = (c_int * len(Db))(*Db)                                                                               
    Dq = (c_int * len(Dq))(*Dq)                                                                               
    Ep = (c_int * len(Ep))(*Ep) 
    Eb = (c_int * len(Eb))(*Eb)
    Ep2 = (c_int * len(Ep2))(*Ep2)                                                                               
    Eb2 = (c_int * len(Eb2))(*Eb2)
    mesh = bem_lib.new_mesh                                                                                 
    mesh.restype = mesh_pointer 
    return mesh(V, Vb, D, Db, Dq, T, Tb, Ep, Eb, Ep2, Eb2, c_int(numpy.shape(vertices)[0]), c_int(numpy.shape(verticesb)[0]), c_int(numpy.shape(dofs)[0]), c_int(numpy.shape(dofsb)[0]), c_int(numpy.shape(dofsq)[0]), c_int(numpy.shape(triangles)[0]),c_int(numpy.shape(trianglesb)[0]), c_int(numpy.shape(edges)[0]), c_int(numpy.shape(edgesb)[0]))


def new_space(vertices, verticesb, dofs, dofsb, dofsq, triangles, edges, trialbf, testbf): 
    T = list(itertools.chain(*triangles))
    V = list(itertools.chain(*vertices))
    D = list(itertools.chain(*dofs))
    Vb = list(itertools.chain(*verticesb))
    Db = list(itertools.chain(*dofsb))
    Dq = list(itertools.chain(*dofsq))
    Ep = list(itertools.chain(*edges)) 
    T = (c_double * len(T))(*T)
    V = (c_double * len(V))(*V)
    D = (c_int * len(D))(*D)
    Vb = (c_double * len(Vb))(*Vb)
    Db = (c_int * len(Db))(*Db)
    Dq = (c_int * len(Dq))(*Dq)
    Ep = (c_int * len(Ep))(*Ep)
    space = bem_lib.new_space
    space.restype = meshes_pointer
    return space(V, Vb, D, Db, Dq, T, Ep, c_int(numpy.shape(vertices)[0]), c_int(numpy.shape(verticesb)[0]), c_int(numpy.shape(dofs)[0]), c_int(numpy.shape(dofsb)[0]), c_int(numpy.shape(dofsq)[0]), c_int(numpy.shape(triangles)[0]), c_int(numpy.shape(edges)[0]), c_char_p(trialbf), c_char_p(testbf))

def new_space2(Mesh, stype):
    space = bem_lib.new_space2
    space.restype = space_pointer
    return space(Mesh,c_char_p(stype))

def hierarchical_matrix_slp(space, eta, leafsize, eps, order, kappa):
    hmat = bem_lib.hierarchical_mat_slp
    hmat.restype = hmat_pointer
    t = numpy.zeros(1).tolist()                                                                               
    t = (c_double * len(t))(*t)
    order = (c_int * 3)(*order)
    h = hmat(space, c_double(eta), c_int(leafsize),c_double(eps), order,c_double( kappa), t)
    return h, t[0]/1000, h.contents.depth

def hierarchical_matrix_hlp(space, eta, leafsize, eps, order, kappa):                                         
    hmat = bem_lib.hierarchical_mat_hlp                                                                           
    hmat.restype = hmat_pointer                                                                               
    t = numpy.zeros(1).tolist()                                                                               
    t = (c_double * len(t))(*t) 
    order = (c_int * 3)(*order) 
    h = hmat(space, c_double(eta), c_int(leafsize),c_double(eps), order,c_double( kappa), t) 
    return h, t[0]/1000, h.contents.depth  

def del_hierarchical_matrix(H):
    del_hmat = bem_lib.del_hierarchical_mat
    del_hmat.restype = None
    del_hmat(H)

def gmres_hmat(A,b,tol,HM,MM):
    gmres = bem_lib.gmres_prec_hmat
    gmres.restype = gmres_pointer
    xr = numpy.zeros(len(b)).tolist()
    xi = numpy.zeros(len(b)).tolist()
    br = numpy.real(b).tolist()                                            
    bi = numpy.imag(b).tolist()
    AR = list(itertools.chain(*(numpy.real(A).tolist())))
    AI = list(itertools.chain(*(numpy.imag(A).tolist())))
    MMR = list(itertools.chain(*(numpy.real(MM).tolist())))                                                     
    MMI = list(itertools.chain(*(numpy.imag(MM).tolist()))) 
    xr = (c_double * len(xr))(*xr)
    xi = (c_double * len(xi))(*xi)
    br = (c_double * len(br))(*br)
    bi = (c_double * len(bi))(*bi)
    AR = (c_double * len(AR))(*AR)
    AI = (c_double * len(AI))(*AI)
    MMR = (c_double * len(MMR))(*MMR)
    MMI = (c_double * len(MMI))(*MMI)
    t = numpy.zeros(5).tolist()                                                                               
    t = (c_double * len(t))(*t)
    start_time = time.time()
    r = gmres( AR, AI, MMR, MMI, HM, br, bi, xr, xi,c_double( tol), c_int(1), c_int(len(b)), t) 
    end_time = time.time()
    R = numpy.zeros(len(b), dtype=numpy.double)
    X = numpy.zeros(len(b), dtype ='complex128')
    for i in range(len(b)):
        R[i] = r[i]
        X[i] = xr[i]+1j*xi[i]
    return R[numpy.nonzero(R)],X, t

#@profile
def gmres_full(A,b,tol):
    xr = numpy.zeros(len(b)).tolist()                                                                         
    xi = numpy.zeros(len(b)).tolist()                                                                         
    br = numpy.real(b).tolist()                                                                               
    bi = numpy.imag(b).tolist()                                                                               
    AR = list(itertools.chain(*(numpy.real(A).tolist())))                                                     
    AI = list(itertools.chain(*(numpy.imag(A).tolist())))                                                     
    xr = (c_double * len(xr))(*xr)                                                                            
    xi = (c_double * len(xi))(*xi)                                                                            
    br = (c_double * len(br))(*br)                                                                            
    bi = (c_double * len(bi))(*bi)                                                                            
    AR = (c_double * len(AR))(*AR)                                                                            
    AI = (c_double * len(AI))(*AI)                                                                            
    t = numpy.zeros(1).tolist()
    t = (c_double * len(t))(*t) 
    gmres = bem_lib.gmres
    gmres.restype = gmres_pointer
    start_time = time.time()
    r = gmres( AR, AI, br, bi, xr, xi, c_double(tol), c_int(1), c_int(len(b)), t)
    end_time = time.time()
    R = numpy.zeros(len(b), dtype = numpy.double)
    X = numpy.zeros(len(b), dtype ='complex128')
    for i in range(len(b)):
        R[i] = r[i]
        X[i] = xr[i] + 1j * xi[i]
    return X,R[numpy.nonzero(R)], t[0]/1000
    

def helmholtz_slp(space, kappa, order):                                                                       
    assemble_matrix = bem_lib.helmholtz_assemble_matrix_sl                                                    
    assemble_matrix.restype = None                                                                            
    t = numpy.zeros(1).tolist()                                                                               
    t = (c_double * len(t))(*t)                                        
    order = (c_int * 3)(*order)
    if space.contents.trialbf == "P0" and space.contents.testbf == "P0":                                      
        M1 = (c_double * (space.contents.ndp * space.contents.ndp))()                                         
        M2 = (c_double * (space.contents.ndp * space.contents.ndp))()                                         
        start_time = time.time()                                                                              
        assemble_matrix(M1, M2, space, c_double(kappa), order, c_int(0), t)                            
        end_time = time.time()                                                                                
        M = numpy.zeros(space.contents.ndp*space.contents.ndp, dtype = 'complex128')                          
        M.real = M1                                                                                           
        M.imag = M2                                                                                           
        M = M.reshape((space.contents.ndp, space.contents.ndp))                                               
    if space.contents.trialbf == "P0b" and space.contents.testbf == "P0b":                                    
        M1 = (c_double * (space.contents.ndb * space.contents.ndb))()                                         
        M2 = (c_double * (space.contents.ndb * space.contents.ndb))()                                         
        start_time = time.time()                                                                              
        assemble_matrix(M1, M2, space, c_double(kappa), order,  t)                                            
        end_time = time.time()                                                                                
        M = numpy.zeros(space.contents.ndb*space.contents.ndb, dtype = 'complex128')                          
        M.real = M1                                                                                           
        M.imag = M2                                                                                           
        M = M.reshape((space.contents.ndb, space.contents.ndb))                                               
    if (space.contents.trialbf == "P1" and space.contents.testbf == "P1"):                                    
        M1 = (c_double * (space.contents.nvp * space.contents.nvp))()                                         
        M2 = (c_double * (space.contents.nvp * space.contents.nvp))()                                         
        start_time = time.time()                                                                              
        assemble_matrix(M1, M2, space, c_double(kappa), order, c_int(0), t)                            
        end_time = time.time()                                                                                
        M = numpy.zeros(space.contents.nvp*space.contents.nvp, dtype = 'complex128')                          
        M.real = M1                                                                                           
        M.imag = M2                                                                                           
        M = M.reshape((space.contents.nvp, space.contents.nvp))                                               
    if (space.contents.trialbf == "P1b" and space.contents.testbf == "P1b"):                                  
        M1 = (c_double * (space.contents.nvb * space.contents.nvb))()                                         
        M2 = (c_double * (space.contents.nvb * space.contents.nvb))()                                         
        start_time = time.time()                                                                              
        assemble_matrix(M1, M2, space, c_double(kappa), order, c_int(0), t)                            
        end_time = time.time()                                                                                
        M = numpy.zeros(space.contents.nvb*space.contents.nvb, dtype = 'complex128')                          
        M.real = M1                                                                                           
        M.imag = M2                                                                                           
        M = M.reshape((space.contents.nvb, space.contents.nvb))                                               
    return M, t[0]/1000          

def helmholtz_hlp(space, kappa, order, opt):                                                                  
    assemble_matrix = bem_lib.helmholtz_assemble_matrix_hl                                                    
    assemble_matrix.restype = None                                                                            
    t = numpy.zeros(1).tolist()                                                                               
    t = (c_double * len(t))(*t)                                                                               
    order = (c_int * 3)(*order)                                                                               
    if space.contents.trialbf == "P0" and space.contents.testbf == "P0":                                      
        M1 = (c_double * (space.contents.ndp * space.contents.ndp))()                                         
        M2 = (c_double * (space.contents.ndp * space.contents.ndp))()                                         
        assemble_matrix(M1, M2, space, c_double(kappa), order,  t, c_int(opt))                                
        M = numpy.zeros(space.contents.ndp*space.contents.ndp, dtype = 'complex128')                          
        M.real = M1                                                                                           
        M.imag = M2                                                                                           
        M = M.reshape((space.contents.ndp, space.contents.ndp))                                               
    if space.contents.trialbf == "P0b" and space.contents.testbf == "P0b":                                    
        M1 = (c_double * (space.contents.ndb * space.contents.ndb))()                                         
        M2 = (c_double * (space.contents.ndb * space.contents.ndb))()                                         
        assemble_matrix(M1, M2, space, c_double(kappa), order,  t,c_int(opt))                                 
        M = numpy.zeros(space.contents.ndb*space.contents.ndb, dtype = 'complex128')                          
        M.real = M1                                                                                           
        M.imag = M2                                                                                           
        M = M.reshape((space.contents.ndb, space.contents.ndb))                                               
    if (space.contents.trialbf == "P1" and space.contents.testbf == "P1"):                                    
        M1 = (c_double * (space.contents.nvp * space.contents.nvp))()                                         
        M2 = (c_double * (space.contents.nvp * space.contents.nvp))()                                         
        assemble_matrix(M1, M2, space, c_double(kappa), order, t, c_int(opt))                                 
        M = numpy.zeros(space.contents.nvp*space.contents.nvp, dtype = 'complex128')                          
        M.real = M1                                                                                           
        M.imag = M2                                                                                           
        M = M.reshape((space.contents.nvp, space.contents.nvp))                                               
    if (space.contents.trialbf == "P1b" and space.contents.testbf == "P1b"):                                  
        M1 = (c_double * (space.contents.nvb * space.contents.nvb))()                                         
        M2 = (c_double * (space.contents.nvb * space.contents.nvb))()                                         
        assemble_matrix(M1, M2, space, c_double(kappa), order, t,c_int(opt))                                  
        M = numpy.zeros(space.contents.nvb*space.contents.nvb, dtype = 'complex128')                          
        M.real = M1                                                                                           
        M.imag = M2                                                                                           
        M = M.reshape((space.contents.nvb, space.contents.nvb))                                               
    if (space.contents.trialbf == "P1d" and space.contents.testbf == "P1d"):                                  
        M1 = (c_double * (space.contents.ntp * space.contents.ntp))()                                         
        M2 = (c_double * (space.contents.ntp * space.contents.ntp))()                                         
        assemble_matrix(M1, M2, space, c_double(kappa), order, t, c_int(opt))                                 
        M = numpy.zeros(space.contents.ntp*space.contents.ntp, dtype = 'complex128')                          
        M.real = M1                                                                                           
        M.imag = M2                                                                                           
        M = M.reshape((space.contents.ntp, space.contents.ntp))                                               
    return M, t[0]/1000



#@profile
def helmholtz_sld(space, kappa, order):                                                                       
    assemble_matrix = bem_lib.helmholtz_assemble_matrix_sl                                             
    assemble_matrix.restype = None
    order = (c_int * 3)(*order)
    t = numpy.zeros(1).tolist()                                                                               
    t = (c_double * len(t))(*t)
    M1 = (c_double * (space.contents.nvp * space.contents.nvp))()                                         
    M2 = (c_double * (space.contents.nvp * space.contents.nvp))() 
    start_time = time.time()
    assemble_matrix(M1, M2, space, c_double(kappa), order, c_int(1), t) 
    end_time = time.time()
    M = numpy.zeros(space.contents.nvp*space.contents.nvp, dtype = 'complex128')                              
    M.real = M1                                                                                               
    M.imag = M2                                                                                               
    M = M.reshape((space.contents.nvp, space.contents.nvp))
    return M , t[0]/1000

def maxwell_T(trial, test, kappa,order):                                                                        
    assemble_matrix = bem_lib.maxwell_T                                                                
    assemble_matrix.restype = None                                                                            
    t = numpy.zeros(1).tolist()                                                                               
    t = (c_double * len(t))(*t)                                                                               
    M1 = (c_double * (trial.contents.M.contents.nep * test.contents.M.contents.nep))()                                               
    M2 = (c_double * (trial.contents.M.contents.nep * test.contents.M.contents.nep))()                                               
    M3 = (c_double * (trial.contents.M.contents.nep * test.contents.M.contents.nep))()                                               
    M4 = (c_double * (trial.contents.M.contents.nep * test.contents.M.contents.nep))() 
    order = (c_int * 4)(*order)                                                                               
    start_time = time.time()                                                                                  
    assemble_matrix(M1, M2, M3, M4, trial, test, c_double(kappa), order, t, c_int(1))                                       
    end_time = time.time()                                                                                    
    M5 = numpy.zeros(trial.contents.M.contents.nep * test.contents.M.contents.nep, dtype = 'complex128')                                
    M5.real = M1                                                                                               
    M5.imag = M2                    
    M6 = numpy.zeros(trial.contents.M.contents.nep * test.contents.M.contents.nep, dtype = 'complex128')
    M6.real = M3                                                                                              
    M6.imag = M4
    M5 = M5.reshape((trial.contents.M.contents.nep, test.contents.M.contents.nep))
    M6 = M6.reshape((trial.contents.M.contents.nep, test.contents.M.contents.nep))
    M = (-1/(1j*kappa))*M5-(1j*kappa)*M6
    return M , t[0]/1000 

def maxwell_weakly(space,kappa,order):
    assemble_matrix = bem_lib.maxwell_weakly
    assemble_matrix.restype = None
    t = numpy.zeros(1).tolist()                                                                               
    t = (c_double * len(t))(*t)                                                                               
    M1 = (c_double * (space.contents.ne * space.contents.ne))()                                             
    M2 = (c_double * (space.contents.ne * space.contents.ne))()
    order = (c_int * 4)(*order) 
    start_time = time.time()
    assemble_matrix(M1, M2, space, c_double(kappa), order, t, c_int(1))                                
    end_time = time.time()                                                                                    
    M = numpy.zeros(space.contents.nep*space.contents.nep, dtype = 'complex128')                              
    M.real = M1                                                                                               
    M.imag = M2                                                                                               
    M = M.reshape((space.contents.nep, space.contents.nep))  
    M = -(1j*kappa)*M
    return M , t[0]/1000

def maxwell_hypersingular(trial, test, kappa,order):                                                                        
    assemble_matrix = bem_lib.maxwell_hypersingular                                                                         
    assemble_matrix.restype = None                                                                            
    t = numpy.zeros(1).tolist()                                                                               
    t = (c_double * len(t))(*t)                                                                               
    M1 = (c_double * (trial.contents.M.contents.nep * test.contents.M.contents.nep))()                                               
    M2 = (c_double * (trial.contents.M.contents.nep * test.contents.M.contents.nep))()                                               
    order = (c_int * 4)(*order)                                                                               
    start_time = time.time()                                                                                  
    assemble_matrix(M1, M2, trial, test, c_double(kappa), order, t, c_int(1))                                       
    end_time = time.time()                                                                                    
    M = numpy.zeros(trial.contents.M.contents.nep*test.contents.M.contents.nep, dtype = 'complex128')                                
    M.real = M1                                                                                               
    M.imag = M2                                                                                               
    M = M.reshape((trial.contents.M.contents.nep,test.contents.M.contents.nep))  
    M = (-1/(1j*kappa))*M
    return M , t[0]/1000 

def maxwell_Id(space,kappa,order):                                                                             
    assemble_matrix = bem_lib.maxwell_id                                                                        
    assemble_matrix.restype = None                                                                            
    t = numpy.zeros(1).tolist()                                                                               
    t = (c_double * len(t))(*t)
    if space.contents.trialbf == "RWG" and space.contents.testbf == "RWG": 
        M1 = (c_double * (space.contents.ne * space.contents.ne))()                                                                                             
        order = (c_int * 3)(*order)                                                                               
        start_time = time.time()                                                                                  
        assemble_matrix(M1, space, c_double(kappa), order, t, c_int(1))                                       
        end_time = time.time()                                                                                    
        M = numpy.zeros(space.contents.ne*space.contents.ne, dtype = 'complex128')                                
        M.real = M1                                                                                                                                                                                     
        M = M.reshape((space.contents.ne, space.contents.ne))
    if space.contents.trialbf == "RWG" and space.contents.testbf == "P0":
        M1 = (c_double * (space.contents.ne * space.contents.ntp))()                                           
        order = (c_int * 3)(*order)                                                                           
        start_time = time.time()                                                                              
        assemble_matrix(M1, space, c_double(kappa), order, t, c_int(1))                                       
        end_time = time.time()                                                                                
        M = numpy.zeros(space.contents.ne*space.contents.ntp, dtype = 'complex128')                            
        M.real = M1                                                                                           
        M = M.reshape((space.contents.ne, space.contents.ntp))  
    return M , t[0]/1000 


def maxwell_lb(space,kappa,order):
    assemble_matrix = bem_lib.maxwell_lb
    assemble_matrix.restype = None
    t = numpy.zeros(1).tolist()
    t = (c_double * len(t))(*t)
    M1 = (c_double * (space.contents.ne * space.contents.ne))()
    order = (c_int * 3)(*order)
    start_time = time.time()
    assemble_matrix(M1, space, c_double(kappa), order, t, c_int(1))
    end_time = time.time()
    M = numpy.zeros(space.contents.ne*space.contents.ne, dtype = 'complex128')
    M.real = M1 
    M = M.reshape((space.contents.ne, space.contents.ne))
    return M, t[0]/1000   

def project_matrix(space,M):
    project = bem_lib.helmholtz_project_mat
    project.restype = None                                                                            
    t = numpy.zeros(1).tolist()                                                                               
    t = (c_double * len(t))(*t) 
    MRb = list(itertools.chain(*(numpy.real(M).tolist())))                                                     
    MIb = list(itertools.chain(*(numpy.imag(M).tolist())))
    MRb = (c_double * len(MRb))(*MRb)                                                                            
    MIb = (c_double * len(MIb))(*MIb)
    MR = (c_double * (space.contents.ntp * space.contents.ntp))()                                         
    MI = (c_double * (space.contents.ntp * space.contents.ntp))()
    project(MRb, MIb, MR, MI, space, t)
    A = numpy.zeros(space.contents.ntp*space.contents.ntp, dtype = 'complex128')                              
    A.real = MR                                                                                               
    A.imag = MI                                                                                               
    A = A.reshape((space.contents.ntp, space.contents.ntp)) 
    return A , t[0]/1000 

def average_matrix(space):                                                                       
    assemble_matrix = bem_lib.helmholtz_assemble_matrix_id
    assemble_matrix.restype = None             
    t = numpy.zeros(1).tolist()                                                                               
    t = (c_double * len(t))(*t)
    M1 = (c_double * (space.contents.ntp * space.contents.nvb))()                                                                                   
    assemble_matrix(M1, space, c_int(0), c_int(5), t)                                                                             
    M = numpy.zeros(space.contents.ntp*space.contents.ntp, dtype='complex128')                                
    M.real = M1                                                                                               
    M = M.reshape((space.contents.ntp, space.contents.ntp)) 
    return M , t[0]/1000

def coupling_matrix(space):                                                                                    
    assemble_matrix = bem_lib.helmholtz_assemble_matrix_id                                                    
    assemble_matrix.restype = None                                                                            
    t = numpy.zeros(1).tolist()                                                                               
    t = (c_double * len(t))(*t)                                                                               
    M1 = (c_double * (space.contents.ntp * space.contents.ndb))()                                             
    assemble_matrix(M1, space, c_int(0), c_int(6), t)       
    M = numpy.zeros(space.contents.ntp*space.contents.ndb, dtype='complex128')                                
    M.real = M1                                                                                               
    M = M.reshape((space.contents.ntp, space.contents.ndb))                                                    
    return M , t[0]/1000 

def identity_dual_primal(space):
    assemble_matrix = bem_lib.helmholtz_assemble_matrix_id
    assemble_matrix.restype = None
    t = numpy.zeros(1).tolist()                                                                               
    t = (c_double * len(t))(*t)
    if(space.contents.trialbf == "P0d" and space.contents.testbf == "P1"):
        M1 = (c_double * (space.contents.nvp * space.contents.nvp))()                                                                                 
        assemble_matrix(M1, space, c_int(1), c_int(3),t)    
        M = numpy.zeros(space.contents.nvp*space.contents.nvp, dtype='complex128')                             
        M.real = M1                                                                                           
        M = M.reshape((space.contents.nvp, space.contents.nvp))   
    else:
        M1 = (c_double * (space.contents.ntp * space.contents.ntp))()                                         
        assemble_matrix(M1, space, c_int(2), c_int(4),t)          
        M = numpy.zeros(space.contents.ntp*space.contents.ntp, dtype='complex128')                            
        M.real = M1                                                                                                
        M = M.reshape((space.contents.ntp, space.contents.ntp))                                             
    return M , t[0]/1000

def identity_primal(space):
    assemble_matrix = bem_lib.helmholtz_assemble_matrix_id 
    assemble_matrix.restype = None                                                                            
    t = numpy.zeros(1).tolist()                                                                               
    t = (c_double * len(t))(*t)
    if(space.contents.trialbf == "P1" and space.contents.testbf == "P0"):                                    
        M1 = (c_double * (space.contents.nvp * space.contents.ntp))()                                         
        assemble_matrix(M1, space, c_int(1), c_int(0),t)                                                      
        M = numpy.zeros(space.contents.nvp*space.contents.ntp, dtype='complex128')                            
        M.real = M1                                                                                           
        M = M.reshape((space.contents.nvp, space.contents.ntp))  
    if(space.contents.trialbf == "P1b" and (space.contents.testbf == "P0b" or space.contents.testbf == "P0")):                                     
        M1 = (c_double * (space.contents.nvb * space.contents.ndb))()                                         
        assemble_matrix(M1, space, c_int(2), c_int(0),t)         
        M = numpy.zeros(space.contents.nvb*space.contents.ndb, dtype='complex128')                            
        M.real = M1                                                                                           
        M = M.reshape((space.contents.nvb, space.contents.ndb))                                               
    return M , t[0]/1000

def mat_mat(HM,A):
    mtmt = bem_lib.mtmt                                                                                       
    mtmt.restype = gsl_mat_pointer                                                                            
    t = numpy.zeros(3).tolist()                                                                               
    t = (c_double * len(t))(*t)                                                                               
    sz = numpy.shape(A)[0]                                                                                    
    AR = numpy.array(list(itertools.chain(*(numpy.real(A).tolist()))))                                        
    AI = numpy.array(list(itertools.chain(*(numpy.imag(A).tolist()))))                                        
    M = mtmt(HM, AR.ctypes.data_as(ctypes.POINTER(ctypes.c_double)), AI.ctypes.data_as(ctypes.POINTER(ctypes.c_double)), t, c_int(sz))
    for i in range(sz):                                                                                       
        for j in range(sz):                                                                                   
            A.real[i][j] = ge_real(M,c_int(i),c_int(j))                                                       
            A.imag[i][j] = ge_imag(M,c_int(i),c_int(j))                                                       
    del AR                                                                                                    
    del AI                                                                                                    
    del_gsl_mat(M)                                                                                            
    t = [t[0]/1000,t[1]/1000, t[2]/1000]  
    return   t

def mat_vec(HM,b):                                                                                            
    mv = bem_lib.mv                                                                                       
    mv.restype = None                                                                                       
    t = numpy.zeros(3).tolist()                                                                               
    t = (c_double * len(t))(*t)                                                                               
    sz = len(b)                                                                                       
    br = numpy.real(b).tolist()                                                                               
    bi = numpy.imag(b).tolist()
    br = (c_double * len(br))(*br)                                                                            
    bi = (c_double * len(bi))(*bi)                                                                            
    mv(HM, br, bi, t, c_int(sz))                                                                            
    b2 = numpy.zeros( (sz), dtype='complex128')                                                                                                                        
    b2.real = br                                                                                              
    b2.imag = bi 
    t = [t[0]/1000,t[1]/1000, t[2]/1000]
    return b2, t


