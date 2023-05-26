from ctypes import *
import numpy


class Triangle(Structure):
    _fields_=[("v",c_double * 9),("a",c_double),("dofs",c_int *3)]

class Quadrilateral(Structure):
    _fields_=[("v",c_double * 12),("a",c_double),("dofs",c_int *4), ("dofsp",c_int *4)]

class sP0P(Structure):
    _fields_=[("T",POINTER(Triangle))]

class sP1P(Structure):
    _fields_=[("T",POINTER(POINTER(Triangle))), ("len",c_int), ("num",POINTER(c_int))]

class sPD0(Structure):
    _fields_=[("Q",POINTER(POINTER(Quadrilateral))), ("num",POINTER(c_int)), ("len",c_int), ("dv",c_int)]

class sPD1(Structure):
    _fields_=[("Q",POINTER(POINTER(Quadrilateral))), ("num",POINTER(c_int)), ("len",c_int)]

class P0Pmesh(Structure):
    _fields_=[("sbf",POINTER(POINTER(sP0P))), ("nt",c_int),("type", POINTER(c_char_p)), ("type2",c_int) ]

class P1Pmesh(Structure):
    _fields_=[("sbf",POINTER(POINTER(sP1P))), ("ndof",c_int),("type", POINTER(c_char_p)), ("type2",c_int) ]

class PD0mesh(Structure):
    _fields_=[("sbf",POINTER(POINTER(sPD0))), ("ndof",c_int),("type", POINTER(c_char_p)), ("type2",c_int) ]

class PD1mesh(Structure):
    _fields_=[("sbf",POINTER(POINTER(sPD1))), ("ndof",c_int),("type", POINTER(c_char_p)), ("type2",c_int) ]

class Meshes(Structure):
    _fields_=[("P0P",POINTER(P0Pmesh)), ("P1P",POINTER(P1Pmesh)), ("PD0",POINTER(PD0mesh)), ("PD1",POINTER(PD1mesh)), ("ntp",c_int), ("ntb",c_int), ("ndp",c_int), ("ndb",c_int), ("nvp",c_int), ("nvb",c_int),("Vp",POINTER(POINTER(c_double))),("Vb",POINTER(POINTER(c_double))),("Dp",POINTER(POINTER(c_double))),("Db",POINTER(POINTER(c_double))),("Tp",POINTER(POINTER(c_double))),("Tb",POINTER(POINTER(c_double))) ]



class cluster(Structure):
    pass

cluster._fields_=[("start",c_int), ("size",c_int), ("sons",c_int), ("d",c_int), ("bmin",POINTER(c_double)), ("bmax",POINTER(c_double)), ("M",POINTER(Meshes)), ("son",POINTER(POINTER(cluster))) ]

class clustertree(Structure):
    _fields_=[("ndof",c_int), ("nidx",c_int), ("dof2idx",POINTER(c_int)), ("idx2dof",POINTER(c_int)), ("root",POINTER(cluster)) ]

class gsl_block_complex(Structure):
    _fields_=[("size",c_int), ("data",POINTER(c_double))]

class gsl_matrix_complex(Structure):
    _fields_=[("size21",c_int), ("size2",c_int), ("tda",c_int), ("data",POINTER(c_double)), ("block",POINTER(gsl_block_complex)), ("owner",c_int) ]

class fullmatrix(Structure):
    _fields_=[("e",c_int), ("cols",c_int), ("root",POINTER(gsl_matrix_complex)) ]

class blockcluster(Structure):
    _fields_=[("row",POINTER(cluster)), ("col",POINTER(cluster)), ("type",c_int), ("son",POINTER(POINTER(cluster))), ("block_rows",c_int), ("block_cols",c_int)  ]

class rkmatrix(Structure):
    _fields_=[("k",c_int), ("kt",c_int),("rows",c_int), ("cols",c_int), ("a",POINTER(gsl_matrix_complex)), ("b",POINTER(gsl_matrix_complex)) ]

class supermatrix(Structure):
    pass

supermatrix._fields_=[("rows",POINTER(cluster)), ("cols",POINTER(cluster)), ("block_rows",c_int), ("block_cols",c_int), ("rk",POINTER(rkmatrix)), ("full",POINTER(fullmatrix)),("s",POINTER(POINTER(supermatrix)))  ]

bct_lib=CDLL("/home/mfierrop/BCT27test/build/lib.linux-x86_64-2.7/_bctlib.so")
mesh_pointer = POINTER(Meshes)
hmat_pointer = POINTER(supermatrix)
gmres_pointer = POINTER(c_double)

def new_space(vertices, verticesb, dofs, dofsb, triangles, trianglesb, trialbf, testbf):
    
    T = (c_double * (numpy.shape(triangles)[1] * numpy.shape(triangles)[0]))()
    V = (c_double * (numpy.shape(vertices)[1] * numpy.shape(vertices)[0]))()
    D = (c_int * (numpy.shape(dofs)[1] * numpy.shape(dofs)[0]))()
    Tb = (c_double *(numpy.shape(trianglesb)[1] * numpy.shape(trianglesb)[0]))()
    Vb = (c_double * (numpy.shape(verticesb)[1] * numpy.shape(verticesb)[0]))()
    Db = (c_int * (numpy.shape(dofsb)[1] * numpy.shape(dofsb)[0]))()
    
    for i in range(numpy.shape(triangles)[0]):
        for j in range(numpy.shape(triangles)[1]):
            T[i*numpy.shape(triangles)[1]+j] = triangles[i][j]

    for i in range(numpy.shape(vertices)[0]):
        for j in range(numpy.shape(vertices)[1]):
            V[i*numpy.shape(vertices)[1]+j] = vertices[i][j]

    for i in range(numpy.shape(dofs)[0]):
        for j in range(numpy.shape(dofs)[1]):
            D[i*numpy.shape(dofs)[1]+j] = dofs[i][j]

    for i in range(numpy.shape(trianglesb)[0]):
        for j in range(numpy.shape(trianglesb)[1]):
            Tb[i*numpy.shape(trianglesb)[1]+j] = trianglesb[i][j]

    for i in range(numpy.shape(verticesb)[0]):
        for j in range(numpy.shape(verticesb)[1]):
            Vb[i*numpy.shape(verticesb)[1]+j] = verticesb[i][j]

    for i in range(numpy.shape(dofsb)[0]):
        for j in range(numpy.shape(dofsb)[1]):
            Db[i*numpy.shape(dofsb)[1]+j] = dofsb[i][j]

    space = bct_lib.new_space
    space.restype = mesh_pointer
    return space(V, Vb, D, Db, T, c_int(numpy.shape(vertices)[0]), c_int(numpy.shape(verticesb)[0]), c_int(numpy.shape(dofs)[0]), c_int(numpy.shape(dofsb)[0]), c_int(numpy.shape(triangles)[0]), c_int(numpy.shape(trianglesb)[0]), c_char_p(trialbf), c_char_p(testbf))

def hierarchical_matrix(space, eta, leafsize, eps, order, kappa):
    hmat = bct_lib.hierarchical_mat
    hmat.restype = hmat_pointer
    return hmat(space, c_double(eta), c_int(leafsize),c_double(eps), c_int(order),c_double( kappa))

def GMRES(A,b,tol,HM,MM):
    gmres = bct_lib.gmres_prec2
    gmres.restype = gmres_pointer
    xr = (c_double * len(b))()
    xi = (c_double * len(b))() 
    br = (c_double * len(b))()
    bi = (c_double * len(b))()  
    AR = (c_double * (numpy.shape(A)[1] * numpy.shape(A)[0]))() 
    AI = (c_double * (numpy.shape(A)[1] * numpy.shape(A)[0]))()
    MMR = (c_double * (numpy.shape(A)[1] * numpy.shape(A)[0]))()                                               
    MMI = (c_double * (numpy.shape(A)[1] * numpy.shape(A)[0]))()
    for i in range(len(b)):
        br[i] = numpy.real(b)[i]
        bi[i] = numpy.imag(b)[i] 
        xr[i] = 0.0
        xi[i] = 0.0 
    for i in range(len(b)):                                                                 
        for j in range(len(b)):
            AR[i*len(b)+j] = numpy.real(A)[i][j] 
            AI[i*len(b)+j] = numpy.imag(A)[i][j]
            MMR[i*len(b)+j] = numpy.real(MM)[i][j] 
            MMI[i*len(b)+j] = numpy.imag(MM)[i][j] 

    r = gmres( AR, AI, MMR, MMI, HM, br, bi, xr, xi,c_double( tol), c_int(len(b)), c_int(len(b))) 
    R = numpy.zeros(len(b))
    X = numpy.zeros(len(b), dtype ='complex')
    for i in range(len(b)):
        R[i] = r[i]
        X[i] = xr[i]+1j*xi[i]
    return R,X

def mvec(HM,x):
    TEST = bct_lib.test2
    xr = (c_double * (len(x)))()
    xi = (c_double * (len(x)))()
    for i in range(len(x)):
        xr[i] = numpy.real(x)[i]
        xi[i] = numpy.imag(x)[i]
    TEST(HM, xr, xi, c_int(len(x)))

#triangles= numpy.load('triangles.npy').tolist()
#vertices = numpy.load('vertices.npy').tolist()
#dofs = numpy.load('dofs.npy').tolist()
#trianglesb= numpy.load('trianglesb.npy').tolist()
#verticesb = numpy.load('verticesb.npy').tolist()
#dofsb = numpy.load('dofsb.npy').tolist()

#for i in range(numpy.shape(triangles)[0]):
#    for j in range(numpy.shape(triangles)[1]):
#        T[i*numpy.shape(triangles)[1]+j] = triangles[i][j]

#for i in range(numpy.shape(vertices)[0]):
#    for j in range(numpy.shape(vertices)[1]):
#        V[i*numpy.shape(vertices)[1]+j] = vertices[i][j]

#for i in range(numpy.shape(dofs)[0]):
#    for j in range(numpy.shape(dofs)[1]):
#        D[i*numpy.shape(dofs)[1]+j] = dofs[i][j]

#for i in range(numpy.shape(trianglesb)[0]):
#    for j in range(numpy.shape(trianglesb)[1]):
#        Tb[i*numpy.shape(trianglesb)[1]+j] = trianglesb[i][j]

#for i in range(numpy.shape(verticesb)[0]):
#    for j in range(numpy.shape(verticesb)[1]):
#        Vb[i*numpy.shape(verticesb)[1]+j] = verticesb[i][j]

#for i in range(numpy.shape(dofsb)[0]):
#    for j in range(numpy.shape(dofsb)[1]):
#        Db[i*numpy.shape(dofsb)[1]+j] = dofsb[i][j]

#SP = new_space(vertices, verticesb, dofs, dofsb, triangles, trianglesb, "P0d", "P0")

#H = hierarchical_matrix(SP, 0.5, 10, 0.0001, 1, 1.0)

#I = numpy.eye(26)

#b = numpy.random.rand(26)

#print "GMRES"
#res = GMRES(I,b, 1e-5, H, I)
