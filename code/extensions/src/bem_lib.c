#include <Python.h>                                                                                           
#define NPY_NO_DEPRECATED_API   NPY_1_7_API_VERSION                                                           
#include <numpy/arrayobject.h>
#include "ext.h"
//
static char module_docstring[] = "This is a module used to create BCT.";

static PyMethodDef module_methods[] =
{
    {NULL,NULL,0,NULL}
    
};

PyMODINIT_FUNC init_bctlib( void )
{
    //                      This is the name of the module
    //                               |
    //                              \ /
    PyObject *m = Py_InitModule3("_bemlib", module_methods, module_docstring);
    
    // Some error checking, because I'm nice like that
    if (m == NULL)
    {
        fprintf( stderr, "%s at %d : Could not initialize fucking helmholtzlib module\n", __FILE__, __LINE__ );
        return;
    }
    
    // Told you this is useful to load things in advance. Lets load numpy functionality
    import_array();
}


