
#include <list>
#include <sstream>

#include "scsi/base.h"
#include "pyscsi.h"

#define PY_ARRAY_UNIQUE_SYMBOL USCSI_PyArray_API
#include <numpy/ndarrayobject.h>

namespace {
static
PyMethodDef modmethods[] = {
    {"_GLPSParse", (PyCFunction)&PyGLPSParse, METH_VARARGS,
     "Parse a GLPS lattice file to AST form"},
    {"GLPSPrinter", (PyCFunction)&PyGLPSPrint, METH_VARARGS,
     "Print a GLPS AST to string"},
    {NULL, NULL, 0, NULL}
};

#if PY_MAJOR_VERSION >= 3
static struct PyModuleDef module = {
   PyModuleDef_HEAD_INIT,
   "uscsi._internal",
   NULL,
   -1,
   modmethods
};
#endif
}

PyMODINIT_FUNC
#if PY_MAJOR_VERSION >= 3
PyInit__internal(void)
#else
init_internal(void)
#endif
{
    try {
        if (_import_array() < 0)
            throw std::runtime_error("Failed to import numpy");

#if PY_MAJOR_VERSION >= 3
        // w/ py3 we own the module object and return NULL on import failure
        PyRef<> modref(PyModule_Create(&module));
        PyObject *mod = modref.py();
#else
        // w/ py2 we get a borrowed ref. and indicate import
        // failure by returning w/ a python exception set
        PyObject *mod = Py_InitModule("uscsi._internal", modmethods);
#endif

        if(registerModMachine(mod))
            throw std::runtime_error("Failed to initialize Machine");
        if(registerModState(mod))
            throw std::runtime_error("Failed to initialize State");

        // add States and Elements
        registerLinear();
        registerMoment();

#if PY_MAJOR_VERSION >= 3
        modref.release();
        return mod;
    }CATCH()
#else
    }CATCH2V(std::exception, RuntimeError)
#endif
}
