
#include <list>
#include <sstream>

#include "flame/base.h"
#include "pyflame.h"

#define PY_ARRAY_UNIQUE_SYMBOL FLAME_PyArray_API

#define NPY_NO_DEPRECATED_API NPY_1_6_API_VERSION
#include <numpy/ndarrayobject.h>

namespace {

#define FLAME_LOGGER_NAME "flame.machine"

// redirect flame logging to python logger
// Assumes interpreter lock is held
struct PyLogger : public Machine::Logger
{
    PyRef<> logger;
    virtual void log(const Machine::LogRecord &r)
    {
        if(logger.py()) {
            std::string msg(r.strm.str());
            size_t pos = msg.find_last_not_of('\n');
            if(pos!=msg.npos)
                msg = msg.substr(0, pos);

            // name, level, function name, lineno, message, args, exc_info
            PyRef<> rec(PyObject_CallMethod(logger.py(), "makeRecord", "sHsHsOO",
                                            FLAME_LOGGER_NAME, r.level, r.fname, r.lnum,
                                            msg.c_str(), Py_None, Py_None));
            PyRef<> junk(PyObject_CallMethod(logger.py(), "handle", "O", rec.py()));
        }
    }
    static void noopdtor(Machine::Logger*) {}
    static void unreg()
    {
        singleton.logger.clear();
    }

    static PyLogger singleton;
};

PyLogger PyLogger::singleton;

static
PyObject* py_set_log(PyObject *unused, PyObject *args)
{
    unsigned short lvl;
    if(!PyArg_ParseTuple(args, "H", &lvl))
        return NULL;
    // no locking here as races only result in messages being displayed (or not)
    Machine::log_detail = lvl;
    Py_RETURN_NONE;
}

static
PyObject* py_get_log(PyObject *unused)
{
    return PyString_FromString(FLAME_LOGGER_NAME);
}

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wcast-function-type"
static
PyMethodDef modmethods[] = {
    {"_GLPSParse", (PyCFunction)&PyGLPSParse, METH_VARARGS|METH_KEYWORDS,
     "Parse a GLPS lattice file to AST form"},
    {"GLPSPrinter", (PyCFunction)&PyGLPSPrint, METH_VARARGS,
     "Print a dictionary in GLPS format to string"},
    {"setLogLevel", (PyCFunction)&py_set_log, METH_VARARGS,
     "Set the FLAME logging level"
    },
    {"getLoggerName", (PyCFunction)&py_get_log, METH_NOARGS,
     "Returns the logger name used by the FLAME C extensions"
    },
    {NULL, NULL, 0, NULL}
};
#pragma GCC diagnostic pop

#if PY_MAJOR_VERSION >= 3
static struct PyModuleDef module = {
   PyModuleDef_HEAD_INIT,
   "flame._internal",
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

        try {
            PyRef<> logging(PyImport_ImportModule("logging"));
            PyLogger::singleton.logger.reset(PyObject_CallMethod(logging.py(), "getLogger", "s", FLAME_LOGGER_NAME));
            if(Py_AtExit(&PyLogger::unreg)){
                std::cerr<<"Failed to add atexit PyLogger::unreg\n";
            } else {
                boost::shared_ptr<Machine::Logger> log(&PyLogger::singleton, &PyLogger::noopdtor);
                Machine::set_logger(log);
            }
        } catch(std::runtime_error& e){
            std::cerr<<"Failed to connect flame logging to python logging : "<<typeid(e).name()<<" : "<<e.what()<<"\n";
        }

#if PY_MAJOR_VERSION >= 3
        // w/ py3 we own the module object and return NULL on import failure
        PyRef<> modref(PyModule_Create(&module));
        PyObject *mod = modref.py();
#else
        // w/ py2 we get a borrowed ref. and indicate import
        // failure by returning w/ a python exception set
        PyObject *mod = Py_InitModule("flame._internal", modmethods);
#endif

        // python API version
        PyModule_AddIntConstant(mod, "_pyapi_version", 0);
        // C API version
        PyModule_AddIntConstant(mod, "_capi_version", FLAME_API_VERSION);

        PyModule_AddIntMacro(mod, FLAME_ERROR);
        PyModule_AddIntMacro(mod, FLAME_WARN);
        PyModule_AddIntMacro(mod, FLAME_INFO);
        PyModule_AddIntMacro(mod, FLAME_DEBUG);
        PyModule_AddIntMacro(mod, FLAME_FINE);

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
