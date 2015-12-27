
#include <sstream>

#include "scsi/base.h"
#include "pyscsi.h"

#define NO_IMPORT_ARRAY
#define PY_ARRAY_UNIQUE_SYMBOL USCSI_PyArray_API
#include <numpy/ndarrayobject.h>

#define TRY PyMachine *machine = (PyMachine*)raw; try

namespace {

struct PyMachine {
    PyObject_HEAD

    PyObject *weak;
    Machine *machine;
};

static
int PyMachine_init(PyObject *raw, PyObject *args, PyObject *kws)
{
    TRY {
        PyObject *conf;
        const char *pnames[] = {"config", NULL};
        if(!PyArg_ParseTupleAndKeywords(args, kws, "|O!", (char**)pnames, &PyDict_Type, &conf))
            return -1;

        assert(!machine->weak);

        std::auto_ptr<Config> C(dict2conf(conf));
        machine->machine = new Machine(*C);

        return 0;
    } CATCH3(key_error, KeyError, -1)
      CATCH3(std::invalid_argument, ValueError, -1)
      CATCH3(std::exception, RuntimeError, -1)
}

static
void PyMachine_free(PyObject *raw)
{
    TRY {
        std::auto_ptr<Machine> S(machine->machine);
        machine->machine = NULL;

        if(machine->weak)
            PyObject_ClearWeakRefs(raw);

        Py_TYPE(raw)->tp_free(raw);
    } CATCH2V(std::exception, RuntimeError)
}

static
PyObject *PyMachine_str(PyObject *raw)
{
    TRY {
        std::ostringstream strm;
        strm << *(machine->machine);
        return PyString_FromString(strm.str().c_str());
    } CATCH()
}

static
PyObject *PyMachine_allocState(PyObject *raw, PyObject *args, PyObject *kws)
{
    TRY {
        PyObject *d;
        const char *pnames[] = {"config", NULL};
        if(!PyArg_ParseTupleAndKeywords(args, kws, "O|", (char**)pnames, &d))
            return NULL;

        std::auto_ptr<Config> C(dict2conf(d));
        std::auto_ptr<StateBase> state(machine->machine->allocState(*C));
        PyObject *ret = wrapstate(state.get());
        state.release();
        return ret;
    } CATCH()
}

static
PyObject *PyMachine_propagate(PyObject *raw, PyObject *args, PyObject *kws)
{

    TRY {
        PyObject *state;
        unsigned long start = 0, max = (unsigned long)-1;
        const char *pnames[] = {"state", "start", "max", NULL};
        if(!PyArg_ParseTupleAndKeywords(args, kws, "O|kk", (char**)pnames, &state, &start, &max))
            return NULL;

        machine->machine->propagate(unwrapstate(state), start, max);
        Py_RETURN_NONE;
    } CATCH2(std::invalid_argument, ValueError)
    CATCH()
}

static PyMethodDef PyMachine_methods[] = {
    {"allocState", (PyCFunction)&PyMachine_allocState, METH_VARARGS|METH_KEYWORDS,
     "Allocate a new State based on this Machine's configuration"},
    {"propagate", (PyCFunction)&PyMachine_propagate, METH_VARARGS|METH_KEYWORDS,
     "Propagate the provided State through the simulation"},
    {NULL, NULL, 0, NULL}
};

static PyTypeObject PyMachineType = {
#if PY_MAJOR_VERSION >= 3
    PyVarObject_HEAD_INIT(NULL, 0)
#else
    PyObject_HEAD_INIT(NULL)
    0,
#endif
    "uscsi._internal.Machine",
    sizeof(PyMachine),
};

} // namespace

static const char pymdoc[] =
        "Simulation execution engine.\n"
        "\n"
        "A Machine object is responsible for carrying out calculations\n"
        "based on a Config provided when it was constructed.\n"
        "\n"
        "See the allocState and propagate methods.\n"
        ;

int registerModMachine(PyObject *mod)
{
    PyMachineType.tp_doc = pymdoc;
    PyMachineType.tp_str = &PyMachine_str;

    PyMachineType.tp_new = &PyType_GenericNew;
    PyMachineType.tp_init = &PyMachine_init;
    PyMachineType.tp_dealloc = &PyMachine_free;

    PyMachineType.tp_weaklistoffset = offsetof(PyMachine, weak);

    PyMachineType.tp_flags = Py_TPFLAGS_DEFAULT|Py_TPFLAGS_BASETYPE;
    PyMachineType.tp_methods = PyMachine_methods;

    if(PyType_Ready(&PyMachineType))
        return -1;

    Py_INCREF(&PyMachineType);
    if(PyModule_AddObject(mod, "Machine", (PyObject*)&PyMachineType)) {
        Py_DECREF(&PyMachineType);
        return -1;
    }

    return 0;
}
