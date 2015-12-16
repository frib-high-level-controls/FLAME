
#include <string>

#include "scsi/base.h"
#include "pyscsi.h"

#define NO_IMPORT_ARRAY
#define PY_ARRAY_UNIQUE_SYMBOL USCSI_PyArray_API
#include <numpy/ndarrayobject.h>

#define TRY PyState *state = (PyState*)raw; try

namespace {

struct PyState {
    PyObject_HEAD
    PyObject *dict, *weak; //  __dict__ and __weakref__
    PyObject *attrs; // lookup name to attribute index (for StateBase)
    StateBase *state;
};

static
int PyState_traverse(PyObject *raw, visitproc visit, void *arg)
{
    PyState *state = (PyState*)raw;
    Py_VISIT(state->attrs);
    Py_VISIT(state->dict);
    return 0;
}

static
int PyState_clear(PyObject *raw)
{
    PyState *state = (PyState*)raw;
    Py_CLEAR(state->dict);
    Py_CLEAR(state->attrs);
    return 0;
}

static
void PyState_free(PyObject *raw)
{
    TRY {
        std::auto_ptr<StateBase> S(state->state);
        state->state = NULL;

        if(state->weak)
            PyObject_ClearWeakRefs(raw);

        PyState_clear(raw);

        Py_TYPE(raw)->tp_free(raw);
    } CATCH2V(std::exception, RuntimeError)
}

static
PyObject *PyState_getattro(PyObject *raw, PyObject *attr)
{
    TRY {
        PyObject *val = PyDict_GetItem(state->attrs, attr);
        int i = PyInt_AsLong(val);

        StateBase::ArrayInfo info;

        if(!state->state->getArray(i, info))
            return PyErr_Format(PyExc_AttributeError, "invalid attribute name");

        npy_intp dims[5];
        memcpy(dims, info.dim, sizeof(dims));

        PyRef<> obj(PyArray_SimpleNewFromData(info.ndim, dims, NPY_DOUBLE, info.ptr));

        Py_INCREF(state);
        PyArray_BASE(obj.py()) = (PyObject*)state;

        return obj.releasePy();
    } CATCH()
}

static
int PyState_setattro(PyObject *, PyObject *, PyObject *)
{
    PyErr_SetString(PyExc_NotImplementedError, "Can't set attributes");
    return -1;
}

static
PyObject* PyState_str(PyObject *raw)
{
    TRY {
        std::ostringstream strm;
        state->state->show(strm);
        return PyString_FromString(strm.str().c_str());
    } CATCH()
}

static PyMethodDef PyState_methods[] = {
    {"__str__", (PyCFunction)&PyState_str, METH_NOARGS,
     "Render as string"},
    {NULL, NULL, 0, NULL}
};

static PyTypeObject PyStateType = {
#if PY_MAJOR_VERSION >= 3
    PyVarObject_HEAD_INIT(NULL, 0)
#else
    PyObject_HEAD_INIT(NULL)
    0,
#endif
    "uscsi._internal.State",
    sizeof(PyState),
};

} // namespace

PyObject* wrapstate(StateBase* b)
{
    try {

        PyRef<PyState> state((PyState*)PyStateType.tp_alloc(&PyStateType, 0));

        state->state = b;
        state->attrs = state->weak = state->dict = 0;

        state->attrs = PyDict_New();
        if(!state->attrs)
            return NULL;

        for(unsigned i=0; true; i++)
        {
            StateBase::ArrayInfo info;

            if(!b->getArray(i, info))
                break;

            PyRef<> name(PyInt_FromLong(i));
            if(PyDict_SetItemString(state->attrs, info.name.c_str(), name.py()))
                throw std::runtime_error("Failed to insert into Dict");

        }

        return state.releasePy();
    } CATCH()
}


StateBase* unwrapstate(PyObject* raw)
{
    if(!PyObject_TypeCheck(raw, &PyStateType))
        throw std::invalid_argument("Argument is not a State");
    PyState *state = (PyState*)raw;
    return state->state;
}

int registerModState(PyObject *mod)
{
    import_array1(-1);

    PyStateType.tp_dealloc = &PyState_free;

    PyStateType.tp_weaklistoffset = offsetof(PyState, weak);
    PyStateType.tp_traverse = &PyState_traverse;
    PyStateType.tp_clear = &PyState_clear;

    PyStateType.tp_dictoffset = offsetof(PyState, dict);
    PyStateType.tp_getattro = &PyState_getattro;
    PyStateType.tp_setattro = &PyState_setattro;

    PyStateType.tp_flags = Py_TPFLAGS_DEFAULT|Py_TPFLAGS_BASETYPE|Py_TPFLAGS_HAVE_GC;
    PyStateType.tp_methods = PyState_methods;

    if(PyType_Ready(&PyStateType))
        return -1;

    Py_INCREF(&PyStateType);
    if(PyModule_AddObject(mod, "State", (PyObject*)&PyStateType)) {
        Py_DECREF(&PyStateType);
        return -1;
    }

    return 0;
}
