
#include <string>
#include <sstream>

#include "flame/base.h"
#include "pyflame.h"

#define NO_IMPORT_ARRAY
#define PY_ARRAY_UNIQUE_SYMBOL FLAME_PyArray_API
#include <numpy/ndarrayobject.h>

#if SIZE_MAX==NPY_MAX_UINT32
#define NPY_SIZE_T NPY_UINT32
#elif SIZE_MAX==NPY_MAX_UINT64
#define NPY_SIZE_T NPY_UINT64
#else
#error logic error with SIZE_MAX
#endif

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
        PyObject *idx = PyDict_GetItem(state->attrs, attr);
        if(!idx) {
            return PyObject_GenericGetAttr(raw, attr);
        }
        int i = PyInt_AsLong(idx);


        StateBase::ArrayInfo info;

        if(!state->state->getArray(i, info))
            return PyErr_Format(PyExc_RuntimeError, "invalid attribute name (sub-class forgot %d)", i);

        if(info.ndim==0) { // Scalar
            switch(info.type) {
            case StateBase::ArrayInfo::Double:
                return PyFloat_FromDouble(*(double*)info.ptr);
            case StateBase::ArrayInfo::Sizet:
                return PyLong_FromSize_t(*(size_t*)info.ptr);
            }
            return PyErr_Format(PyExc_TypeError, "unsupported type code %d", info.type);
        } else if(info.ndim>2) {
            return PyErr_Format(PyExc_ValueError, "parameters with ndim>2 not supported %s has %u",
                                info.name, info.ndim);
        }

        int pytype;
        switch(info.type) {
        case StateBase::ArrayInfo::Double: pytype = NPY_DOUBLE; break;
        case StateBase::ArrayInfo::Sizet: pytype = NPY_SIZE_T; break;
        default:
            return PyErr_Format(PyExc_TypeError, "unsupported type code %d", info.type);
        }

        npy_intp dims[StateBase::ArrayInfo::maxdims];
        std::copy(info.dim,
                  info.dim+StateBase::ArrayInfo::maxdims,
                  dims);

        // Alloc new array and copy in

        PyRef<> obj(PyArray_SimpleNew(info.ndim, dims, pytype));

        if(info.ndim==1) {
            for(size_t i=0; i<info.dim[0]; i++) {
                void *dest = PyArray_GETPTR1(obj.py(), i);
                const void *src  = ((char*)info.ptr) + info.stride[0]*i;
                switch(info.type) {
                case StateBase::ArrayInfo::Double: *(double*)dest = *(double*)src; break;
                case StateBase::ArrayInfo::Sizet:  *(size_t*)dest = *(size_t*)src; break;
                }
            }
        } else if(info.ndim==2) {
            for(size_t i=0; i<info.dim[0]; i++) {
                for(size_t j=0; j<info.dim[1]; j++) {
                    void *dest = PyArray_GETPTR2(obj.py(), i, j);
                    const void *src  = ((char*)info.ptr) + info.stride[0]*i + info.stride[1]*j;
                    switch(info.type) {
                    case StateBase::ArrayInfo::Double: *(double*)dest = *(double*)src; break;
                    case StateBase::ArrayInfo::Sizet:  *(size_t*)dest = *(size_t*)src; break;
                    }
                }
            }
        }

        return obj.releasePy();
    } CATCH()
}

static
int PyState_setattro(PyObject *raw, PyObject *attr, PyObject *val)
{
    TRY {
        PyObject *idx = PyDict_GetItem(state->attrs, attr);
        if(!idx)
            return PyObject_GenericSetAttr(raw, attr, val);
        int i = PyInt_AsLong(idx);

        StateBase::ArrayInfo info;

        if(!state->state->getArray(i, info)) {
            PyErr_Format(PyExc_RuntimeError, "invalid attribute name (sub-class forgot %d)", i);
            return -1;
        }

        if(info.ndim==0) {
            switch(info.type) {
            case StateBase::ArrayInfo::Double: {
                double *dest = (double*)info.ptr;
                if(PyFloat_Check(val))
                    *dest = PyFloat_AsDouble(val);
                else if(PyLong_Check(val))
                    *dest = PyLong_AsDouble(val);
                else if(PyInt_Check(val))
                    *dest = PyInt_AsLong(val);
                else
                    PyErr_Format(PyExc_ValueError, "Can't assign to double field");
            }
                break;
            case StateBase::ArrayInfo::Sizet: {
                size_t *dest = (size_t*)info.ptr;
                if(PyFloat_Check(val))
                    *dest = PyFloat_AsDouble(val);
                else if(PyLong_Check(val))
                    *dest = PyLong_AsUnsignedLongLong(val);
                else if(PyInt_Check(val))
                    *dest = PyInt_AsLong(val);
                else
                    PyErr_Format(PyExc_ValueError, "Can't assign to double field");
            }
                break;
            default:
                PyErr_Format(PyExc_TypeError, "unsupported type code %d", info.type);
            }

            return PyErr_Occurred() ? -1 : 0;
        }

        int pytype;
        switch(info.type) {
        case StateBase::ArrayInfo::Double: pytype = NPY_DOUBLE; break;
        case StateBase::ArrayInfo::Sizet: pytype = NPY_SIZE_T; break;
        default:
            PyErr_Format(PyExc_TypeError, "unsupported type code %d", info.type);
            return -1;
        }

        PyRef<> arr(PyArray_FromObject(val, pytype, 1, 2));

        if(info.ndim!=(size_t)PyArray_NDIM(arr.py())) {
            PyErr_Format(PyExc_ValueError, "dimension don't match");
            return -1;
        }

        if(info.ndim==1) {
            for(size_t i=0; i<info.dim[0]; i++) {
                const void *src = PyArray_GETPTR1(arr.py(), i);
                void *dest  = ((char*)info.ptr) + info.stride[0]*i;
                switch(info.type) {
                case StateBase::ArrayInfo::Double: *(double*)dest = *(double*)src; break;
                case StateBase::ArrayInfo::Sizet:  *(size_t*)dest = *(size_t*)src; break;
                }
            }
        } else if(info.ndim==2) {
            for(size_t i=0; i<info.dim[0]; i++) {
                for(size_t j=0; j<info.dim[1]; j++) {
                    const void *src = PyArray_GETPTR2(arr.py(), i, j);
                    void *dest  = ((char*)info.ptr) + info.stride[0]*i + info.stride[1]*j;
                    switch(info.type) {
                    case StateBase::ArrayInfo::Double: *(double*)dest = *(double*)src; break;
                    case StateBase::ArrayInfo::Sizet:  *(size_t*)dest = *(size_t*)src; break;
                    }
                }
            }
        }

        return 0;
    } CATCH3(std::exception, RuntimeError, -1)
}

static
PyObject* PyState_str(PyObject *raw)
{
    TRY {
        std::ostringstream strm;
        state->state->show(strm, 0);
        return PyString_FromString(strm.str().c_str());
    } CATCH()
}

static
PyObject* PyState_iter(PyObject *raw)
{
    TRY {
        return PyObject_GetIter(state->attrs);
    }CATCH()
}

static
Py_ssize_t PyState_len(PyObject *raw)
{
    TRY{
        return PyObject_Length(state->attrs);
    }CATCH1(-1)
}

static PySequenceMethods PyState_seq = {
    &PyState_len
};

static
PyObject* PyState_clone(PyObject *raw, PyObject *unused)
{
    TRY {
        std::auto_ptr<StateBase> newstate(state->state->clone());

        PyObject *ret = wrapstate(newstate.get());
        newstate.release();
        return ret;
    } CATCH()
}

static PyMethodDef PyState_methods[] = {
    {"clone", (PyCFunction)&PyState_clone, METH_NOARGS,
     "clone()\n\n"
     "Returns a new State instance which is a copy of this one"
    },
    {NULL, NULL, 0, NULL}
};

static PyTypeObject PyStateType = {
#if PY_MAJOR_VERSION >= 3
    PyVarObject_HEAD_INIT(NULL, 0)
#else
    PyObject_HEAD_INIT(NULL)
    0,
#endif
    "flame._internal.State",
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
            if(PyDict_SetItemString(state->attrs, info.name, name.py()))
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

static const char pymdoc[] =
        "The interface to a sub-class of C++ StateBase.\n"
        "Can't be constructed from python, see Machine.allocState()\n"
        "\n"
        "Provides access to some C++ member variables via the Machine::getArray() interface.\n"
        ;

int registerModState(PyObject *mod)
{
    PyStateType.tp_doc = pymdoc;

    PyStateType.tp_str = &PyState_str;
    PyStateType.tp_repr = &PyState_str;
    PyStateType.tp_dealloc = &PyState_free;

    PyStateType.tp_iter = &PyState_iter;
    PyStateType.tp_as_sequence = &PyState_seq;

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
