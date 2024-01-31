#include <climits>
#include <sstream>

#include "flame/base.h"
#include "pyflame.h"


#define TRY PyMachine *machine = reinterpret_cast<PyMachine*>(raw); try

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
        assert(!machine->weak);

        std::unique_ptr<Config> C(PyGLPSParse2Config(raw, args, kws));

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
        std::unique_ptr<Machine> S(machine->machine);
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
PyObject *PyMachine_conf(PyObject *raw, PyObject *args, PyObject *kws)
{
    TRY {
        PyObject *pyindex = Py_None;
        const char *pnames[] = {"index", NULL};
        if(!PyArg_ParseTupleAndKeywords(args, kws, "|O", (char**)pnames, &pyindex))
            return NULL;

        Config C;
        if(pyindex==Py_None) {
            C = machine->machine->conf();
        } else if(PyNumber_Check(pyindex)) {
            PyRef<> pylong(PyNumber_Long(pyindex));
            long index = PyLong_AsLong(pylong.py());
            if(index<0 || (unsigned long)index>=machine->machine->size())
                return PyErr_Format(PyExc_IndexError, "Element index out of range");
            C = (*machine->machine)[index]->conf();
        } else {
            return PyErr_Format(PyExc_ValueError, "'index' must be an integer or None");
        }
        C.flatten();

        return conf2dict(&C);
    } CATCH()
}

static
PyObject *PyMachine_allocState(PyObject *raw, PyObject *args, PyObject *kws)
{
    TRY {
        PyObject *d = Py_None, *W = Py_False;
        const char *pnames[] = {"config", "inherit", NULL};
        if(!PyArg_ParseTupleAndKeywords(args, kws, "|OO", (char**)pnames, &d, &W))
            return NULL;

        Config C;
        if(d==Py_None) {
            C = machine->machine->conf();
        } else if(PyDict_Check(d)) {
            if(PyObject_IsTrue(W)) {
                C = machine->machine->conf();
                C.push_scope();
            }
            PyRef<> list(PyMapping_Items(d));
            List2Config(C, list.py());
        } else {
            return PyErr_Format(PyExc_ValueError, "allocState() needs config=None or {}");
        }
        std::unique_ptr<StateBase> state(machine->machine->allocState(C));
        PyObject *ret = wrapstate(state.get());
        state.release();
        return ret;
    } CATCH()
}

struct PyStoreObserver : public Observer
{
    PyRef<> list;
    PyStoreObserver()
        :list(PyList_New(0))
    {}
    virtual ~PyStoreObserver() {}
    virtual void view(const ElementVoid* elem, const StateBase* state)
    {
        PyRef<> tuple(PyTuple_New(2));
        std::unique_ptr<StateBase> tmpstate(state->clone());
        PyRef<> statecopy(wrapstate(tmpstate.get()));
        tmpstate.release();

        PyTuple_SET_ITEM(tuple.py(), 0, PyInt_FromSize_t(elem->index));
        PyTuple_SET_ITEM(tuple.py(), 1, statecopy.release());
        if(PyList_Append(list.py(), tuple.py()))
            throw std::runtime_error(""); // a py exception is active
    }
};

struct PyScopedObserver
{
    Machine *machine;
    std::vector<size_t> observed;
    PyScopedObserver(Machine *m) : machine(m) {}
    ~PyScopedObserver() {
        for(size_t i=0; i<machine->size(); i++) {
            (*machine)[i]->set_observer(NULL);
        }
    }
    void observe(size_t i, Observer *o)
    {
        if(i>=machine->size())
            throw std::runtime_error("element index out of range");
        observed.push_back(i);
        (*machine)[i]->set_observer(o);
    }
};

static
PyObject *PyMachine_propagate(PyObject *raw, PyObject *args, PyObject *kws)
{

    TRY {
        PyObject *state, *toobserv = Py_None, *pymax = Py_None;
        unsigned long start = 0;
        int max = INT_MAX;
        const char *pnames[] = {"state", "start", "max", "observe", NULL};
        if(!PyArg_ParseTupleAndKeywords(args, kws, "O|kOO", (char**)pnames, &state, &start, &pymax, &toobserv))
            return NULL;

        if (pymax!=Py_None) max = (int) PyLong_AsLong(pymax);

        PyStoreObserver observer;
        PyScopedObserver observing(machine->machine);

        if(toobserv!=Py_None) {
            PyRef<> iter(PyObject_GetIter(toobserv)), item;

            while(item.reset(PyIter_Next(iter.py()), PyRef<>::allow_null())) {
                Py_ssize_t num = PyNumber_AsSsize_t(item.py(), PyExc_ValueError);
                if(PyErr_Occurred())
                    throw std::runtime_error(""); // caller will get active python exception
                observing.observe(num, &observer);

            }
        }

        machine->machine->propagate(unwrapstate(state), start, max);
        if(toobserv) {
            return observer.list.release();
        } else {
            Py_RETURN_NONE;
        }
    } CATCH2(std::invalid_argument, ValueError)
    CATCH()
}

static
PyObject *PyMachine_reconfigure(PyObject *raw, PyObject *args, PyObject *kws)
{
    TRY{
        unsigned long idx;
        PyObject *conf, *replace = Py_False;
        const char *pnames[] = {"index", "config", "replace", NULL};
        if(!PyArg_ParseTupleAndKeywords(args, kws, "kO!|O", (char**)pnames, &idx, &PyDict_Type, &conf, &replace))
            return NULL;

        if(idx>=machine->machine->size())
            return PyErr_Format(PyExc_ValueError, "invalid element index %lu", idx);

        Config newconf;
        if(!PyObject_IsTrue(replace))
            newconf = (*machine->machine)[idx]->conf();

        PyRef<> list(PyMapping_Items(conf));
        List2Config(newconf, list.py(), 3); // set depth=3 to prevent recursion

        machine->machine->reconfigure(idx, newconf);

        Py_RETURN_NONE;
    } CATCH2(std::invalid_argument, ValueError)
    CATCH()
}

static
PyObject *PyMachine_find(PyObject *raw, PyObject *args, PyObject *kws)
{
    TRY{
        const char *ename = NULL, *etype = NULL;
        const char *pnames[] = {"name", "type", NULL};
        if(!PyArg_ParseTupleAndKeywords(args, kws, "|zz", (char**)pnames, &ename, &etype))
            return NULL;

        PyRef<> ret(PyList_New(0));

        std::pair<Machine::lookup_iterator, Machine::lookup_iterator> range;

        if(ename && etype) {
            return PyErr_Format(PyExc_ValueError, "only one of 'name' or 'type' may be given");
        } else if(ename) {
            range = machine->machine->equal_range(ename);
        } else if(etype) {
            range = machine->machine->equal_range_type(etype);
        } else {
            range = machine->machine->all_range();
        }

        for(; range.first!=range.second; ++range.first) {
            ElementVoid *elem = *range.first;

            PyRef<> pyidx(PyInt_FromLong(elem->index));

            if(PyList_Append(ret.py(), pyidx.py()))
                return NULL;
        }

        return ret.release();
    }CATCH()
}

static
Py_ssize_t PyMachine_len(PyObject *raw)
{
    TRY{
        return machine->machine->size();
    }CATCH1(-1)
}

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wcast-function-type"
static PyMethodDef PyMachine_methods[] = {
    {"conf", (PyCFunction)&PyMachine_conf, METH_VARARGS|METH_KEYWORDS,
     "conf() -> {} Machine config\n"
     "conf(index) -> {} Element config"},
    {"allocState", (PyCFunction)&PyMachine_allocState, METH_VARARGS|METH_KEYWORDS,
     "allocState() -> State\n"
     "allocState({'variable':int|str}) -> State\n"
     "Allocate a new State based on this Machine's configuration."
     "  Optionally provide additional configuration"},
    {"propagate", (PyCFunction)&PyMachine_propagate, METH_VARARGS|METH_KEYWORDS,
     "propagate(State, start=0, max=INT_MAX, observe=None)\n"
     "propagate(State, start=0, max=INT_MAX, observe=[1,4,...]) -> [(index,State), ...]\n"
     "Propagate the provided State through the simulation.\n"
     "\n"
     "start and max selects through which element the State will be passed.\n"
     "\n"
     "observe may be None or an iterable yielding element indicies.\n"
     "In the second form propagate() returns a list of tuples with the output State of the selected elements."
    },
    {"reconfigure", (PyCFunction)&PyMachine_reconfigure, METH_VARARGS|METH_KEYWORDS,
     "reconfigure(index, {'variable':int|str})\n"
     "Change the configuration of an element."},
    {"find", (PyCFunction)&PyMachine_find, METH_VARARGS|METH_KEYWORDS,
    "find(name=None, type=None) -> [int]\n"
    "Return a list of element indices for element name or type matching the given string."},
    {NULL, NULL, 0, NULL}
};
#pragma GCC diagnostic pop

static PySequenceMethods PyMachine_seq = {
    &PyMachine_len
};

static PyTypeObject PyMachineType = {
#if PY_MAJOR_VERSION >= 3
    PyVarObject_HEAD_INIT(NULL, 0)
#else
    PyObject_HEAD_INIT(NULL)
    0,
#endif
    "flame._internal.Machine",
    sizeof(PyMachine),
};

} // namespace

static const char pymdoc[] =
        "Machine(config, path=None, extra=None)\n"
        "Machine(config, path='/directry/', extra={'variable':float|str}})\n"
        "\n"
        "A Machine() the primary interface to the FLAME simulation engine.\n"
        "\n"
        "The 'config' argument may be a file-like object (with read())"
        " or a buffer which will be parsed with the GLPS parser (see GLPSParser::parse).\n"
        " Or it may be a dictionary.\n"
        "\n"
        ">>> with open('some.lat', 'rb') as F:\n"
        "      M = Machine(F)\n"
        ">>>\n"
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
    PyMachineType.tp_as_sequence = &PyMachine_seq;

    if(PyType_Ready(&PyMachineType))
        return -1;

    Py_INCREF((PyObject*)&PyMachineType);
    if(PyModule_AddObject(mod, "Machine", (PyObject*)&PyMachineType)) {
        Py_DECREF((PyObject*)&PyMachineType);
        return -1;
    }

    return 0;
}
