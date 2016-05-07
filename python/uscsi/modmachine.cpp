
#include <sstream>

#include "scsi/base.h"
#include "pyscsi.h"


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
        assert(!machine->weak);

        std::auto_ptr<Config> C(PyGLPSParse2Config(raw, args, kws));

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
        } else {
            long index = PyInt_AsLong(pyindex);
            if(index<0 || (unsigned long)index>=machine->machine->size())
                return PyErr_Format(PyExc_IndexError, "Element index out of range");
            C = (*machine->machine)[index]->conf();
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
            Dict2Config(C, d);
        }
        std::auto_ptr<StateBase> state(machine->machine->allocState(C));
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
        std::auto_ptr<StateBase> tmpstate(state->clone());
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
        for(size_t i=0; i<observed.size(); i++) {
            (*machine)[i]->set_observer(NULL);
        }
    }
    void observe(size_t i, Observer *o)
    {
        if(i>=machine->size())
            throw std::runtime_error("element index out of range");
        if((*machine)[i]->observer())
            throw std::runtime_error("element already observed");
        observed.push_back(i);
        (*machine)[i]->set_observer(o);
    }
};

static
PyObject *PyMachine_propagate(PyObject *raw, PyObject *args, PyObject *kws)
{

    TRY {
        PyObject *state, *toobserv = NULL;
        unsigned long start = 0, max = (unsigned long)-1;
        const char *pnames[] = {"state", "start", "max", "observe", NULL};
        if(!PyArg_ParseTupleAndKeywords(args, kws, "O|kkO", (char**)pnames, &state, &start, &max, &toobserv))
            return NULL;

        PyStoreObserver observer;
        PyScopedObserver observing(machine->machine);

        if(toobserv) {
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
        PyObject *conf;
        const char *pnames[] = {"index", "config", NULL};
        if(!PyArg_ParseTupleAndKeywords(args, kws, "kO!|", (char**)pnames, &idx, &PyDict_Type, &conf))
            return NULL;

        if(idx>=machine->machine->size())
            return PyErr_Format(PyExc_ValueError, "invalid element index %lu", idx);

        // TODO: only allow existing elements to be changed.
        //       allow unassign
        Config newconf((*machine->machine)[idx]->conf());

        Dict2Config(newconf, conf, 3); // set depth=3 to prevent recursion

        machine->machine->reconfigure(idx, newconf);

        Py_RETURN_NONE;
    } CATCH2(std::invalid_argument, ValueError)
    CATCH()
}

static
Py_ssize_t PyMachine_len(PyObject *raw)
{
    TRY{
        return machine->machine->size();
    }CATCH1(-1)
}

static PyMethodDef PyMachine_methods[] = {
    {"conf", (PyCFunction)&PyMachine_conf, METH_VARARGS|METH_KEYWORDS,
     "Return configuration used to construct the Machine"},
    {"allocState", (PyCFunction)&PyMachine_allocState, METH_VARARGS|METH_KEYWORDS,
     "Allocate a new State based on this Machine's configuration"},
    {"propagate", (PyCFunction)&PyMachine_propagate, METH_VARARGS|METH_KEYWORDS,
     "Propagate the provided State through the simulation"},
    {"reconfigure", (PyCFunction)&PyMachine_reconfigure, METH_VARARGS|METH_KEYWORDS,
     "Change the configuration of an element."},
    {NULL, NULL, 0, NULL}
};

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
    "uscsi._internal.Machine",
    sizeof(PyMachine),
};

} // namespace

static const char pymdoc[] =
        "Machine(config, path=\"/directory/\", extra={})\n"
        "\n"
        "A Machine() the primary interface to the FLAME simulation engine.\n"
        "\n"
        "The 'config' argument may be a file-like object (with read())"
        " or buffer which will be parsed with the GLPS parser (see GLPSParser::parse)."
        "  Or it may be a dictionary.\n"
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
