#include <exception>

#ifndef PYSCSI_H
#define PYSCSI_H
#define PY_SSIZE_T_CLEAN
#include <Python.h>

struct Config;
struct StateBase;

Config* dict2conf(PyObject *dict);
void Dict2Config(Config& ret, PyObject *dict, unsigned depth=0);

PyObject* wrapstate(StateBase*); // takes ownership of argument from caller
StateBase* unwrapstate(PyObject*); // ownership of returned pointer remains with argument

PyObject* PyGLPSPrint(PyObject *, PyObject *args);
PyObject* PyGLPSParse(PyObject *, PyObject *args, PyObject *kws);

int registerModMachine(PyObject *mod);
int registerModState(PyObject *mod);

#define CATCH2V(CXX, PYEXC) catch(CXX& e) { if(!PyErr_Occurred()) PyErr_SetString(PyExc_##PYEXC, e.what()); return; }
#define CATCH3(CXX, PYEXC, RET) catch(CXX& e) { if(!PyErr_Occurred()) PyErr_SetString(PyExc_##PYEXC, e.what()); return RET; }
#define CATCH2(CXX, PYEXC) CATCH3(CXX, PYEXC, NULL)
#define CATCH1(RET) CATCH3(std::exception, RuntimeError, RET)
#define CATCH() CATCH2(std::exception, RuntimeError)

#if PY_MAJOR_VERSION >= 3
#define PyInt_Type PyLong_Type
#define PyInt_Check PyLong_Check
#define PyInt_FromLong PyLong_FromLong
#define PyInt_FromSize_t PyLong_FromSize_t
#define PyInt_AsLong PyLong_AsLong
#define PyString_Check PyUnicode_Check
#define PyString_FromString PyUnicode_FromString
#define MODINIT_RET(VAL) return (VAL)

#else
#define MODINIT_RET(VAL) return

#endif

template<typename T = PyObject>
struct PyRef {
    T* _ptr;
    PyRef() :_ptr(NULL) {}
    PyRef(T* p) : _ptr(p) {
        if(!p)
            throw std::bad_alloc(); // TODO: probably already a python exception
    }
    ~PyRef() {Py_CLEAR(_ptr);}
    T* release() {
        T* ret = _ptr;
        assert(ret);
        _ptr = NULL;
        return ret;
    }
    void clear() {
        Py_CLEAR(_ptr);
        _ptr = NULL;
    }
    void reset(T* p) {
        if(!p)
            throw std::bad_alloc(); // TODO: probably already a python exception
        Py_CLEAR(_ptr);
        _ptr = p;
    }
    struct allow_null {};
    T* reset(T* p, allow_null) {
        Py_CLEAR(_ptr);
        _ptr = p;
        return p;
    }
    PyObject* releasePy() {
        return (PyObject*)release();
    }
    PyObject* py() const {
        return (PyObject*)_ptr;
    }
    template<typename E>
    E* as() {
        return (E*)_ptr;
    }
    T& operator*() {
        return *_ptr;
    }
    T* operator->() {
        return _ptr;
    }
private:
    PyRef(const PyRef&);
    PyRef& operator =(const PyRef&);
};

#endif // PYSCSI_H
