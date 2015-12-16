#ifndef PYSCSI_H
#define PYSCSI_H

#include <Python.h>

struct Config;
struct StateBase;

Config* dict2conf(PyObject *dict);

PyObject* wrapstate(StateBase*); // takes ownership of argument from caller
StateBase* unwrapstate(PyObject*); // ownership of returned pointer remains with argument

PyObject* PyGLPSPrint(PyObject *, PyObject *args);
PyObject* PyGLPSParse(PyObject *, PyObject *args);

int registerModMachine(PyObject *mod);
int registerModState(PyObject *mod);

#define CATCH2V(CXX, PYEXC) catch(CXX& e) { std::cerr<<"Exception during "<<__PRETTY_FUNCTION__<<" :"<<e.what()<<"\n"; return; }
#define CATCH3(CXX, PYEXC, RET) catch(CXX& e) { PyErr_SetString(PyExc_##PYEXC, e.what()); return RET; }
#define CATCH2(CXX, PYEXC) CATCH3(CXX, PYEXC, NULL)
#define CATCH() CATCH2(std::exception, RuntimeError)

#if PY_MAJOR_VERSION >= 3
#define PyInt_FromLong PyLong_FromLong
#define PyInt_AsLong PyLong_AsLong
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
