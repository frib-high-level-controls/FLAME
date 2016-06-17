#include <stdexcept>

#ifndef PYFLAME_H
#define PYFLAME_H
#define PY_SSIZE_T_CLEAN
#include <Python.h>

struct Config;
struct StateBase;

Config* dict2conf(PyObject *dict);
void Dict2Config(Config& ret, PyObject *dict, unsigned depth=0);
PyObject* conf2dict(const Config *conf);

PyObject* wrapstate(StateBase*); // takes ownership of argument from caller
StateBase* unwrapstate(PyObject*); // ownership of returned pointer remains with argument

PyObject* PyGLPSPrint(PyObject *, PyObject *args);
Config* PyGLPSParse2Config(PyObject *, PyObject *args, PyObject *kws);
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

struct borrow {};

template<typename T = PyObject>
struct PyRef {
    T* _ptr;
    PyRef() :_ptr(NULL) {}
    //! copy an existing reference (calls Py_XINCREF)
    PyRef(const PyRef& o)
        :_ptr(o._ptr)
    {
        Py_XINCREF(_ptr);
    }
    //! Construct around an existing reference (does not call Py_XINCREF).
    //! Explicitly take ownership of a reference which must already
    //! be implicitly owned by the caller.
    //! @throws std::bad_allo if p==NULL (assumed to be a python exception)
    explicit PyRef(T* p) : _ptr(p) {
        if(!p)
            throw std::bad_alloc(); // TODO: probably already a python exception
    }
    //! Construct around a borrowed reference (calls Py_XINCREF)
    PyRef(T* p, borrow) : _ptr(p) {
        if(!p)
            throw std::bad_alloc();
        Py_INCREF(p);
    }
    ~PyRef() {Py_CLEAR(_ptr);}
    PyRef& operator =(const PyRef& o)
    {
        if(&o!=this) {
            PyObject *tmp = _ptr;
            _ptr = o._ptr;
            Py_XINCREF(_ptr);
            Py_XDECREF(tmp);
        }
        return *this;
    }
    T* release() {
        T* ret = _ptr;
        assert(ret);
        _ptr = NULL;
        return ret;
    }
    void clear() {
        Py_CLEAR(_ptr);
    }
    void reset(T* p) {
        if(!p)
            throw std::bad_alloc(); // TODO: probably already a python exception
        Py_CLEAR(_ptr);
        _ptr = p;
    }
    void reset(T* p, borrow) {
        if(!p)
            throw std::bad_alloc(); // TODO: probably already a python exception
        PyObject *tmp = _ptr;
        _ptr = p;
        Py_INCREF(p);
        Py_XDECREF(tmp);
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
    E* as() const {
        return (E*)_ptr;
    }
    T& operator*() const {
        return *_ptr;
    }
    T* operator->() const {
        return _ptr;
    }
};

// Extract C string from python object (py2 str or py3 unicode)
struct PyCString
{
#if PY_MAJOR_VERSION >= 3
    PyRef<> ascii;
#else
    PyRef<> pystr;
#endif
    PyCString() {}
    PyCString(PyRef<>& o)
#if PY_MAJOR_VERSION >= 3
        :ascii(PyUnicode_AsASCIIString(o.py()))
#else
        :pystr(o)
#endif
    {}
    const char *c_str(PyObject *obj) {
        if(!obj)
            throw std::bad_alloc();
#if PY_MAJOR_VERSION >= 3
        ascii.reset(PyUnicode_AsASCIIString(obj));
#else
        pystr.reset(obj, borrow());
#endif
        return c_str();
    }
    const char *c_str() const {
        const char *ret;
#if PY_MAJOR_VERSION >= 3
        ret = PyBytes_AsString(ascii.py());
#else
        ret = PyString_AsString(pystr.py());
#endif
        if(!ret)
            throw std::invalid_argument("Can't extract string from object");
        return ret;

    }
private:
    PyCString(const PyCString&);
    PyCString& operator=(const PyCString&);
};

struct PyGetBuf
{
    Py_buffer buf;
    bool havebuf;
    PyGetBuf() : havebuf(false) {}
    ~PyGetBuf() {
        if(havebuf) PyBuffer_Release(&buf);
    }
    bool get(PyObject *obj) {
        if(havebuf) PyBuffer_Release(&buf);
        havebuf = PyObject_GetBuffer(obj, &buf, PyBUF_SIMPLE)==0;
        if(!havebuf) PyErr_Clear();
        return havebuf;
    }
    size_t size() const { return havebuf ? buf.len : 0; }
    void * data() const { return buf.buf; }
};

#endif // PYFLAME_H
