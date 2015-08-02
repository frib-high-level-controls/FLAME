
#include <string>
#include <map>

#include <boost/python/object.hpp>
#include <boost/python/class.hpp>
#include <boost/python/errors.hpp>

#define NO_IMPORT_ARRAY
#define PY_ARRAY_UNIQUE_SYMBOL USCSI_PyArray_API
#include <numpy/ndarrayobject.h>

#include "scsi/base.h"
#include "pyscsi.h"

namespace bp = boost::python;

static PyObject *statetype;

struct PyState : public boost::noncopyable
{
    std::auto_ptr<StateBase> state;
    typedef std::map<std::string, bp::object> attrs_t;
    attrs_t attrs;

    PyState() {}
    ~PyState() {}

    void attach(StateBase* b, PyObject *pyobj)
    {
        attrs_t newattrs;
        for(unsigned i=0; true; i++)
        {
            StateBase::ArrayInfo info;

            if(!b->getArray(i, info))
                break;

            npy_intp dims[5];
            memcpy(dims, info.dim, sizeof(dims));

            PyObject *obj = PyArray_New(&PyArray_Type, info.ndim, dims, NPY_DOUBLE,
                                        NULL, info.ptr, sizeof(*info.ptr), NPY_CARRAY, pyobj);
            if(!obj)
                throw std::runtime_error("Failed to wrap array");
            Py_INCREF(pyobj);

            newattrs[info.name] = bp::object(bp::handle<>(obj));
        }

        attrs.swap(newattrs);
        state.reset(b);
    }

    std::string tostring() const
    {
        std::ostringstream strm;
        state->show(strm);
        return strm.str();
    }

    static PyObject* getattro(PyObject *obj, PyObject *name)
    {
        char *sname;
#if PY_MAJOR_VERSION >= 3
        PyObject *data = PyUnicode_AsEncodedString(name, "ascii", "Encoding error:");
        if(!data)
            return NULL;
        sname = PyUnicode_AS_DATA(data);
#else
        sname = PyString_AsString(name);
#endif
        PyObject *ret = getattr(obj, sname);
#if PY_MAJOR_VERSION >= 3
        Py_DECREF(data);
#endif
        return ret;
    }

    static PyObject* getattr(PyObject *obj, char *name)
    {
        PyState *self=bp::extract<PyState*>(obj);
        attrs_t::const_iterator it = self->attrs.find(name);
        if(it!=self->attrs.end())
        {
            PyObject *ret = it->second.ptr();
            Py_INCREF(ret);
            return ret;
        }
        return PyErr_Format(PyExc_AttributeError, "Unknown attribute %s", name);
    }
};

PyObject* wrapstate(StateBase* b)
{
    if(b->pyptr) {
        Py_INCREF(b->pyptr);
    } else {
        static char fmt[1];
        b->pyptr = PyObject_CallFunction(statetype, fmt);
        if(!b->pyptr)
            return NULL;
        bp::pytype_check((PyTypeObject*)statetype, (PyObject*)b->pyptr);
        PyState *state = bp::extract<PyState*>((PyObject*)b->pyptr);
        try{
            state->attach(b, (PyObject*)b->pyptr);
        }catch(...){
            Py_DECREF(b->pyptr);
            b->pyptr = NULL;
            throw;
        }
    }
    return (PyObject*)b->pyptr;
    Py_RETURN_NONE;
}

StateBase* unwrapstate(PyObject* py)
{
    bp::pytype_check((PyTypeObject*)statetype, py);
    PyState *self=bp::extract<PyState*>(py);
    return self->state.get();
}

void registerModState(void)
{
    using namespace boost::python;

    object so = class_<PyState, boost::noncopyable>("State")
            .def("__str__", &PyState::tostring)
            //.def_readwrite("next_elem", &StateBase::next_elem)
            ;

    statetype = so.ptr();
    Py_INCREF(statetype);

    PyTypeObject *ptype = (PyTypeObject*)statetype;

    // class_ has ready'd this type
    assert((ptype->tp_flags&Py_TPFLAGS_READY)==Py_TPFLAGS_READY);

    ptype->tp_getattr = PyState::getattr;
    ptype->tp_getattro = PyState::getattro;
}
