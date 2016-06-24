
#include <list>
#include <sstream>

#include "flame/base.h"

#include "pyflame.h"

#define NO_IMPORT_ARRAY
#define PY_ARRAY_UNIQUE_SYMBOL FLAME_PyArray_API
#include <numpy/ndarrayobject.h>

/** Translate python list of tuples to Config
 *
 *  [(K,V)]   -> Config
 *    float   -> double
 *    str     -> string
 *    [[()]]    -> vector<Config>  (recurse)
 *    ndarray -> vector<double>
 *    TODO: [0.0]   -> vector<double>
 */
void List2Config(Config& ret, PyObject *list, unsigned depth)
{
    if(depth>3)
        throw std::runtime_error("too deep for Dict2Config");

    PyRef<> iter(PyObject_GetIter(list));
    while(true) {
        PyObject *item = PyIter_Next(iter.py());
        if(!item) break;
        PyRef<> itemref(item);
        const char *kname;
        PyObject *value;
        if(!PyArg_ParseTuple(item, "sO", &kname, &value))
            throw std::runtime_error("list item is not a tuple?");

        if(PyArray_Check(value)) { // array as vector<double>
            PyRef<> arr(PyArray_ContiguousFromAny(value, NPY_DOUBLE, 0, 2));
            double *buf = (double*)PyArray_DATA(arr.py());
            std::vector<double> temp(PyArray_SIZE(arr.py()));
            std::copy(buf, buf+temp.size(), temp.begin());

            ret.swap<std::vector<double> >(kname, temp);

        } else if(PyNumber_Check(value)) { // scalar as double
            PyRef<> dval(PyNumber_Float(value));
            double val = PyFloat_AsDouble(dval.py());
            ret.set<double>(kname, val);

        } else if(PyUnicode_Check(value) || (PY_MAJOR_VERSION < 3 && PyBytes_Check(value))) { // string
            PyRef<> valref(value, borrow());
            PyCString sval(valref);
            const char *val = sval.c_str();

            ret.set<std::string>(kname, val);

        } else if(PySequence_Check(value)) { // list of dict
            Py_ssize_t N = PySequence_Size(value);

            Config::vector_t output;
            output.reserve(N);

            for(Py_ssize_t i=0; i<N; i++) {
                PyObject *elem = PySequence_GetItem(value, i);
                assert(elem);

                PyRef<> list;
                if(PyDict_Check(elem)) {
                    list.reset(PyMapping_Items(elem));
                    elem = list.py();
                }
                if(!PyList_Check(elem)) {
                    PyTypeObject *valuetype = (PyTypeObject*)PyObject_Type(elem);
                    throw std::invalid_argument(SB()<<"lists must contain only dict or list of tuples, not "<<valuetype->tp_name);
                }

                output.push_back(ret.new_scope());

                List2Config(output.back(), elem, depth+1); // inheirt parent scope
            }

            ret.set<Config::vector_t>(kname, output);

        } else {
            PyTypeObject *valuetype = (PyTypeObject*)PyObject_Type(value);
            throw std::invalid_argument(SB()<<"Must be a dict, not "<<valuetype->tp_name);
        }
    }
}

namespace {

PyObject* pydirname(PyObject *obj)
{
    if(!obj) return NULL;

    PyRef<> ospath(PyImport_ImportModule("os.path"));
    return PyObject_CallMethod(ospath.py(), "dirname", "O", obj);
}

struct confval : public boost::static_visitor<PyRef<> >
{
    PyRef<> operator()(double v) const
    {
        return PyRef<>(PyFloat_FromDouble(v));
    }

    PyRef<> operator()(const std::string& v) const
    {
        return PyRef<>(PyString_FromString(v.c_str()));
    }

    PyRef<> operator()(const std::vector<double>& v) const
    {
        npy_intp dims[]  = {(npy_intp)v.size()};
        PyRef<> obj(PyArray_SimpleNew(1, dims, NPY_DOUBLE));
        std::copy(v.begin(), v.end(), (double*)PyArray_DATA(obj.py()));
        return obj;
    }

    PyRef<> operator()(const Config::vector_t& v) const
    {
        PyRef<> L(PyList_New(v.size()));

        for(size_t i=0, N=v.size(); i<N; i++)
        {
            PyList_SetItem(L.py(), i, conf2dict(&v[i]));
        }
        return L;
    }
};

} // namespace

PyObject* conf2dict(const Config *conf)
{
    PyRef<> list(PyList_New(0));

    for(Config::const_iterator it=conf->begin(), end=conf->end();
        it!=end; ++it)
    {
        PyRef<> val(boost::apply_visitor(confval(), it->second));
        PyRef<> tup(Py_BuildValue("sO", it->first.c_str(), val.py()));
        if(PyList_Append(list.py(), tup.py()))
            throw std::runtime_error("Failed to insert into dictionary from conf2dict");
    }

    return list.releasePy();
}

Config* list2conf(PyObject *dict)
{
    if(!PyList_Check(dict))
        throw std::invalid_argument("Not a list");
    std::auto_ptr<Config> conf(new Config);
    List2Config(*conf, dict);
    return conf.release();
}

PyObject* PyGLPSPrint(PyObject *, PyObject *args)
{
    try {
        PyObject *inp;
        if(!PyArg_ParseTuple(args, "O", &inp))
            return NULL;

        PyRef<> list;
        if(PyDict_Check(inp)) {
            list.reset(PyMapping_Items(inp));
            inp = list.py();
        }
        if(!PyList_Check(inp))
            return PyErr_Format(PyExc_ValueError, "argument must be dict or list of tuples");

        Config conf;
        List2Config(conf, inp);
        std::ostringstream strm;
        GLPSPrint(strm, conf);
        return PyString_FromString(strm.str().c_str());
    }CATCH()
}

#ifndef PY_SSIZE_T_CLEAN
#error the following assumes ssize_t is used
#endif

Config* PyGLPSParse2Config(PyObject *, PyObject *args, PyObject *kws)
{
    PyObject *conf = NULL, *extra_defs = Py_None;
    const char *path = NULL;
    const char *pnames[] = {"config", "path", "extra", NULL};
    if(!PyArg_ParseTupleAndKeywords(args, kws, "O|zO", (char**)pnames, &conf, &path, &extra_defs))
        return NULL;

    GLPSParser parser;

    if(extra_defs==Py_None) {
        // no-op
    } else if(PyDict_Check(extra_defs)) {
        PyObject *key, *value;
        Py_ssize_t pos = 0;

        while(PyDict_Next(extra_defs, &pos, &key, &value)) {
            PyRef<> keyx(key, borrow());
            PyCString keystr(keyx);

            Config::value_t curval;

            if(PyNumber_Check(value)) {
                PyRef<> pyf(PyNumber_Float(value));
                curval = PyFloat_AsDouble(pyf.py());

            } else if(PyString_Check(value)) {
                PyRef<> valuex(value, borrow());
                PyCString valstr(valuex);

                curval = valstr.c_str();

            } else {
                PyErr_SetString(PyExc_ValueError, "extra {} can contain only numbers or strings");
                return NULL;
            }

            parser.setVar(keystr.c_str(), curval);
        }
    } else {
        PyErr_SetString(PyExc_ValueError, "'extra' must be a dict");
        return NULL;
    }

    PyGetBuf buf;
    std::auto_ptr<Config> C;

    PyRef<> listref;

    if(PyObject_HasAttrString(conf, "read")) { // file-like
        PyCString pyname;

        if(!path && PyObject_HasAttrString(conf, "name")) {
            path = pyname.c_str(pydirname(PyObject_GetAttrString(conf, "name")));
        }

        PyRef<> pybytes(PyObject_CallMethod(conf, "read", ""));
        if(!buf.get(pybytes.py())) {
            PyErr_SetString(PyExc_TypeError, "read() must return a buffer");
            return NULL;
        }
        C.reset(parser.parse_byte((const char*)buf.data(), buf.size(), path));

    } else if(buf.get(conf)) {
        C.reset(parser.parse_byte((const char*)buf.data(), buf.size(), path));

#if PY_MAJOR_VERSION >= 3
    } else if(PyUnicode_Check(conf)) { // py3 str (aka unicode) doesn't implement buffer iface
        PyCString buf;
        const char *cbuf = buf.c_str(conf);

        C.reset(parser.parse_byte(cbuf, strlen(cbuf), path));
#endif

    } else {
        if(PyDict_Check(conf)) {
            listref.reset(PyMapping_Items(conf));
            conf = listref.py();
        }
        if(PyList_Check(conf)) {
            C.reset(list2conf(conf));

        } else {
            throw std::invalid_argument("'config' must be dict, list of tuples, or byte buffer");
        }
    }

    return C.release();
}

PyObject* PyGLPSParse(PyObject *unused, PyObject *args, PyObject *kws)
{
    try{
        return conf2dict(PyGLPSParse2Config(unused, args, kws));
    }CATCH()
}
