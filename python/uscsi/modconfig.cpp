
#include <list>
#include <sstream>

#include "scsi/base.h"

#include "pyscsi.h"

#define NO_IMPORT_ARRAY
#define PY_ARRAY_UNIQUE_SYMBOL USCSI_PyArray_API
#include <numpy/ndarrayobject.h>

namespace {

/** Translate python dict to Config
 *
 *  {}      -> Config
 *    float   -> double
 *    str     -> string
 *    [{}]    -> vector<Config>  (recurse)
 *    ndarray -> vector<double>
 *    TODO: [0.0]   -> vector<double>
 */
static
void Dict2Config(Config& ret, PyObject *dict, unsigned depth=0)
{
    if(depth>3)
        throw std::runtime_error("too deep for Dict2Config");

    PyObject *key, *value;
    Py_ssize_t pos = 0;

    while(PyDict_Next(dict, &pos, &key, &value)) {
#if PY_MAJOR_VERSION >= 3
        PyRef<> ascii(PyUnicode_AsASCIIString(key));
        const char *kname = PyBytes_AsString(ascii.py());
#else
        const char *kname = PyString_AsString(key);
#endif

        if(!kname)
            throw std::invalid_argument("dict() keys must be string");

        PyTypeObject *valuetype = (PyTypeObject*)PyObject_Type(value);
        if(valuetype==&PyFloat_Type) { // scalar double
            double val = PyFloat_AsDouble(value);
            ret.set<double>(kname, val);

        } else if(valuetype==&PyInt_Type) { // scalar integer (treated as double)
            long val = PyInt_AsLong(value);
            ret.set<double>(kname, val);

        } else if(PyString_Check(value)) { // string
#if PY_MAJOR_VERSION >= 3
            PyRef<> kascii(PyUnicode_AsASCIIString(value));
            const char *val = PyBytes_AsString(kascii.py());
#else
            const char *val = PyString_AsString(value);
#endif

            ret.set<std::string>(kname, val);

        } else if(PyArray_Check(value)) { // array (ndarray)
            PyRef<> arr(PyArray_ContiguousFromAny(value, NPY_DOUBLE, 0, 2));
            double *buf = (double*)PyArray_DATA(arr.py());
            std::vector<double> temp(PyArray_SIZE(arr.py()));
            std::copy(buf, buf+temp.size(), temp.begin());

            ret.set<std::vector<double> >(kname, temp);

        } else if(PySequence_Check(value)) { // list of dict
            Py_ssize_t N = PySequence_Size(value);

            Config::vector_t output;
            output.reserve(N);

            for(Py_ssize_t i=0; i<N; i++) {
                PyObject *elem = PySequence_GetItem(value, i);
                assert(elem);

                if(!PyDict_Check(elem))
                    throw std::invalid_argument("lists must contain only dict()s");

                output.push_back(ret.new_scope());

                Dict2Config(output.back(), elem, depth+1); // inheirt parent scope
            }

            ret.set<Config::vector_t>(kname, output);

        } else {
            std::ostringstream msg;
            msg<<"Must be a dict, not "<<valuetype->tp_name;
            throw std::invalid_argument(msg.str());
        }
    }
}

static
PyObject* conf2dict(const Config *conf);

namespace {
struct confval : public boost::static_visitor<PyObject*>
{
    PyObject* operator()(double v) const
    {
        return PyFloat_FromDouble(v);
    }

    PyObject* operator()(const std::string& v) const
    {
        return PyString_FromString(v.c_str());
    }

    PyObject* operator()(const std::vector<double>& v) const
    {
        npy_intp dims[]  = {(npy_intp)v.size()};
        PyRef<> obj(PyArray_SimpleNew(1, dims, NPY_DOUBLE));
        std::copy(v.begin(), v.end(), (double*)PyArray_DATA(obj.py()));
        return obj.release();
    }

    PyObject* operator()(const Config::vector_t& v) const
    {
        PyRef<> L(PyList_New(v.size()));

        for(size_t i=0, N=v.size(); i<N; i++)
        {
            PyList_SetItem(L.py(), i, conf2dict(&v[i]));
        }
        return L.release();
    }
};
}

static
PyObject* conf2dict(const Config *conf)
{
    PyRef<> ret(PyDict_New());

    for(Config::const_iterator it=conf->begin(), end=conf->end();
        it!=end; ++it)
    {
        if(PyDict_SetItemString(ret.py(), it->first.c_str(),
                                boost::apply_visitor(confval(), it->second)
                                ))
            throw std::runtime_error("Failed to insert into dictionary from conf2dict");
    }

    return ret.release();
}

} // namespace

Config* dict2conf(PyObject *dict)
{
    if(!PyDict_Check(dict))
        throw std::invalid_argument("Not a dict");
    std::auto_ptr<Config> conf(new Config);
    Dict2Config(*conf, dict);
    return conf.release();
}

PyObject* PyGLPSPrint(PyObject *, PyObject *args)
{
    try {
        PyObject *dict;
        if(!PyArg_ParseTuple(args, "O!", &PyDict_Type, &dict))
            return NULL;

        Config conf;
        Dict2Config(conf, dict);
        std::ostringstream strm;
        GLPSPrint(strm, conf);
        return PyString_FromString(strm.str().c_str());
    }CATCH()
}

#ifndef PY_SSIZE_T_CLEAN
#error the following assumes ssize_t is used
#endif

PyObject* PyGLPSParse(PyObject *, PyObject *args)
{
    try{
        const char *buf;
        Py_ssize_t blen;
        if(!PyArg_ParseTuple(args, "s#", &buf, &blen))
            return NULL;

        GLPSParser parser;
        std::auto_ptr<Config> conf(parser.parse(buf, blen));
        return conf2dict(conf.get());

    }CATCH()
}
