
#include <list>
#include <sstream>

//#include <boost/python.hpp>
#include <boost/python/def.hpp>
#include <boost/python/class.hpp>
#include <boost/python/extract.hpp>
#include <boost/python/dict.hpp>
#include <boost/python/list.hpp>
#include <boost/python/str.hpp>
#include <boost/python/return_value_policy.hpp>
#include <boost/python/manage_new_object.hpp>

#define NO_IMPORT_ARRAY
#define PY_ARRAY_UNIQUE_SYMBOL USCSI_PyArray_API
#include <numpy/ndarrayobject.h>

#include "scsi/base.h"

#include "pyscsi.h"

namespace bp = boost::python;

static
std::string
showConfig(const Config* c)
{
    std::ostringstream strm;
    strm << *c;
    return strm.str();
}

template<typename T>
static inline
bool XnS(Config& C, const std::string& n, bp::object& O)
{
    bp::extract<T> E(O);
    if(E.check()) {
        C.set<T>(n, E());
        return true;
    }
    return false;
}

template<typename T>
static inline
bool XnL(Config::vector_t& L, bp::object& O)
{
    bp::extract<T> E(O);
    if(E.check()) {
        L.push_back(E());
        return true;
    }
    return false;
}

/** Translate python dict to Config
 *
 *  {}      -> Config
 *  float   -> double
 *  str     -> string
 *  [{}]    -> vector<Config>
 *  ndarray -> vector<double>
 *  TODO: [0.0]   -> vector<double>
 */
static
void Dict2Config(Config& ret, const bp::dict& O, unsigned depth=0)
{
    if(depth>3)
        throw std::runtime_error("too deep for Dict2Config");

    bp::object it;
#if PY_MAJOR_VERSION >= 3
    {
        // it seems that in python 3 land it is necessary to create an
        // intermediate "set-like" object to iterate a dict as tuples...
        bp::object temp(O.items());
        it = bp::object(bp::handle<>(PyObject_GetIter(temp.ptr())));
    }
#else
    it = bp::object(O.iteritems());
#endif
    PyObject *item;
    while( !!(item=PyIter_Next(it.ptr())) )
    {
        std::string name;
        bp::object value;
        {
            bp::object K(bp::handle<>(bp::borrowed(PyTuple_GetItem(item, 0))));
            value = bp::object(bp::handle<>(bp::borrowed(PyTuple_GetItem(item, 1))));

            Py_DECREF(item);
            bp::extract<std::string> X(K);
            if(!X.check())
                throw std::runtime_error("keys must be strings");

            name = X();
        }

        if(XnS<double>(ret, name, value))
            continue;
        else if(XnS<std::string>(ret, name, value))
            continue;
        else if(PyArray_Check(value.ptr()))
        {
            bp::object arr(bp::handle<>(PyArray_ContiguousFromAny(value.ptr(), NPY_DOUBLE, 0, 2)));
            double *buf = (double*)PyArray_DATA(arr.ptr());
            std::vector<double> temp(PyArray_SIZE(arr.ptr()));
            std::copy(buf, buf+temp.size(), temp.begin());

            ret.set<std::vector<double> >(name, temp);
            continue;
        }

        {
            bp::extract<bp::list> E(value);
            if(E.check()) {
                Config::vector_t output;
                ssize_t L = bp::len(value);
                output.reserve(L>=0 ? L : -L);

                for(size_t i=0, len = L; i<len; i++)
                {
                    bp::object ent(value[i]);

                    {
                        bp::extract<bp::dict> X(ent);
                        if(X.check()) {
                            Config recurse;
                            Dict2Config(recurse, X(), depth+1);
                            output.push_back(recurse);
                        } else {
                            throw std::invalid_argument("list entry");
                        }
                    }
                }
                ret.set<Config::vector_t>(name, output);
                continue;
            }
        }

        throw key_error(name);
    }
}

Config* dict2conf(const bp::dict& O)
{
    std::auto_ptr<Config> conf(new Config);
    Dict2Config(*conf, O);
    return conf.release();
}

static
bp::dict conf2dict(const Config *conf);

namespace {
struct confval : public boost::static_visitor<bp::object>
{
    bp::object operator()(double v) const
    {
        return bp::object(bp::handle<>(PyFloat_FromDouble(v)));
    }

    bp::object operator()(const std::string& v) const
    {
        return bp::str(v.c_str());
    }

    bp::object operator()(const std::vector<double>& v) const
    {
        npy_intp dims[]  = {(npy_intp)v.size()};
        PyArrayObject *obj;
        obj = (PyArrayObject *)PyArray_SimpleNew(1, dims, NPY_DOUBLE);
        bp::object ret(bp::handle<>((PyObject*)obj)); // throws if obj==NULL
        std::copy(v.begin(), v.end(), (double*)PyArray_DATA(obj));
        return ret;
    }

    bp::object operator()(const Config::vector_t& v) const
    {
        bp::list L;
        for(size_t i=0, N=v.size(); i<N; i++)
        {
            L.append(conf2dict(&v[i]));
        }
        return L;
    }
};
}

static
bp::dict conf2dict(const Config *conf)
{
    bp::dict ret;

    for(Config::const_iterator it=conf->begin(), end=conf->end();
        it!=end; ++it)
    {
        bp::str key(it->first.c_str());
        ret[key] = boost::apply_visitor(confval(), it->second);
    }

    return ret;
}

namespace {
struct PyGLPSParser : public boost::noncopyable
{
    GLPSParser parser;

    Config* parse(const bp::str& s)
    {
        bp::extract<const char *> X(s);
        if(!X.check())
            throw std::invalid_argument("string required");
        return parser.parse(X());
    }
};
}

bp::str PyGLPSPrint(const Config* c)
{
    std::ostringstream strm;
    GLPSPrint(strm, *c);

    return bp::str(strm.str().c_str());
}

void registerModConfig(void)
{
    using namespace boost::python;

    def("dictshow", showConfig);
    def("dict2conf", dict2conf, return_value_policy<bp::manage_new_object>());
    def("conf2dict", conf2dict);
    def("GLPSPrinter",  &PyGLPSPrint);

    class_<Config, boost::noncopyable>("Config", no_init);

    class_<PyGLPSParser, boost::noncopyable>("GLPSParser")
            .def("parse", &PyGLPSParser::parse, return_value_policy<bp::manage_new_object>())
            ;
}
