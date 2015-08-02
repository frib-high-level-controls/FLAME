
#include <list>
#include <sstream>

//#include <boost/python.hpp>
#include <boost/python/def.hpp>
#include <boost/python/class.hpp>
#include <boost/python/extract.hpp>
#include <boost/python/dict.hpp>
#include <boost/python/list.hpp>
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
bool XnL(std::list<boost::any>& L, bp::object& O)
{
    bp::extract<T> E(O);
    if(E.check()) {
        L.push_back(E());
        return true;
    }
    return false;
}

static
void Dict2Config(Config& ret, const bp::dict& O, unsigned depth=0)
{
    if(depth>3)
        throw std::runtime_error("too deep for Dict2Config");

    PyObject* it(O.iteritems().ptr()), *item;
    while( !!(item=PyIter_Next(it)) )
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

        if(XnS<int>(ret, name, value))
            continue;
        else if(XnS<double>(ret, name, value))
            continue;
        else if(XnS<std::string>(ret, name, value))
            continue;
        else if(PyArray_Check(value.ptr()))
        {
            if(PyArray_NDIM(value.ptr())!=2)
                throw std::runtime_error("ndarray ndim!=2");
            bp::object arr(bp::handle<>(PyArray_ContiguousFromAny(value.ptr(), NPY_DOUBLE, 0, 2)));
            double *buf = (double*)PyArray_DATA(arr.ptr());
            std::vector<double> temp(PyArray_SIZE(arr.ptr()));
            std::copy(buf, buf+temp.size(), temp.begin());

            ret.setAny(name, temp);
            continue;
        }

        {
            bp::extract<bp::dict> E(value);
            if(E.check()) {
                Config recurse;
                Dict2Config(recurse, E(), depth+1);
                ret.setAny(name, recurse);
                continue;
            }
        }
        {
            bp::extract<bp::list> E(value);
            if(E.check()) {
                std::list<boost::any> output;

                for(size_t i=0, len = bp::len(value); i<len; i++)
                {
                    bp::object ent(value[i]);

                    if(XnL<int>(output, ent))
                        continue;
                    else if(XnL<double>(output, ent))
                        continue;
                    else if(XnL<std::string>(output, ent))
                        continue;
                    else {
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
                ret.setAny(name, output);
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

void registerModConfig(void)
{
    using namespace boost::python;

    def("dictshow", showConfig);
    def("dict2conf", dict2conf, return_value_policy<bp::manage_new_object>());

    class_<Config, boost::noncopyable>("Config", no_init);
}
