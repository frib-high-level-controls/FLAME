
#include <list>
#include <sstream>

//#include <boost/python.hpp>
#include <boost/python/module.hpp>
#include <boost/python/exception_translator.hpp>

#define PY_ARRAY_UNIQUE_SYMBOL USCSI_PyArray_API
#include <numpy/ndarrayobject.h>

#include "scsi/base.h"
#include "pyscsi.h"

namespace bp = boost::python;

static
void translate_key(const key_error& e)
{
    PyErr_SetString(PyExc_KeyError, e.what());
}

static
void translate_any(const boost::bad_any_cast& e)
{
    PyErr_SetString(PyExc_ValueError, e.what());
}

BOOST_PYTHON_MODULE(_internal)
{
    using namespace boost::python;

    if (_import_array() < 0)
        throw std::runtime_error("Failed to import numpy");

    register_exception_translator<key_error>(translate_key);
    register_exception_translator<boost::bad_any_cast>(translate_any);

    // add python stuff to this module
    registerModConfig();
    registerModMachine();
    registerModState();

    // add States and Elements
    registerLinear();
    registerMoment();
}
