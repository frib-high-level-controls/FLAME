
#include <list>
#include <sstream>

//#include <boost/python.hpp>
#include <boost/python/module.hpp>
#include <boost/python/exception_translator.hpp>

#include "scsi/base.h"
#include "pyscsi.h"

namespace bp = boost::python;

void translate_key(const key_error& e)
{
    PyErr_SetString(PyExc_KeyError, e.what());
}

BOOST_PYTHON_MODULE(_internal)
{
    using namespace boost::python;

    register_exception_translator<key_error>(translate_key);

    // add python stuff to this module
    registerModConfig();
    registerModMachine();
    registerModState();

    // add States and Elements
    registerLinear();
}
