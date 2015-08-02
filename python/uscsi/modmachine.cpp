
#include <list>
#include <sstream>

//#include <boost/python.hpp>
#include <boost/python/class.hpp>
#include <boost/python/def.hpp>

#include "scsi/base.h"
#include "pyscsi.h"

namespace bp = boost::python;

struct PyMachine : public boost::noncopyable
{
    std::auto_ptr<Machine> machine;

    PyMachine(const bp::dict& c)
    {
        std::auto_ptr<Config> C(dict2conf(c));
        machine.reset(new Machine(*C));
    }

    std::string tostring() const
    {
        std::ostringstream strm;
        strm << *machine;
        return strm.str();
    }

    PyObject* allocState(const bp::dict& d)
    {
        std::auto_ptr<Config> C(dict2conf(d));
        std::auto_ptr<StateBase> state(machine->allocState(*C));
        PyObject *ret = wrapstate(state.get());
        state.release();
        return ret;
    }

    void propagate(PyObject* state, size_t start, size_t max)
    {
        machine->propagate(unwrapstate(state), start, max);
    }
};
void registerModMachine(void)
{
    using namespace boost::python;

    class_<PyMachine, boost::noncopyable>("Machine", init<const bp::dict&>())
            .def("__str__", &PyMachine::tostring)
            .def("allocState", &PyMachine::allocState)
            .def("propagate", &PyMachine::propagate,
                 (arg("state"), arg("start")=0, arg("max")=(size_t)-1))
            ;

}
