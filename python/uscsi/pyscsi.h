#ifndef PYSCSI_H
#define PYSCSI_H

#include <boost/python/dict.hpp>

struct Config;
struct StateBase;

Config* dict2conf(const boost::python::dict& O);

PyObject* wrapstate(StateBase*);
StateBase* unwrapstate(PyObject*);

void registerModConfig(void);
void registerModMachine(void);
void registerModState(void);

#endif // PYSCSI_H
