#!/usr/bin/env python
"""
Emit a file suitible to include() in a CMakeLists.txt file
with information from a python interpreter

Compatible for python 2.6 -> 3.4
"""

from __future__ import print_function

import sys, os

if len(sys.argv)<2:
    out = sys.stdout
else:
    out = open(sys.argv[1], 'w')

from distutils.sysconfig import get_config_var, get_python_inc, get_python_lib

incdirs = [get_python_inc()]
libdirs = [get_config_var('LIBDIR')]

have_np='NO'
try:
    from numpy import get_include
    incdirs += [get_include()]
    have_np='YES'
except ImportError:
    pass

libdirs = [get_config_var('LIBDIR')]

# prepend introspected numpy directory so that it is checked before
# system python directory, which may contained a different version
# when virtualenv is used.  Debian helpfully symlinks the numpy headers
# as /usr/include/pythonX.Y/numpy :P
libdirs.reverse()
incdirs.reverse()

# location of extension modules relative to prefix (eg. "lib/python3/dist-packages")
moddir = get_python_lib()

print('set(Python_DEFINITIONS, "%s")'%get_config_var('BASECFLAGS'), file=out)

print('set(Python_VERSION "%s")'%get_config_var('VERSION'), file=out)
print('set(Python_VERSION_LD "%s")'%(get_config_var('LDVERSION') or get_config_var('VERSION')), file=out)
print('set(Python_INCLUDE_DIRS "%s")'%';'.join(incdirs), file=out)
print('set(Python_LIBRARY_DIRS "%s")'%';'.join(libdirs), file=out)
print('set(Python_MODULE_DIR "%s")'%moddir, file=out)
print('set(Python_NUMPY_FOUND %s)'%have_np, file=out)

print('set(Python_VERSION_MAJOR %s)'%sys.version_info[0], file=out)
print('set(Python_VERSION_MINOR %s)'%sys.version_info[1], file=out)
print('set(Python_VERSION_PATCH %s)'%sys.version_info[2], file=out)

print('set(Python_FOUND YES)', file=out)
