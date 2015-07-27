#!/usr/bin/env python
from __future__ import print_function
from numpy.distutils.misc_util import get_numpy_include_dirs
print(';'.join(get_numpy_include_dirs()))
