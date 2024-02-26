
from collections import OrderedDict

from ._internal import (Machine as MachineBase,
                        GLPSPrinter, _GLPSParse,
                        _pyapi_version, _capi_version,
                        FLAME_ERROR, FLAME_WARN,
                        FLAME_INFO, FLAME_DEBUG,
                        setLogLevel, getLoggerName)

def _list2odict(L):
    'Recursively turn list of tuples into OrderedDict'
    for i in range(len(L)):
        K,V = L[i]
        if isinstance(V, list):
            L[i] = (K,list(map(OrderedDict, V)))
    return OrderedDict(L)

class GLPSParser(object):
    """GLPS parser context
    """
    def parse(self, *args, **kws):
        """parse(file_or_buf, path=None, extra=None)
        parse(file_or_buf, path="/dir/", extra={'VAR':'value'})

        Parse the provided buffer or file-like object.

        'path' is used to expand relative paths found while parsing.
        If not 'path' is None then either PWD or the .name of
        the file-like object is used.

        'extra' may be used to provide additional variable definitions when parsing.

        Returns an OrderedDict.
        """
        return _GLPSParse(*args, **kws)

class Machine(MachineBase):
    def conf(self, *args, **kws):
        return _list2odict(super(Machine, self).conf(*args, **kws))

# by default pass all but DEBUG to python logger.
# May set to FLAME_WARN for performance
setLogLevel(FLAME_WARN)

__all__ = ['Machine',
    'GLPSPrinter',
    'GLPSParser',
]

__version__ = '1.9.0'
