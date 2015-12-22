
from ._internal import Machine
from ._internal import GLPSPrinter, _GLPSParse

class GLPSParser(object):
    def parse(self, s):
        return _GLPSParse(s)

__all__ = ['Machine',
    'GLPSPrinter',
    'GLPSParser',
]
