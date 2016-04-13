
from ._internal import Machine
from ._internal import GLPSPrinter, _GLPSParse

class GLPSParser(object):
    def parse(self, s, path=None):
        return _GLPSParse(s, path=path)

__all__ = ['Machine',
    'GLPSPrinter',
    'GLPSParser',
]
