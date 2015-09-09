
from ._internal import Machine
from ._internal import GLPSPrinter, GLPSParser
from ._internal import dict2conf, conf2dict

Machine.__doc__ = """Simulation execution engine.

A Machine object is responsible for carrying out calculations
based on a Config provided when it was constructed.

See the allocState and propagate methods.
"""

__all__ = ['Machine',
    'GLPSPrinter',
    'GLPSParser',
    'dict2conf',
    'conf2dict',
]
