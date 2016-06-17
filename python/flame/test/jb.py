import sys
import math
import numpy

sys.path.append('/home/johan/git_repos/jmbgsddb/build/python')

from uscsi import Machine
from uscsi import GLPSPrinter
from uscsi import _GLPSParse

class GLPSParser(object):
    def parse(self, s):
        return _GLPSParse(s)

__all__ = ['Machine',
    'GLPSPrinter',
    'GLPSParser',
]

# Must be executed from: ../build/src

PS_X = 0; PS_PX = 1; PS_Y = 2; PS_PY = 3; PS_S = 4; PS_PS = 5

MeVtoeV = 1e6


def sqr(a): return a*a; 


home_dir = '/home/johan/git_repos/jmbgsddb/build/python/uscsi/test/'

P = GLPSParser()

#file_name = home_dir+'LS1.lat'
#file_name = home_dir+'to_strl.lat'
file_name = home_dir+'to_strl_err.lat'

with open(file_name, 'rb') as inf:
    M = Machine(inf)

    S = M.allocState({})
    M.propagate(S, 0)

    print S.moment0
    print S.moment0_rms
    print S.state
