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

file_name = home_dir+'latticeout_IMP_withPV_consolidate2.lat'

with open(file_name, 'rb') as inf:
    M = Machine(inf.read())

    S = M.allocState({})
    M.propagate(S, 0, 1)
    S.moment0[:] = [
        -0.0007886,   1.08371e-05,
         0.01337343,  6.678534e-06,
        -0.0001847729, 0.000309995, 1.0
    ];

    S.real_gamma     = S.real_IonW/S.real_IonEs;
    S.real_beta      = math.sqrt(1e0-1e0/sqr(S.real_gamma));
    S.real_bg        = S.real_beta*S.real_gamma;

    S.real_phis      = S.moment0[PS_S];
    S.real_Ekinetic += S.moment0[PS_PS]*MeVtoeV;

    print S

    M.propagate(S, 1, len(M))

    print S
