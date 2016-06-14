import sys
import math
import numpy

home_dir  = '/home/johan/git_repos/jmbgsddb/build/python'

sys.path.append(home_dir)

from flame import Machine
from flame import GLPSPrinter
from flame import _GLPSParse

class GLPSParser(object):
    def parse(self, s):
        return _GLPSParse(s)

__all__ = ['Machine',
    'GLPSPrinter',
    'GLPSParser',
]


PS_X = 0; PS_PX = 1; PS_Y = 2; PS_PY = 3; PS_S = 4; PS_PS = 5

MeVtoeV = 1e6


def sqr(a): return a*a;


P = GLPSParser()


def prt_lat_prms(M):
    print '\nIon Charge States = ', M.conf()['IonChargeStates']
    print 'IonEs [MeV]       = ', M.conf()['IonEs']/1e6
    print 'IonEk [MeV]       = ', M.conf()['IonEk']/1e6
    print '\nBaryCenter 1:\n', M.conf()['BaryCenter0']
    print '\nBaryCenter 2:\n', M.conf()['BaryCenter1']
    print '\nBeam Envelope 1:\n', M.conf()['S0']
    print '\nBeam Envelope 2:\n', M.conf()['S1']


def track_chg_state(IonZ_ref, IonZ, m1, m2):
    S = M.allocState({})
    M.propagate(S, 0, 1)

    S.ref_IonZ    = IonZ_ref
    S.real_IonZ   = IonZ

    S.moment0[:]  = m1
    S.state[:]    = numpy.split(m2, 7)

    S.real_gamma  = S.real_IonW/S.real_IonEs;
    S.real_beta   = math.sqrt(1e0-1e0/sqr(S.real_gamma));
    S.real_bg     = S.real_beta*S.real_gamma;

    S.real_phis   = S.moment0[PS_S];
    S.real_IonEk += S.moment0[PS_PS]*MeVtoeV;

    M.propagate(S, 1, len(M))
    print S


file_name = home_dir+'/flame/test/to_strl.lat'

with open(file_name, 'rb') as inf:
    M = Machine(inf.read())

    # print '\n', dir(M)
    prt_lat_prms(M)
    # prt_lat(M)

    Stripper_GetMat(conf, sim, ST)

    track_chg_state(M.conf()['IonChargeStates'][0],
          M.conf()['IonChargeStates'][0],
          M.conf()['BaryCenter1'], M.conf()['S1'])

