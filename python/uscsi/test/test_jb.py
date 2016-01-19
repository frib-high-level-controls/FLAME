
import sys
import unittest
import math
import numpy
from numpy import testing as NT
#from numpy.testing import assert_array_almost_equal_nulp as assert_aequal
from numpy.testing import assert_array_almost_equal as assert_aequal

from .. import Machine

from .._internal import (GLPSPrinter as dictshow, _GLPSParse)

import os
datadir = os.path.dirname(__file__)

def print_state(S):
  n = 6
  for i in range(n):
    sys.stdout.write('[')
    for j in range(n):
      sys.stdout.write('%9.6f' % (S.state[i, j]))
      if j != n-1: sys.stdout.write(', ')
    if i != n-1:
      sys.stdout.write('],\n')
    else:
      sys.stdout.write(']\n')

class GLPSParser(object):
    def parse(self, s):
        return _GLPSParse(s)

class testBasic(unittest.TestCase):
  def setUp(self):
    # Called before each test.

    # T = self.expect = numpy.asfarray([
    #   [1, 0, 0, 0, 0, 0],
    #   [0, 1, 0, 0, 0, 0],
    #   [0, 0, 1, 0, 0, 0],
    #   [0, 0, 0, 1, 0, 0],
    #   [0, 0, 0, 0, 1, 0],
    #   [0, 0, 0, 0, 0, 1]
    # ])

   # self.M = Machine({
   #   'sim_type':'MomentMatrix',
   #   'elements':[
   #     {'name':'elem0', 'type':'source', 'initial':T},
   #     {'name':'elem2', 'type':'dipole', 'L':1.234},
   #     {'name':'elem1', 'type':'generic', 'transfer':numpy.identity(6)},
   #   ]
   # })

   return

  def test_generic(self):
    "Propagate a state matrix through a generic section"

    T = self.expect = numpy.asfarray([
      [1, 0, 0, 0, 0, 0],
      [0, 1, 0, 0, 0, 0],
      [0, 0, 1, 0, 0, 0],
      [0, 0, 0, 1, 0, 0],
      [0, 0, 0, 0, 1, 0],
      [0, 0, 0, 0, 0, 1]
    ])

    self.M = Machine({
      'sim_type':'MomentMatrix',
      'elements':[
        {'name':'elem0', 'type':'source', 'initial':T},
        {'name':'elem1', 'type':'generic', 'transfer':numpy.identity(6)}
      ]
    })

    S = self.M.allocState({})

#    S.state[:] = [
#      [1, 0, 0, 0, 0, 0],
#      [0, 1, 0, 0, 0, 0],
#      [0, 0, 1, 0, 0, 0],
#      [0, 0, 0, 1, 0, 0],
#      [0, 0, 0, 0, 1, 0],
#      [0, 0, 0, 0, 0, 1]
#    ]

    self.M.propagate(S)

    assert_aequal(S.state, [
      [1, 0, 0, 0, 0, 0],
      [0, 1, 0, 0, 0, 0],
      [0, 0, 1, 0, 0, 0],
      [0, 0, 0, 1, 0, 0],
      [0, 0, 0, 0, 1, 0],
      [0, 0, 0, 0, 0, 1]
    ])

  def test_drift(self):
    "Propagate a state matrix through a drift"

    T = self.expect = numpy.asfarray([
      [0.1, 0,   0,   0,   0,   0  ],
      [0,   0.2, 0,   0,   0,   0  ],
      [0,   0,   0.3, 0,   0,   0  ],
      [0,   0,   0,   0.4, 0,   0  ],
      [0,   0,   0,   0,   0.5, 0  ],
      [0,   0,   0,   0,   0,   0.6]
    ])

    self.M = Machine({
      'sim_type':'MomentMatrix',
      'elements':[
        {'name':'elem0', 'type':'source', 'initial':T},
        {'name':'elem1', 'type':'drift', 'L':1.234},
      ]
    })

    S = self.M.allocState({})

    self.M.propagate(S)

    assert_aequal(S.state, [
      [ 0.404551,  0.246800,  0.000000,  0.000000,  0.000000,  0.000000],
      [ 0.246800,  0.200000,  0.000000,  0.000000,  0.000000,  0.000000],
      [ 0.000000,  0.000000,  0.909102,  0.493600,  0.000000,  0.000000],
      [ 0.000000,  0.000000,  0.493600,  0.400000,  0.000000,  0.000000],
      [ 0.000000,  0.000000,  0.000000,  0.000000,  0.761378,  0.000000],
      [ 0.000000,  0.000000,  0.000000,  0.000000,  0.000000,  0.600000]
    ])

  def test_sbend1(self):
    "Propagate a state matrix through a focusing sector bend"

    T = self.expect = numpy.asfarray([
      [0.1, 0,   0,   0,   0,   0  ],
      [0,   0.2, 0,   0,   0,   0  ],
      [0,   0,   0.3, 0,   0,   0  ],
      [0,   0,   0,   0.4, 0,   0  ],
      [0,   0,   0,   0,   0.5, 0  ],
      [0,   0,   0,   0,   0,   0.6]
    ])

    self.M = Machine({
      'sim_type':'MomentMatrix',
      'elements':[
        {'name':'elem0', 'type':'source', 'initial':T},
        {'name':'elem1', 'type':'sbend', 'L':2.0, 'phi':math.radians(25.0),
         'K':1.0},
      ]
    })

    S = self.M.allocState({})

    self.M.propagate(S)

    assert_aequal(S.state, [
      [ 0.171806, -0.037912,  0.000000,  0.000000,  0.000000,  0.000000],
      [-0.037912,  0.124777,  0.000000,  0.000000,  0.000000,  0.000000],
      [ 0.000000,  0.000000,  9.507881,  9.551471,  0.000000,  0.000000],
      [ 0.000000,  0.000000,  9.551471,  9.607881,  0.000000,  0.000000],
      [ 0.000000,  0.000000,  0.000000,  0.000000,  2.000000,  0.000000],
      [ 0.000000,  0.000000,  0.000000,  0.000000,  0.000000,  0.600000]
   ], decimal=6)

  def test_sbend2(self):
    "Propagate a state matrix through a defocusing sector bend"

    T = self.expect = numpy.asfarray([
      [0.1, 0,   0,   0,   0,   0  ],
      [0,   0.2, 0,   0,   0,   0  ],
      [0,   0,   0.3, 0,   0,   0  ],
      [0,   0,   0,   0.4, 0,   0  ],
      [0,   0,   0,   0,   0.5, 0  ],
      [0,   0,   0,   0,   0,   0.6]
    ])

    self.M = Machine({
      'sim_type':'MomentMatrix',
      'elements':[
        {'name':'elem0', 'type':'source', 'initial':T},
        {'name':'elem1', 'type':'sbend', 'L':2.0, 'phi':math.radians(25.0),
         'K':-1.0},
      ]
    })

    S = self.M.allocState({})

    self.M.propagate(S)

    assert_aequal(S.state, [
      [ 3.789181,  3.748527,  0.000000,  0.000000,  0.000000,  0.000000],
      [ 3.748527,  3.713589,  0.000000,  0.000000,  0.000000,  0.000000],
      [ 0.000000,  0.000000,  0.382682, -0.037840,  0.000000,  0.000000],
      [ 0.000000,  0.000000, -0.037840,  0.317318,  0.000000,  0.000000],
      [ 0.000000,  0.000000,  0.000000,  0.000000,  2.000000,  0.000000],
      [ 0.000000,  0.000000,  0.000000,  0.000000,  0.000000,  0.600000]
    ], decimal=6)

  def test_quad1(self):
    "Propagate a state matrix through a focusing quadrupole"

    T = self.expect = numpy.asfarray([
      [0.1, 0,   0,   0,   0,   0  ],
      [0,   0.2, 0,   0,   0,   0  ],
      [0,   0,   0.3, 0,   0,   0  ],
      [0,   0,   0,   0.4, 0,   0  ],
      [0,   0,   0,   0,   0.5, 0  ],
      [0,   0,   0,   0,   0,   0.6]
    ])

    self.M = Machine({
      'sim_type':'MomentMatrix',
      'elements':[
        {'name':'elem0', 'type':'source', 'initial':T},
        {'name':'elem1', 'type':'quadrupole', 'L':2.0, 'K':1.1},
      ]
    })

    S = self.M.allocState({})

    self.M.propagate(S)

    assert_aequal(S.state, [
      [ 0.161135, -0.037295,  0.000000,  0.000000,  0.000000,  0.000000],
      [-0.037295,  0.132752,  0.000000,  0.000000,  0.000000,  0.000000],
      [ 0.000000,  0.000000, 10.981961, 11.546105,  0.000000,  0.000000],
      [ 0.000000,  0.000000, 11.546105, 12.150157,  0.000000,  0.000000],
      [ 0.000000,  0.000000,  0.000000,  0.000000,  2.000000,  0.000000],
      [ 0.000000,  0.000000,  0.000000,  0.000000,  0.000000,  0.600000]
    ], decimal=6)

  def test_quad2(self):
    "Propagate a state matrix through a defocusing quadrupole"

    T = self.expect = numpy.asfarray([
      [0.1, 0,   0,   0,   0,   0  ],
      [0,   0.2, 0,   0,   0,   0  ],
      [0,   0,   0.3, 0,   0,   0  ],
      [0,   0,   0,   0.4, 0,   0  ],
      [0,   0,   0,   0,   0.5, 0  ],
      [0,   0,   0,   0,   0,   0.6]
    ])

    self.M = Machine({
      'sim_type':'MomentMatrix',
      'elements':[
        {'name':'elem0', 'type':'source', 'initial':T},
        {'name':'elem1', 'type':'quadrupole', 'L':2.0, 'K':-1.1},
      ]
    })

    S = self.M.allocState({})

    self.M.propagate(S)

    assert_aequal(S.state, [
      [ 4.636175,  4.903140,  0.000000,  0.000000,  0.000000,  0.000000],
      [ 4.903140,  5.189793,  0.000000,  0.000000,  0.000000,  0.000000],
      [ 0.000000,  0.000000,  0.347549, -0.029007,  0.000000,  0.000000],
      [ 0.000000,  0.000000, -0.029007,  0.347696,  0.000000,  0.000000],
      [ 0.000000,  0.000000,  0.000000,  0.000000,  2.000000,  0.000000],
      [ 0.000000,  0.000000,  0.000000,  0.000000,  0.000000,  0.600000]
    ], decimal=6)

  def test_solenoid(self):
    "Propagate a state matrix through a defocusing quadrupole"

    T = self.expect = numpy.asfarray([
      [0.1, 0,   0,   0,   0,   0  ],
      [0,   0.2, 0,   0,   0,   0  ],
      [0,   0,   0.3, 0,   0,   0  ],
      [0,   0,   0,   0.4, 0,   0  ],
      [0,   0,   0,   0,   0.5, 0  ],
      [0,   0,   0,   0,   0,   0.6]
    ])

    self.M = Machine({
      'sim_type':'MomentMatrix',
      'elements':[
        {'name':'elem0', 'type':'source', 'initial':T},
        {'name':'elem1', 'type':'solenoid', 'L':1.123, 'K':0.0},
      ]
    })

    S = self.M.allocState({})

    self.M.propagate(S)

    # assert_aequal(S.state, [
    #   [ 4.636175,  4.903140,  0.000000,  0.000000,  0.000000,  0.000000],
    #   [ 4.903140,  5.189793,  0.000000,  0.000000,  0.000000,  0.000000],
    #   [ 0.000000,  0.000000,  0.347549, -0.029007,  0.000000,  0.000000],
    #   [ 0.000000,  0.000000, -0.029007,  0.347696,  0.000000,  0.000000],
    #   [ 0.000000,  0.000000,  0.000000,  0.000000,  2.000000,  0.000000],
    #   [ 0.000000,  0.000000,  0.000000,  0.000000,  0.000000,  0.600000]
    # ], decimal=6)

  def test_lat(self):
    # Propagate a state matrix through a lattice defined by a lattice file.
    P = GLPSParser()

    with open(os.path.join(datadir, 'moment_jb.lat'), 'rb') as inf:
      self.M = Machine(inf.read())
#    print self.M

    S = self.M.allocState({})

    self.M.propagate(S)
#    print_state(S)
    assert_aequal(S.state, [
      [ 5.989936,  0.321261,  0.000000,  0.000000,  0.000000,  0.000000],
      [ 0.321261,  0.020569,  0.000000,  0.000000,  0.000000,  0.000000],
      [ 0.000000,  0.000000, 15.331388,  3.337662,  0.000000,  0.000000],
      [ 0.000000,  0.000000,  3.337662,  0.734440,  0.000000,  0.000000],
      [ 0.000000,  0.000000,  0.000000,  0.000000, 13.430708,  0.000000],
      [ 0.000000,  0.000000,  0.000000,  0.000000,  0.000000,  0.600000]
    ], decimal=6)

  def test_parse(self):
    "Test global scope when parsing"
    # Propagate a state matrix through a lattice defined by a lattice file.
    P = GLPSParser()

    inf = open(os.path.join(datadir, 'moment_jb_2.lat'), 'r')
    self.M = Machine(inf.read())
    print self.M

    S = self.M.allocState({})
    print S

 #    assert_aequal(S.state, [0, 0, 1.000, 1e-3, 0, 0])

    self.M.propagate(S)

#    assert_aequal(S.state, [0, 0, 1.006, 1e-3, 0, 0])

  def test_LS1(self):
    # Propagate a state matrix through a lattice defined by a lattice file.

    # Global parameters:

#    AU        = 931.49432                 # MeV/u.
#    qom_U_238 = 33.0/(238.0*AU)           # Charge over mass ratio for U-238.

#    E_mass = AU                           # Ion type.
#    E_kin  = 0.9149311118819696           # Kinetic energy [MeV/u].
#    E_tot  = sqrt(sqr(E_kin)+sqr(E_Mass)) # Total energy [MeV/u].
#    E_tot  = E_kin + E_Mass

#    gamma = numpy.array([E_tot/E_mass, 0.0, 0.0]);

    P = GLPSParser()

    with open(os.path.join(datadir, 'latticeout_IMP_withPV.lat'), 'rb') as inf:
      self.M = Machine(inf.read())
#    print self.M

    S = self.M.allocState({})

    self.M.propagate(S)
#    print_state(S)
    # assert_aequal(S.state, [
    #   [ 5.367983,  0.300692,  0.000000,  0.000000,  0.000000,  0.000000],
    #   [ 0.300692,  0.020569,  0.000000,  0.000000,  0.000000,  0.000000],
    #   [ 0.000000,  0.000000,  9.390504,  2.603222,  0.000000,  0.000000],
    #   [ 0.000000,  0.000000,  2.603222,  0.734440,  0.000000,  0.000000],
    #   [ 0.000000,  0.000000,  0.000000,  0.000000, 13.430708,  0.000000],
    #   [ 0.000000,  0.000000,  0.000000,  0.000000,  0.000000,  0.600000]
    # ], decimal=6)
