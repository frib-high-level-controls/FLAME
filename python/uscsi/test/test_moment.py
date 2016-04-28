
from __future__ import print_function

import unittest
import numpy
from numpy import testing as NT
from numpy.testing import assert_array_almost_equal_nulp as assert_aequal

from .. import Machine
import sys

class testBasic(unittest.TestCase):
  def setUp(self):
    T = self.expect = numpy.asfarray([
      [1, 0, 1, 0, 1, 0, 0],
      [0, 1, 0, 1, 0, 1, 0],
      [1, 0, 1, 0, 1, 0, 0],
      [0, 1, 0, 1, 0, 1, 0],
      [1, 0, 1, 0, 1, 0, 0],
      [0, 1, 0, 1, 0, 1, 0],
      [0, 0, 0, 0, 0, 0, 0]
    ])
    I = self.expect0 = numpy.asfarray(
      [1, 1, 0, 0, 0, 0, 0]
    )
    self.M = Machine({
      'sim_type':'MomentMatrix',
      'elements':[
        {'name':'elem0', 'type':'source', 'initial':T, 'moment0':I},
        {'name':'elem1', 'type':'generic', 'transfer':numpy.identity(7)},
      ],
    })

  def test_generic(self):
    "Propogate a state matrix through a generic section"

    S = self.M.allocState({})

    self.M.propagate(S)

    assert_aequal(S.moment0, self.expect0)
    assert_aequal(S.state, self.expect)

  def test_reconfig(self):
    self.M.reconfigure(1, {'transfer':numpy.identity(7)*2.0})

    S = self.M.allocState({})

    self.M.propagate(S)

    assert_aequal(S.moment0, self.expect0*2)
    assert_aequal(S.state, self.expect*4)

    self.M.reconfigure(1, {'transfer':numpy.identity(7)*5.0})

    S = self.M.allocState({})

    self.M.propagate(S)

    assert_aequal(S.moment0, self.expect0*5)
    assert_aequal(S.state, self.expect*25)

class testMoment2Single(unittest.TestCase):
    def setUp(self):
        self.M = Machine(b'''
sim_type = "MomentMatrix2";
Frf = 80.5e6;
IonEs = 1.0;
IonEk = 1.0;
IM = [1,0,1,0,1,0,0,
      0,1,0,1,0,1,0,
      1,0,1,0,1,0,0,
      0,1,0,1,0,1,0,
      1,0,1,0,1,0,0,
      0,1,0,1,0,1,0,
      0,0,0,0,0,0,0
      ];
IV = [1, 1, 0, 0, 0, 0, 0];
TM = [1,0,0,0,0,0,0,
      0,1,0,0,0,0,0,
      0,0,1,0,0,0,0,
      0,0,0,1,0,0,0,
      0,0,0,0,1,0,0,
      0,0,0,0,0,1,0,
      0,0,0,0,0,0,1];
elem0 : source, initial = IM, moment0=IV;
elem1 : generic, transfer = TM, L=2.0;
foo : LINE = (elem0, elem1);
''')

    def test_config(self):
        C = self.M.conf()
        self.assertTrue("elements" in C)
        self.assertEqual(C['IonEk'], 1.0)
        assert_aequal(C['IV'], numpy.asfarray([1,1,0,0,0,0,0]))

    expect = numpy.asfarray([
        [1, 0, 1, 0, 1, 0, 0],
        [0, 1, 0, 1, -0.001193, 1, 0],
        [1, 0, 1, 0, 1, 0, 0],
        [0, 1, 0, 1, -0.001193, 1, 0],
        [1, -0.001193, 1.00000142, -0.001193, 1, -0.001193, 0],
        [0, 1, 0, 1, -0.001193, 1, 0],
        [0, 0, 0, 0, 0, 0, 0]
    ])

    def test_generic(self):
        "Propogate a state matrix through a generic section"
        C = self.M.conf()

        S = self.M.allocState({})
        self.M.propagate(S, 0, 1)

        ##S.moment0[:] = self.expect0
        self.assertEqual(S.pos, 0.0)
        self.assertEqual(S.real_Ekinetic, 0.0)
        self.assertEqual(S.real_gamma, 1.0)
        self.assertAlmostEqual(S.real_beta, 0.0)
        self.assertEqual(S.gamma, 2.0)
        self.assertAlmostEqual(S.beta, 0.8660254037844386)

        print("moment0",  S.moment0, C['IV'])
        assert_aequal(S.moment0, C['IV'])
        print("state", S.state, C['IM'])
        assert_aequal(S.state, C['IM'].reshape((7,7)), 1e10)

        self.M.propagate(S, 1, len(self.M))

        self.assertEqual(S.pos, 2.0)
        print("moment0",  S.moment0, C['IV'])
        # m56 is re-computed.
        assert_aequal(S.state, [
          [ 1.00000000e+00,  0.00000000e+00,  1.00000000e+00,  0.00000000e+00,  1.00000000e+00,  0.00000000e+00,  0.00000000e+00],
          [ 0.00000000e+00,  1.00000000e+00,  0.00000000e+00,  1.00000000e+00, -3.37431049e+06,  1.00000000e+00,  0.00000000e+00],
          [ 1.00000000e+00,  0.00000000e+00,  1.00000000e+00,  0.00000000e+00,  1.00000000e+00,  0.00000000e+00,  0.00000000e+00],
          [ 0.00000000e+00,  1.00000000e+00,  0.00000000e+00,  1.00000000e+00, -3.37431049e+06,  1.00000000e+00,  0.00000000e+00],
          [ 1.00000000e+00, -3.37431049e+06,  1.00000000e+00, -3.37431049e+06,  1.13859713e+13, -3.37431049e+06,  0.00000000e+00],
          [ 0.00000000e+00,  1.00000000e+00,  0.00000000e+00,  1.00000000e+00, -3.37431049e+06,  1.00000000e+00,  0.00000000e+00],
          [ 0.00000000e+00,  0.00000000e+00,  0.00000000e+00,  0.00000000e+00,  0.00000000e+00,  0.00000000e+00,  0.00000000e+00]
        ], 1e10)

    def test_source(self):
        "Initialize w/ all zeros, then propagate through source element to overwrite"
        C = self.M.conf()

        S = self.M.allocState({}, inherit=False)

        self.M.propagate(S, max=1) # propagate through source element

        self.assertEqual(S.pos, 0.0)
        self.assertEqual(S.ref_IonEk, 1.0)
        self.assertEqual(S.ref_IonEs, 1.0)
        self.assertEqual(S.ref_gamma, 2.0)
        self.assertAlmostEqual(S.ref_beta, 0.8660254037844386)

        print("moment0",  S.moment0, C['IV'])
        assert_aequal(S.moment0, C['IV'])
        print("state", S.state, C['IM'])
        assert_aequal(S.state, C['IM'].reshape((7,7)), 1e10)

class testMoment2Multi(unittest.TestCase):
    lattice = b'''
sim_type = "MomentMatrix2";
Frf = 80.5e6;
IonEs = 1.0;
IonEk = 1.0;
IM0 = [1,0,1,0,1,0,0,
      0,1,0,1,0,1,0,
      1,0,1,0,1,0,0,
      0,1,0,1,0,1,0,
      1,0,1,0,1,0,0,
      0,1,0,1,0,1,0,
      0,0,0,0,0,0,0
      ];
IM1 = [2,0,1,0,1,0,0,
      0,1,0,1,0,1,0,
      1,0,1,0,1,0,0,
      0,1,0,1,0,1,0,
      1,0,1,0,1,0,0,
      0,1,0,1,0,1,0,
      0,0,0,0,0,0,0
      ];
IV0 = [1, 1, 0, 0, 0, 0, 0];
IV1 = [2, 2, 0, 0, 0, 0, 0];
TM = [1,0,0,0,0,0,0,
      0,1,0,0,0,0,0,
      0,0,1,0,0,0,0,
      0,0,0,1,0,0,0,
      0,0,0,0,1,0,0,
      0,0,0,0,0,1,0,
      0,0,0,0,0,0,1];
IonChargeStates = [42, 43];
elem0 : source, vector_variable="IV", matrix_variable="IM";
elem1 : generic, transfer = TM, L=2.0;
foo : LINE = (elem0, elem1);
'''

    def test_source(self):
        """See that source element initializes correctly.

        Use cstate=1 to select
        S.IonZ = IonChargeStates[1]
        S.moment0 = IV1
        S.state = IM1
        """
        M = Machine(self.lattice, extra={"cstate":1})
        C = M.conf()

        S = M.allocState({}, inherit=False) # defaults

        M.propagate(S, max=1) # propagate through source element

        self.assertEqual(S.pos, 0.0)
        self.assertEqual(S.ref_IonEk, 1.0)
        self.assertEqual(S.ref_IonEs, 1.0)
        self.assertEqual(S.ref_gamma, 2.0)
        self.assertAlmostEqual(S.ref_beta, 0.8660254037844386)

        print("moment0",  S.moment0, C['IV1'])
        assert_aequal(S.moment0, C['IV1'])
        print("state", S.state, C['IM1'])
        assert_aequal(S.state, C['IM1'].reshape((7,7)), 1e10)
