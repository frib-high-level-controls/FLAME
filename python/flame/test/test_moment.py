
from __future__ import print_function

from math import sqrt
from collections import OrderedDict

import unittest
import numpy
from numpy import testing as NT
from numpy.testing import assert_array_almost_equal_nulp as assert_aequal

from .. import Machine, GLPSParser
import sys

class testMomentSingle(unittest.TestCase):
    def setUp(self):
        self.M = Machine(b'''
sim_type = "MomentMatrix";
Frf = 80.5e6;
IonEs = 930e6;
IonEk = 500e3;
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
foo : LINE = (elem0);
''')

    def test_config(self):
        C = self.M.conf()
        self.assertTrue("elements" in C)
        self.assertEqual(C['IonEk'], 500e3)
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

    def test_source(self):
        "Initialize w/ all zeros, then propagate through source element to overwrite"
        C = self.M.conf()

        S = self.M.allocState({}, inherit=False)

        self.M.propagate(S, max=1) # propagate through source element

        self.assertEqual(S.pos, 0.0)
        self.assertEqual(S.real_IonEs, C['IonEs'])
        self.assertEqual(S.real_IonW, C['IonEk']+C['IonEs'])
        self.assertEqual(S.real_IonEk, C['IonEk'])
        self.assertEqual(S.real_gamma, (C['IonEk']+C['IonEs'])/C['IonEs'])
        self.assertAlmostEqual(S.real_beta, sqrt(1-1/(S.real_gamma**2)))

        print("moment0",  S.moment0_env, C['IV'])
        assert_aequal(S.moment0_env, C['IV'])
        print("state", S.moment1_env, C['IM'])
        assert_aequal(S.moment1_env, C['IM'].reshape((7,7)), 1e10)

    def test_modify_state(self):
        S = self.M.allocState({}, inherit=False)

        # assign scalar
        S.pos = 42.2
        self.assertEqual(S.pos, 42.2)

        # assign 1d
        S.moment0_env = numpy.asfarray([1, 2, 3, 4, 5, 6, 7])
        assert_aequal(S.moment0_env, numpy.asfarray([1, 2, 3, 4, 5, 6, 7]))

        # assign 2d
        S.moment1_env = self.expect*2.0
        assert_aequal(S.moment1_env, self.expect*2.0)

        # assign 3d
        X = self.expect.reshape((7,7,1))
        S.moment1 = X*2.0
        assert_aequal(S.moment1, X*2.0)

class testMomentMulti(unittest.TestCase):
    lattice = b'''
sim_type = "MomentMatrix";
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
NCharge = [1000, 1010];
elem0 : source, vector_variable="IV", matrix_variable="IM";
foo : LINE = (elem0);
'''

    def test_through_dict(self):
        "Check that list is translated to Config w/ scope"
        P = GLPSParser()
        L = P.parse(self.lattice)
        print("E", L)
        M = Machine(L)

        # check that 'Frf' is defined in file scope, but not (explicitly) in element scope
        L = OrderedDict(L)
        self.assertIn('Frf', L)
        E = OrderedDict(L['elements'][0])
        self.assertNotIn('Frf', E)

        # however, 'Frf' is still found implicitly at element scope
        self.assertEqual(L['Frf'], M.conf()['Frf'])

    def test_source_single(self):
        """See that source element initializes correctly for a single charge state

        Use cstate=1 to select
        S.IonZ = IonChargeStates[1]
        S.moment0_env = IV1
        S.moment1_env = IM1
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

        assert_aequal(S.IonQ, C['NCharge'][1:])

        print("moment0",  S.moment0_env, C['IV1'])
        assert_aequal(S.moment0_env, C['IV1'])
        print("state", S.moment1_env, C['IM1'])
        assert_aequal(S.moment1_env, C['IM1'].reshape((7,7)), 1e10)

    def test_source_multi(self):
        """See that source element initializes correctly for many (two) charge states

        Use cstate=1 to select
        S.IonZ = IonChargeStates[1]
        S.moment0_env = IV1
        S.moment1_env = IM1
        """
        M = Machine(self.lattice)
        C = M.conf()

        S = M.allocState({}, inherit=False) # defaults

        M.propagate(S, max=1) # propagate through source element

        self.assertEqual(S.pos, 0.0)
        self.assertEqual(S.ref_IonEk, 1.0)
        self.assertEqual(S.ref_IonEs, 1.0)
        self.assertEqual(S.ref_gamma, 2.0)
        self.assertAlmostEqual(S.ref_beta, 0.8660254037844386)

        # check access to array attributes
        assert_aequal(S.IonQ, C['NCharge'])
        assert_aequal(S.moment0[...,0], C['IV0'])
        assert_aequal(S.moment0[...,1], C['IV1'])
        assert_aequal(S.moment0.shape, numpy.asarray([7, 2]))

        IM0 = C['IM0'].reshape((7,7))
        IM1 = C['IM1'].reshape((7,7))
        assert_aequal(S.moment1[...,0], IM0, 1e10)
        assert_aequal(S.moment1[...,1], IM1, 1e10)
        assert_aequal(S.moment1.shape, numpy.asarray([7, 7, 2]))

        # check consistency of *_env and *_rms stats

        # moment0_env is average of moment0 over charge state weighted by charge in each state
        W = (S.IonQ/S.IonQ.sum()).reshape((1,2)).repeat(7,axis=0)
        M0env = (S.moment0*W).sum(axis=1)

        print("moment0_env",  S.moment0_env, M0env)
        assert_aequal(S.moment0_env, M0env)

        # moment1_env is ... more complex
        IM  = numpy.zeros(IM0.shape)

        Qs, M0, M0env, M1 = S.IonQ, S.moment0, S.moment0_env, S.moment1
        # using a loop because outer() doesn't understand axis=
        for i in range(len(Qs)):
            Q = Qs[i]
            m0diff = M0[:,i]-M0env

            IM[:7,:7] += Q*(M1[...,i]+numpy.outer(m0diff, m0diff))
        IM /= Qs.sum()

        print("moment1_env", S.moment1_env, IM)
        assert_aequal(S.moment1_env, IM, 1e10)

        # moment0_rms is derived from moment1_env
        M0rms = numpy.sqrt(numpy.diagonal(S.moment1_env))
        print("moment0_rms", S.moment0_rms, M0rms)
        assert_aequal(S.moment0_rms, M0rms)
