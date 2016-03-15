
from __future__ import print_function

import unittest
import numpy
from numpy import testing as NT
from numpy.testing import assert_array_almost_equal_nulp as assert_aequal

from .. import Machine

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

class testMoment2(unittest.TestCase):
    def setUp(self):
      self.expect0 = numpy.asfarray(
        [1, 1, 0, 0, 0, 0, 0]
      )
      self.M = Machine(b'''
sim_type = "MomentMatrix2";
Frf = 80.5e6;
IonEs = 1.0;
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

        S = self.M.allocState({})
        self.assertEqual(S.pos, 0.0)
        self.assertEqual(S.Ekinetic, 0.0)
        self.assertEqual(S.gamma, 1.0)
        self.assertAlmostEqual(S.beta, numpy.sqrt(2.0))

        self.M.propagate(S)

        self.assertEqual(S.pos, 2.0)
        assert_aequal(S.moment0, self.expect0)
        assert_aequal(S.state, self.expect, 1e10)
