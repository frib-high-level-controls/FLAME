
import unittest
import numpy
from numpy import testing as NT
from numpy.testing import assert_array_almost_equal_nulp as assert_aequal

from .. import Machine

class testBasic(unittest.TestCase):
  def setUp(self):
    self.M = Machine({
      'sim_type':'Vector',
      'elements':[
        {'name':'elem0', 'type':'drift', 'L':1.0},
      ],
    })

  def test_print(self):
    self.assertEqual(str(self.M), """sim_type: Vector
#Elements: 1
Element 0: elem0 (drift)
Transfer: [6,6]((1,1,0,0,0,0),(0,1,0,0,0,0),(0,0,1,1,0,0),(0,0,0,1,0,0),(0,0,0,0,1,0),(0,0,0,0,0,1))
""")

  def test_drift(self):
    "Propogate a state vector through a drift section"

    S = self.M.allocState({})

    S.state[:] = [0, 0, 0, 0, 1, 1e-3]
    assert_aequal(S.state, [0, 0, 0, 0, 1.000, 1e-3])

    self.M.propagate(S)

    assert_aequal(S.state, [0, 0, 0, 0, 1.000, 1e-3])

  def test_gc(self):
    "See that State attributes have appropriate lifetime"

    import weakref, gc

    S = self.M.allocState({})
    state = S.state

    R = weakref.ref(S)

    del S
    gc.collect()

    # S should be kept alive by reference from state
    self.assertIsNot(R(), None)

    del state

    gc.collect()
    S = R()
    self.assertIs(R(), None)

  def test_err(self):
    "Try to propagate the something which is not a State"
    self.assertRaises(ValueError, self.M.propagate, None)

class TestMatrix(unittest.TestCase):
  def setUp(self):
    self.M = Machine({
      'sim_type':'TransferMatrix',
      'elements':[
        {'name':'elem0', 'type':'drift', 'L':1.0},
        {'name':'elem1', 'type':'drift', 'L':1.0},
      ],
    })

  def test_drift(self):
    """Propogate an identity matrix to find
    the cumulative transfer matrix of two drifts
    """

    S = self.M.allocState({})

    assert_aequal(S.state, [
      [1, 0, 0, 0, 0, 0],
      [0, 1, 0, 0, 0, 0],
      [0, 0, 1, 0, 0, 0],
      [0, 0, 0, 1, 0, 0],
      [0, 0, 0, 0, 1, 0],
      [0, 0, 0, 0, 0, 1],
    ])

    self.M.propagate(S)

    assert_aequal(S.state, [
      [1, 2, 0, 0, 0, 0],
      [0, 1, 0, 0, 0, 0],
      [0, 0, 1, 2, 0, 0],
      [0, 0, 0, 1, 0, 0],
      [0, 0, 0, 0, 1, 0],
      [0, 0, 0, 0, 0, 1],
    ])

class TestGlobal(unittest.TestCase):
    def test_parse(self):
        "Test global scope when parsing"
        M = Machine(b"""
sim_type = "Vector";
L = 2.0;
elem0: drift;
elem1: drift;
elem2: drift;
foo: LINE = (elem0, elem1, elem2);
""")

        S = M.allocState({})

        S.state[:] = [0, 0, 1, 1e-3, 0, 0]
        assert_aequal(S.state, [0, 0, 1.000, 1e-3, 0, 0])

        M.propagate(S)

        assert_aequal(S.state, [0, 0, 1.006, 1e-3, 0, 0])

class testGeneric(unittest.TestCase):
    def test_generic(self):
        T = numpy.asfarray([
          [1, 0, 1, 0, 1, 0],
          [0, 1, 0, 1, 0, 1],
          [1, 0, 1, 0, 1, 0],
          [0, 1, 0, 1, 0, 1],
          [1, 0, 1, 0, 1, 0],
          [0, 1, 0, 1, 0, 1],
        ])
        M = Machine({
          'sim_type':'TransferMatrix',
          'elements':[
            {'name':'elem0', 'type':'generic', 'transfer':T},
          ],
        })

        S = M.allocState({})

        M.propagate(S)

        assert_aequal(S.state, T)

class TestSource(unittest.TestCase):
    def test_matrix(self):
        T = numpy.asfarray([1, 0, 1, 0, 1, 0])
        M = Machine({
          'sim_type':'Vector',
          'elements':[
            {'name':'elem0', 'type':'source', 'initial':T},
          ],
        })

        S = M.allocState({})

        M.propagate(S)

        assert_aequal(S.state, T)

    def test_matrix(self):
        T = numpy.asfarray([
          [1, 0, 1, 0, 1, 0],
          [0, 1, 0, 1, 0, 1],
          [1, 0, 1, 0, 1, 0],
          [0, 1, 0, 1, 0, 1],
          [1, 0, 1, 0, 1, 0],
          [0, 1, 0, 1, 0, 1],
        ])
        M = Machine({
          'sim_type':'TransferMatrix',
          'elements':[
            {'name':'elem0', 'type':'source', 'initial':T},
          ],
        })

        S = M.allocState({})

        M.propagate(S)

        assert_aequal(S.state, T)
