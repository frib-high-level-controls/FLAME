
import unittest
from numpy import testing as NT
from numpy.testing import assert_array_almost_equal_nulp as assert_aequal

from .. import Machine

class testBasic(unittest.TestCase):
  def setUp(self):
    self.M = Machine({
      'sim-type':'Vector',
      'elements':[
        {'name':'elem0', 'type':'drift', 'length':1.0},
      ],
    })

  def test_drift(self):
    "Propogate a state vector through a drift section"

    S = self.M.allocState({})

    S.state[:] = [0, 0, 0, 0, 1, 1e-3]
    assert_aequal(S.state, [0, 0, 0, 0, 1.000, 1e-3])

    self.M.propogate(S)

    assert_aequal(S.state, [0, 0, 0, 0, 1.001, 1e-3])

  def test_gc(self):
    "See that State attributes have appropriate lifetime"

    import weakref, gc
    S = self.M.allocState({})
    R = weakref.ref(S)
    state = S.state
    del S
    gc.collect()
    self.assertIsNot(R(), None)

  def test_err(self):
    "Try to propogate the something which is not a State"
    self.assertRaises(TypeError, self.M.propogate, None)

class TestMatrix(unittest.TestCase):
  def setUp(self):
    self.M = Machine({
      'sim-type':'TransferMatrix',
      'elements':[
        {'name':'elem0', 'type':'drift', 'length':1.0},
        {'name':'elem1', 'type':'drift', 'length':1.0},
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

    self.M.propogate(S)

    assert_aequal(S.state, [
      [1, 0, 0, 0, 0, 0],
      [0, 1, 0, 0, 0, 0],
      [0, 0, 1, 0, 0, 0],
      [0, 0, 0, 1, 0, 0],
      [0, 0, 0, 0, 1, 2],
      [0, 0, 0, 0, 0, 1],
    ])
