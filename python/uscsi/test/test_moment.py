
import unittest
import numpy
from numpy import testing as NT
from numpy.testing import assert_array_almost_equal_nulp as assert_aequal

from .. import Machine

class testBasic(unittest.TestCase):
  def setUp(self):
    T = self.expect = numpy.asfarray([
      [1, 0, 1, 0, 1, 0],
      [0, 1, 0, 1, 0, 1],
      [1, 0, 1, 0, 1, 0],
      [0, 1, 0, 1, 0, 1],
      [1, 0, 1, 0, 1, 0],
      [0, 1, 0, 1, 0, 1],
    ])
    I = self.expect0 = numpy.asfarray(
      [1, 1, 0, 0, 0, 0]
    )
    self.M = Machine({
      'sim_type':'MomentMatrix',
      'elements':[
        {'name':'elem0', 'type':'source', 'initial':T, 'moment0':I},
        {'name':'elem1', 'type':'generic', 'transfer':numpy.identity(6)},
      ],
    })

  def test_generic(self):
    "Propogate a state matrix through a generic section"

    S = self.M.allocState({})

    self.M.propagate(S)

    assert_aequal(S.moment0, self.expect0)
    assert_aequal(S.state, self.expect)

  def test_reconfig(self):
    self.M.reconfigure(1, {'transfer':numpy.identity(6)*2.0})

    S = self.M.allocState({})

    self.M.propagate(S)

    assert_aequal(S.moment0, self.expect0*2)
    assert_aequal(S.state, self.expect*4)

    self.M.reconfigure(1, {'transfer':numpy.identity(6)*5.0})

    S = self.M.allocState({})

    self.M.propagate(S)

    assert_aequal(S.moment0, self.expect0*5)
    assert_aequal(S.state, self.expect*25)
