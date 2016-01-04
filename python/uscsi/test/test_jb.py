
import unittest
import numpy
from numpy import testing as NT
#from numpy.testing import assert_array_almost_equal_nulp as assert_aequal
from numpy.testing import assert_array_almost_equal as assert_aequal

from .. import Machine

class testBasic(unittest.TestCase):
  def setUp(self):
    T = self.expect = numpy.asfarray([
      [1, 0, 0, 0, 0, 0],
      [0, 1, 0, 0, 0, 0],
      [0, 0, 1, 0, 0, 0],
      [0, 0, 0, 1, 0, 0],
      [0, 0, 0, 0, 1, 0],
      [0, 0, 0, 0, 0, 1],
    ])
    self.M = Machine({
      'sim_type':'MomentMatrix',
      'elements':[
        #        {'name':'elem0', 'type':'source', 'initial':T},
        #        {'name':'elem2', 'type':'dipole', 'length':1.234},
        #        {'name':'elem1', 'type':'generic', 'transfer':numpy.identity(6)},
      ]
    })
#    print self.M

  def test_generic(self):
    "Propagate a state matrix through a generic section"

    self.M = Machine({
      'sim_type':'MomentMatrix',
      'elements':[
        {'name':'elem1', 'type':'generic', 'transfer':numpy.identity(6)}
      ]
    })

    S = self.M.allocState({})

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
    "Propagate a state vector through a drift section"

    self.M = Machine({
      'sim_type':'MomentMatrix',
      'elements':[
        {'name':'elem1', 'type':'drift', 'length':1.234},
      ]
    })

    S = self.M.allocState({})

    S.state[:] = [
      [0.1, 0,   0,   0,   0,   0  ],
      [0,   0.2, 0,   0,   0,   0  ],
      [0,   0,   0.3, 0,   0,   0  ],
      [0,   0,   0,   0.4, 0,   0  ],
      [0,   0,   0,   0,   0.5, 0  ],
      [0,   0,   0,   0,   0,   0.6]
    ]

    self.M.propagate(S)

    assert_aequal(S.state, [
      [0.4045512, 0.2468, 0,         0,      0,        0  ],
      [0.2468,    0.2,    0,         0,      0,        0  ],
      [0,         0,      0.9091024, 0.4936, 0,        0  ],
      [0,         0,      0.4936,    0.4,    0,        0  ],
      [0,         0,      0,         0,      0.761378, 0  ],
      [0,         0,      0,         0,      0,        0.6]
    ])

  def test_dipole(self):
    "Propagate a state vector through a dipole section"

    # Vertical is 0 or 1.
    self.M = Machine({
      'sim_type':'MomentMatrix',
      'elements':[
        {'name':'elem1', 'type':'dipole', 'vertical':0, 'radius':1.0,
         'angle':25.0*numpy.pi/180.0},
      ]
    })

    S = self.M.allocState({})

    S.state[:] = [
      [0.1, 0,   0,   0,   0,   0  ],
      [0,   0.2, 0,   0,   0,   0  ],
      [0,   0,   0.3, 0,   0,   0  ],
      [0,   0,   0,   0.4, 0,   0  ],
      [0,   0,   0,   0,   0.5, 0  ],
      [0,   0,   0,   0,   0,   0.6]
    ]

    self.M.propagate(S)

    assert_aequal(S.state, [
      [0.11786062, 0.03830222, 0,          0,          0,          0  ],
      [0.03830222, 0.18213938, 0,          0,          0,          0  ],
      [0,          0,          0.37615435, 0.17453293, 0,          0  ],
      [0,          0,          0.17453293, 0.4,        0,          0  ],
      [0,          0,          0,          0,          0.5,        0  ],
      [0,          0,          0,          0,          0,          0.6]
    ], decimal=6)
