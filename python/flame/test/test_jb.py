
from __future__ import print_function

import sys
import unittest
import math
import numpy
from numpy import asfarray
from numpy import testing as NT
from numpy.testing import assert_array_almost_equal as assert_aequal

from .. import Machine

from .. import GLPSParser

import os
datadir = os.path.dirname(__file__)

def print_state(S):
  n = 7
  sys.stdout.write('\n')
  for i in range(n):
    sys.stdout.write('[')
    for j in range(n):
      sys.stdout.write('%17.10e' % (S.moment1_env[i, j]))
      if j != n-1: sys.stdout.write(', ')
    if i != n-1:
      sys.stdout.write('],\n')
    else:
      sys.stdout.write(']\n')

class MomentTest(object):
    'Helper for testing moment2 sim'

    def assertConsistent(self, S, msg=None):
        'Check internal consistencies of a State'
        def checkConsist(self, P='real'):
            self.assertEqual(getattr(S, P+'_IonW')    , getattr(S, P+'_IonEs')+getattr(S, P+'_IonEk'), msg)
            self.assertEqual(getattr(S, P+'_gamma')   , getattr(S, P+'_IonW')/getattr(S, P+'_IonEs'), msg)
            self.assertEqual(getattr(S, P+'_beta')    , math.sqrt(1e0-1e0/sqr(getattr(S, P+'_gamma'))), msg)
            self.assertEqual(getattr(S, P+'_bg')      , getattr(S, P+'_beta')*getattr(S, P+'_gamma'), msg)
        checkConsist('ref')
        checkConsist('real')

    def assertStateEqual(self, expect, actual, msg=None, decimal=10):
        'Assert that two moment2 State instances are equal'
        for k,v in expect.items():
            #self.assertIn(k, actual, "%s %s missing"%(msg or "",k)) #TODO: make State iterable
            A = getattr(actual, k, None)
            if A is None:
                self.assertTrue(False, "%s %s missing"%(msg or "",k))
            if isinstance(v, numpy.ndarray):
                assert_aequal(v, getattr(actual, k), decimal=decimal, err_msg="%s %s doesn't match"%(msg or "",k))
            else:
                self.assertAlmostEqual(v, getattr(actual, k), places=decimal, msg="%s %s doesn't match"%(msg or "",k))


    def checkPropagate(self, elem, instate, outstate, max=1):
        '''Pass given input state through the named element

        Propagate the same input twice in succession to verify that element(s)
        give the same output.
        The repeat with 'transfer' caching disabled
        '''

        S1 = self.M.allocState(instate, inherit=False)
        S2 = self.M.allocState(instate, inherit=False)

        self.M.propagate(state=S1, start=elem, max=max)
        self.M.propagate(state=S2, start=elem, max=max)

#        print_state(S1)

        self.assertStateEqual(outstate, S1, 'first pass')
        self.assertStateEqual(outstate, S2, 'second pass') # fails if some Element has a caching problem...

        S1 = self.ICM.allocState(instate, inherit=False)
        S2 = self.ICM.allocState(instate, inherit=False)

        self.ICM.propagate(state=S1, start=elem, max=max)
        self.ICM.propagate(state=S2, start=elem, max=max)

        self.assertStateEqual(outstate, S1, 'third pass')
        self.assertStateEqual(outstate, S2, 'fourth pass')

class TestToStrl(unittest.TestCase, MomentTest):
    """Strategy is to test the state after the first instance of each element type.
    """
    lattice = 'to_strl.lat'

    def setUp(self):
        with open(os.path.join(datadir, self.lattice), 'rb') as F:
            self.M = Machine(F)
            F.seek(0)
            self.ICM = Machine(F, extra={'skipcache':1.0})

    def test_source(self):
        self.checkPropagate(0, {}, {
            'next_elem':1,
            'ref_IonZ':0.13865546218487396,
            'ref_IonEs':931494320.0,
            'ref_IonEk':500000.0,
            'real_IonZ':0.13865546218487396,
            'real_IonEs':931494320.0,
            'real_IonEk':500309.99500000,
            'moment0_env':asfarray([3.3444511955e-03,1.2841341839e-05,8.3106826155e-03,-5.0277881002e-07,1.1632697213e-02,1.2028520756e-03,1.0000000000e+00]),
            'moment1_env':asfarray([
                [ 2.7784895410e+00,  2.4291087928e-04, -2.8819146986e-03, -1.3098890047e-05,  5.3594209519e-04, -1.7328482801e-05,  0.0000000000e+00],
                [ 2.4291087928e-04,  4.0675140732e-06,  9.1913901431e-06, -2.2427664520e-08,  3.5190231444e-07, -1.1570581361e-10,  0.0000000000e+00],
                [-2.8819146986e-03,  9.1913901431e-06,  2.5433968050e+00,  2.3922914294e-04, -8.7109112521e-04, -2.3018796436e-06,  0.0000000000e+00],
                [-1.3098890047e-05, -2.2427664520e-08,  2.3922914294e-04,  4.6600708351e-06, -8.7385463546e-07,  1.7247530508e-08,  0.0000000000e+00],
                [ 5.3594209519e-04,  3.5190231444e-07, -8.7109112521e-04, -8.7385463546e-07,  1.0267518940e-03,  1.0056757657e-05,  0.0000000000e+00],
                [-1.7328482801e-05, -1.1570581361e-10, -2.3018796436e-06,  1.7247530508e-08,  1.0056757657e-05,  2.6314343481e-06,  0.0000000000e+00],
                [ 0.0000000000e+00,  0.0000000000e+00,  0.0000000000e+00,  0.0000000000e+00,  0.0000000000e+00,  0.0000000000e+00,  0.0000000000e+00]
            ]),
        })

    def test_drift1(self):
        # drift_1
        self.checkPropagate(0, {}, {
            'next_elem':2,
            'ref_IonZ':1.3865546218e-01,
            'ref_IonEs':9.3149432000e+08,
            'ref_IonEk':5.0000000000e+05,
            'ref_phis':3.7089624016e+00,
            'real_IonZ':1.3865546218e-01,
            'real_IonEs':9.3149432000e+08,
            'real_IonEk':5.0030999500e+05,
            'real_phis':3.7076293282e+00,
            'moment0_env':asfarray([4.2690278079e-03,  1.2841341839e-05,  8.2744825412e-03, -5.0277881002e-07,  7.1994323704e-03,  1.2028520756e-03,  1.0000000000e+00]),
            'moment1_env':asfarray([
                [ 2.8345547006e+00,  5.3577189256e-04, -3.2795197045e-03, -1.4713681892e-05,  6.2542645808e-04, -1.7336813620e-05,  0.0000000000e+00],
                [ 5.3577189256e-04,  4.0675140732e-06,  7.5765982977e-06, -2.2427664520e-08,  3.5228427569e-07, -1.1570581361e-10,  0.0000000000e+00],
                [-3.2795197045e-03,  7.5765982977e-06,  2.6020036088e+00,  5.7475424307e-04, -9.3022235764e-04, -1.0600574470e-06,  0.0000000000e+00],
                [-1.4713681892e-05, -2.2427664520e-08,  5.7475424307e-04,  4.6600708351e-06, -9.3792991960e-07,  1.7247530508e-08,  0.0000000000e+00],
                [ 6.2542645808e-04,  3.5228427569e-07, -9.3022235764e-04, -9.3792991960e-07,  9.8879108726e-04,  3.4841141434e-07,  0.0000000000e+00],
                [-1.7336813620e-05, -1.1570581361e-10, -1.0600574470e-06,  1.7247530508e-08,  3.4841141434e-07,  2.6314343481e-06,  0.0000000000e+00],
                [ 0.0000000000e+00,  0.0000000000e+00,  0.0000000000e+00,  0.0000000000e+00,  0.0000000000e+00,  0.0000000000e+00,  0.0000000000e+00]
            ]),
        }, max=2)

    def test_rfcav_41(self):
        # ls1_ca01_cav1_d1127  cavtype = "0.041QWR"
        self.checkPropagate(0, {}, {
            'next_elem':4,
            'ref_IonZ':1.3865546218e-01,
            'ref_IonEs':9.3149432000e+08,
            'ref_IonEk':553528.7992770672,
            'ref_phis':22.74717198685346,
            'real_IonZ':1.3865546218e-01,
            'real_IonEs':9.3149432000e+08,
            'real_IonEk':553638.1760196686,
            'real_phis':22.741352002135365,
            'moment0_env':asfarray([9.8156783003e-03,  2.0041768075e-05, -5.2898622977e-03, -1.8638562712e-04, -1.7223152071e-02,  1.7187912211e-03,  1.0000000000e+00]),
            'moment1_env':asfarray([
                [ 4.5737889068e+00,  6.0669648277e-03, -1.0792259945e-02, -3.1759168651e-05,  9.1002516921e-04,  3.0358922803e-05,  0.0000000000e+00],
                [ 6.0669648277e-03,  1.0279529426e-05, -1.0736164288e-05, -4.6551002950e-08,  1.0542984436e-06,  4.2234819235e-08,  0.0000000000e+00],
                [-1.0792259945e-02, -1.0736164288e-05,  4.4609155107e+00,  6.1807101898e-03, -1.1255393909e-03, -6.3650656953e-05,  0.0000000000e+00],
                [-3.1759168651e-05, -4.6551002950e-08,  6.1807101898e-03,  1.0961176841e-05, -1.8057731754e-06, -9.3605511068e-08,  0.0000000000e+00],
                [ 9.1002516921e-04,  1.0542984436e-06, -1.1255393909e-03, -1.8057731754e-06,  1.2760917190e-03, -4.6806127721e-08,  0.0000000000e+00],
                [ 3.0358922803e-05,  4.2234819235e-08, -6.3650656953e-05, -9.3605511068e-08, -4.6806127721e-08,  4.2230866649e-06,  0.0000000000e+00],
                [ 0.0000000000e+00,  0.0000000000e+00,  0.0000000000e+00,  0.0000000000e+00,  0.0000000000e+00,  0.0000000000e+00,  5.5079869820e-18]
            ]),
        }, max=4)

    def test_sol(self):
        # ls1_ca01_sol1_d1131_1
        self.checkPropagate(0, {}, {
            'next_elem':9,
            'ref_IonZ':1.3865546218e-01,
            'ref_IonEs':9.3149432000e+08,
            'ref_IonEk':553528.7992770672,
            'ref_phis':37.252760222292466,
            'real_IonZ':1.3865546218e-01,
            'real_IonEs':9.3149432000e+08,
            'real_IonEk':553638.1760196686,
            'real_phis':37.245508581428005,
            'moment0_env':asfarray([-5.8659453275e-03, -4.0712274747e-05, -5.9179321956e-02, -1.1742389615e-04, -3.9532448827e-02,  1.7187912211e-03,  1.0000000000e+00]),
            'moment1_env':asfarray([
                [ 8.0532835806e+00, -1.1996773374e-03, -1.2076556961e-02,  7.8260918481e-05,  4.4722636729e-04,  7.7455375507e-06,  0.0000000000e+00],
                [-1.1996773374e-03,  1.4386139131e-06,  9.9286883346e-05, -3.4739396769e-09, -2.5932535513e-07, -3.8185010165e-09,  0.0000000000e+00],
                [-1.2076556961e-02,  9.9286883346e-05,  8.0927099123e+00, -9.7902740601e-04, -6.4020082998e-04, -9.4677007214e-05,  0.0000000000e+00],
                [ 7.8260918481e-05, -3.4739396769e-09, -9.7902740601e-04,  1.4431734181e-06, -4.0226350812e-08,  1.6715286132e-08,  0.0000000000e+00],
                [ 4.4722636729e-04, -2.5932535513e-07, -6.4020082998e-04, -4.0226350812e-08,  1.9877692144e-03, -5.4891187320e-05,  0.0000000000e+00],
                [ 7.7455375507e-06, -3.8185010165e-09, -9.4677007214e-05,  1.6715286132e-08, -5.4891187320e-05,  4.2230866649e-06,  0.0000000000e+00],
                [ 0.0000000000e+00,  0.0000000000e+00,  0.0000000000e+00,  0.0000000000e+00,  0.0000000000e+00,  0.0000000000e+00,  1.4697180793e-39]
            ]),
        }, max=9)

    def test_orbtrim(self):
        # ls1_ca01_dch_d1131_2
        self.checkPropagate(0, {}, {
            'next_elem':10,
            'ref_IonZ':1.3865546218e-01,
            'ref_IonEs':9.3149432000e+08,
            'ref_IonEk':553528.7992770672,
            'ref_phis':37.252760222292466,
            'real_IonZ':1.3865546218e-01,
            'real_IonEs':9.3149432000e+08,
            'real_IonEk':553638.1760196686,
            'real_phis':37.245508581428005,
            'moment0_env':asfarray([-5.8659453275e-03, -4.0712274747e-05, -5.9179321956e-02, -1.1742389615e-04, -3.9532448827e-02,  1.7187912211e-03,  1.0000000000e+00]),
            'moment1_env':asfarray([
                [ 8.0532835806e+00, -1.1996773374e-03, -1.2076556961e-02,  7.8260918481e-05,  4.4722636729e-04,  7.7455375507e-06,  0.0000000000e+00],
                [-1.1996773374e-03,  1.4386139131e-06,  9.9286883346e-05, -3.4739396769e-09, -2.5932535513e-07, -3.8185010165e-09,  0.0000000000e+00],
                [-1.2076556961e-02,  9.9286883346e-05,  8.0927099123e+00, -9.7902740601e-04, -6.4020082998e-04, -9.4677007214e-05,  0.0000000000e+00],
                [ 7.8260918481e-05, -3.4739396769e-09, -9.7902740601e-04,  1.4431734181e-06, -4.0226350812e-08,  1.6715286132e-08,  0.0000000000e+00],
                [ 4.4722636729e-04, -2.5932535513e-07, -6.4020082998e-04, -4.0226350812e-08,  1.9877692144e-03, -5.4891187320e-05,  0.0000000000e+00],
                [ 7.7455375507e-06, -3.8185010165e-09, -9.4677007214e-05,  1.6715286132e-08, -5.4891187320e-05,  4.2230866649e-06,  0.0000000000e+00],
                [ 0.0000000000e+00,  0.0000000000e+00,  0.0000000000e+00,  0.0000000000e+00,  0.0000000000e+00,  0.0000000000e+00,  1.4697180793e-39]
            ]),
        }, max=10)

    def test_rfcav_85(self):
        # ls1_cb01_cav1_d1229   cavtype = "0.085QWR"
        self.checkPropagate(0, {}, {
            'next_elem':111,
            'ref_IonZ':1.3865546218e-01,
            'ref_IonEs':9.3149432000e+08,
            'ref_IonEk':1556687.5901972055,
            'ref_phis':408.52804841338207,
            'real_IonZ':1.3865546218e-01,
            'real_IonEs':9.3149432000e+08,
            'real_IonEk':1556657.5690774918,
            'real_phis':408.52403688727543,
            'moment0_env':asfarray([-9.3631759606e-02, -1.2341993322e-04, -2.9374835471e-01, -1.1199135936e-03, -3.1818121415e-02,  6.0555868280e-04,  1.0000000000e+00]),
            'moment1_env':asfarray([
                [ 3.5741635753e+00,  3.9046681137e-03,  2.6319553841e-01,  1.8850927362e-04, -3.3171353330e-04, -2.5463826946e-05,  0.0000000000e+00],
                [ 3.9046681137e-03,  5.3877451153e-06,  2.4814486533e-04,  1.4838700828e-07, -7.1082274148e-07, -1.1436375795e-08,  0.0000000000e+00],
                [ 2.6319553841e-01,  2.4814486533e-04,  3.6460516525e+00,  3.2992822533e-03,  2.4128455791e-03, -6.4746360830e-05,  0.0000000000e+00],
                [ 1.8850927362e-04,  1.4838700828e-07,  3.2992822533e-03,  4.0807404978e-06,  3.2058212902e-06, -3.8424368637e-08,  0.0000000000e+00],
                [-3.3171353330e-04, -7.1082274148e-07,  2.4128455791e-03,  3.2058212902e-06,  1.2103473040e-03,  1.1651738943e-05,  0.0000000000e+00],
                [-2.5463826946e-05, -1.1436375795e-08, -6.4746360830e-05, -3.8424368637e-08,  1.1651738943e-05,  5.4283844066e-06,  0.0000000000e+00],
                [ 0.0000000000e+00,  0.0000000000e+00,  0.0000000000e+00,  0.0000000000e+00,  0.0000000000e+00,  0.0000000000e+00,  0.0000000000e+00]
            ]),
        }, max=111)

    def test_quad(self):
        # ls1_bts_qh_d1942
        self.checkPropagate(0, {}, {
            'next_elem':698,
            'ref_IonZ':1.3865546218e-01,
            'ref_IonEs':9.3149432000e+08,
            'ref_IonEk':17090945.594993234,
            'ref_phis':1340.7114319734278,
            'real_IonZ':1.3865546218e-01,
            'real_IonEs':9.3149432000e+08,
            'real_IonEk':17091872.10503781,
            'real_phis':1340.7128064371002,
            'moment0_env':asfarray([2.2103713766e+00,  3.3787339937e-03,  5.6689616131e-01, -1.9481100199e-04, -1.7814120019e-02, -5.4391793335e-03,  1.0000000000e+00]),
            'moment1_env':asfarray([
                [ 2.9112671746e+00,  4.3634419121e-03,  2.2706243347e-01, -4.8249852320e-04, -1.0221104296e-02, -3.5342038702e-03,  0.0000000000e+00],
                [ 4.3634419121e-03,  6.6698566950e-06,  3.4073808752e-04, -7.2120764900e-07, -1.4794010504e-05, -5.0989015094e-06,  0.0000000000e+00],
                [ 2.2706243347e-01,  3.4073808752e-04,  1.3314803329e+00, -1.8941670140e-03, -2.9206161602e-03, -9.8817258354e-04,  0.0000000000e+00],
                [-4.8249852320e-04, -7.2120764900e-07, -1.8941670140e-03,  3.0434962764e-06,  7.1940692081e-06,  2.4835020532e-06,  0.0000000000e+00],
                [-1.0221104296e-02, -1.4794010504e-05, -2.9206161602e-03,  7.1940692081e-06,  4.4492782545e-04,  8.1421698426e-05,  0.0000000000e+00],
                [-3.5342038702e-03, -5.0989015094e-06, -9.8817258354e-04,  2.4835020532e-06,  8.1421698426e-05,  8.1344632398e-05,  0.0000000000e+00],
                [ 0.0000000000e+00,  0.0000000000e+00,  0.0000000000e+00,  0.0000000000e+00,  0.0000000000e+00,  0.0000000000e+00,  0.0000000000e+00]
            ]),
        }, max=698)

    def test_sbend(self):
        # fs1_css_dh_d2163_1
        self.checkPropagate(0, {}, {
            'next_elem':831,
            'ref_IonZ':1.3865546218e-01,
            'ref_IonEs':9.3149432000e+08,
            'ref_IonEk':17090412.218117952,
            'ref_phis':1536.1552870706407,
            'real_IonZ':1.3865546218e-01,
            'real_IonEs':9.3149432000e+08,
            'real_IonEk':17089939.45941341,
            'real_phis':1536.153193681811,
            'moment0_env':asfarray([-1.6829029034e+00,  1.5337535377e-04, -1.6541522987e+00, -7.6642847191e-04,  9.9069891853e-03, -1.8048333267e-03,  1.0000000000e+00]),
            'moment1_env':asfarray([
                [ 1.7289456693e+00,  2.6671946343e-04,  1.5225241352e-01, -5.5867646567e-05, -1.5792405815e-03,  3.4344208786e-05,  0.0000000000e+00],
                [ 2.6671946343e-04,  8.8232084038e-07, -5.4049940765e-04, -5.3701103993e-08,  8.6683922816e-06, -9.8867507881e-07,  0.0000000000e+00],
                [ 1.5225241352e-01, -5.4049940765e-04,  3.1222535733e+00,  2.9140374831e-04, -6.8708500271e-03,  7.3814984344e-04,  0.0000000000e+00],
                [-5.5867646567e-05, -5.3701103993e-08,  2.9140374831e-04,  1.7019622036e-07, -6.9809042619e-08,  9.5400072773e-08,  0.0000000000e+00],
                [-1.5792405815e-03,  8.6683922816e-06, -6.8708500271e-03, -6.9809042619e-08,  6.3726044042e-04,  2.3522916789e-04,  0.0000000000e+00],
                [ 3.4344208786e-05, -9.8867507881e-07,  7.3814984344e-04,  9.5400072773e-08,  2.3522916789e-04,  1.3367141782e-04,  0.0000000000e+00],
                [ 0.0000000000e+00,  0.0000000000e+00,  0.0000000000e+00,  0.0000000000e+00,  0.0000000000e+00,  0.0000000000e+00,  0.0000000000e+00]
            ]),
        }, max=831)

    def test_stripper(self):
        # fs1_strl_strip_d2237_164
        pass

    def test_final(self):
        self.checkPropagate(0, {}, {
            'next_elem':len(self.M),
            'ref_IonZ':0.3277310924369748,
            'ref_IonEs':931494320.0,
            'ref_IonEk':16816951.191958785,
            'ref_phis':2023.6512871901991,
            'real_IonZ':0.31932773109243695,
            'real_IonEs':931494320.0,
            'real_IonEk':16829746.76773572,
            'real_phis':2023.6417753279932,
            'moment0_env':asfarray([1.0923868437e+00,3.5728723827e-04,1.4420500562e+00,-1.1959772575e-04,-1.7378583261e-02,-2.2999692554e-03,1.0000000000e+00]),
            'moment1_env':asfarray([
                [ 4.3986359287e+00, -1.0616512564e-04,  1.4411487761e-02, -3.0296029747e-04, -1.8869234689e-02,  2.0393230824e-03,  0.0000000000e+00],
                [-1.0616512564e-04,  4.9570504271e-07,  3.0901853376e-04,  3.9626728990e-08,  8.2462663320e-07, -9.8476930606e-07,  0.0000000000e+00],
                [ 1.4411487761e-02,  3.0901853376e-04,  6.5727958190e+00,  2.7540054815e-03,  7.7069472866e-03,  5.2455475777e-03,  0.0000000000e+00],
                [-3.0296029747e-04,  3.9626728990e-08,  2.7540054815e-03,  1.8784036222e-06,  9.3079770177e-06,  4.3760834645e-06,  0.0000000000e+00],
                [-1.8869234689e-02,  8.2462663320e-07,  7.7069472866e-03,  9.3079770177e-06,  5.0800488894e-04,  5.1394909319e-04,  0.0000000000e+00],
                [ 2.0393230824e-03, -9.8476930606e-07,  5.2455475777e-03,  4.3760834645e-06,  5.1394909319e-04,  1.3051561248e-03,  0.0000000000e+00],
                [ 0.0000000000e+00,  0.0000000000e+00,  0.0000000000e+00,  0.0000000000e+00,  0.0000000000e+00,  0.0000000000e+00,  0.0000000000e+00]
                ]),
        }, max=len(self.M))

class testComplete(unittest.TestCase):
  'Tests of entire lattice files'

  def test_LS1_cs0(self):
    "Propagate charge state 0 through LS1."
    PS_X = 0; PS_PX = 1; PS_Y = 2; PS_PY = 3; PS_S = 4; PS_PS = 5

    MeVtoeV = 1e6

    def sqr(a): return a*a;

    P = GLPSParser()

    file_name = 'LS1.lat'

    with open(os.path.join(datadir, file_name), 'rb') as inf:
      M = Machine(inf, extra={'cstate':0})

      S = M.allocState({})
      M.propagate(S, 0, 1)

      self.assertAlmostEqual(S.real_IonZ, 0.13865546218487396, 14)
      self.assertAlmostEqual(S.real_IonEs, 931494320.0, 14)

      self.assertAlmostEqual(S.ref_IonZ, 0.13865546218487396, 14)
      self.assertAlmostEqual(S.ref_IonEs, 931494320.0, 14)
      self.assertAlmostEqual(S.ref_IonEk, 500000.0, 14)

      def checkConsist(self, S, P='real'):
          self.assertEqual(getattr(S, P+'_IonW')    , getattr(S, P+'_IonEs')+getattr(S, P+'_IonEk'))
          self.assertEqual(getattr(S, P+'_gamma')   , getattr(S, P+'_IonW')/getattr(S, P+'_IonEs'))
          self.assertEqual(getattr(S, P+'_beta')    , math.sqrt(1e0-1e0/sqr(getattr(S, P+'_gamma'))))
          self.assertEqual(getattr(S, P+'_bg')      , getattr(S, P+'_beta')*getattr(S, P+'_gamma'))

      checkConsist(self, S, 'real')
      checkConsist(self, S, 'ref')
      self.assertAlmostEqual(S.ref_phis,  0.0, 14)

      assert_aequal(S.moment0_env,
        [-7.88600000e-04, 1.08371000e-05, 1.33734300e-02, 6.67853400e-06, -1.84772900e-04, 3.09995000e-04, 1.00000000e+00],
      decimal=8)

      assert_aequal(S.moment1_env, [
        [ 2.7630945017e+00, -4.2824733660e-04,  1.5817856917e-02,  2.1559419100e-05,  1.8638050602e-04, -2.9939448697e-05,  0.0000000000e+00],
        [-4.2824733660e-04,  3.8494660679e-06, -1.3838544007e-06, -1.8540987783e-08,  1.0677795224e-07,  5.2856401597e-09,  0.0000000000e+00],
        [ 1.5817856917e-02, -1.3838544007e-06,  2.3625122600e+00, -6.6932045998e-04, -5.8009993858e-04,  6.7165153406e-06,  0.0000000000e+00],
        [ 2.1559419100e-05, -1.8540987783e-08, -6.6932045998e-04,  4.8971138918e-06, -5.0161488223e-07,  5.5748422180e-08,  0.0000000000e+00],
        [ 1.8638050602e-04,  1.0677795224e-07, -5.8009993858e-04, -5.0161488223e-07,  6.7168718453e-04, -1.2322334153e-05,  0.0000000000e+00],
        [-2.9939448697e-05,  5.2856401597e-09,  6.7165153406e-06,  5.5748422180e-08, -1.2322334153e-05,  1.9952466899e-06,  0.0000000000e+00],
        [ 0.0000000000e+00,  0.0000000000e+00,  0.0000000000e+00,  0.0000000000e+00,  0.0000000000e+00,  0.0000000000e+00,  0.0000000000e+00]
      ], decimal=8)

      M.propagate(S, 1, len(M))

      self.assertAlmostEqual(S.real_IonZ, 0.13865546218487396, 12)
      self.assertAlmostEqual(S.real_IonEs, 931494320.0, 12)

      self.assertAlmostEqual(S.ref_IonZ, 0.13865546218487396, 12)
      self.assertAlmostEqual(S.ref_IonEs, 931494320.0, 12)

      checkConsist(self, S, 'real')
      self.assertAlmostEqual(S.real_phis,  1532.1994551432952, 14)
      self.assertAlmostEqual(S.real_IonEk,  17089939.45941341, 14)

      checkConsist(self, S, 'ref')
      self.assertAlmostEqual(S.ref_bg,  0.19243502172784563, 14)
      self.assertAlmostEqual(S.ref_phis,  1532.2018533339335, 14)
      self.assertAlmostEqual(S.ref_IonEk,  17090412.218117952, 14)

      assert_aequal(S.moment0_env,
        [-1.31487198e+00, -6.14116226e-04, -7.31682439e-01, -7.46618429e-04, -2.39819002e-03, -4.72758705e-04,  1.00000000e+00],
        decimal=8)

      assert_aequal(S.moment1_env, [
        [ 1.28985466e+00,  5.25768610e-04, -1.07818343e-01, -8.75275442e-05, -1.74466591e-04, -4.89638797e-05,  0.00000000e+00],
        [ 5.25768610e-04,  4.63230951e-07, -1.87620874e-04, -4.61354166e-08,  4.10553482e-08,  1.86305659e-08,  0.00000000e+00],
        [-1.07818343e-01, -1.87620874e-04,  2.54172968e+00,  4.50707865e-04,  2.29903897e-04,  1.14226041e-04,  0.00000000e+00],
        [-8.75275442e-05, -4.61354166e-08,  4.50707865e-04,  2.09807522e-07,  1.67123574e-07,  7.00281218e-08,  0.00000000e+00],
        [-1.74466591e-04,  4.10553482e-08,  2.29903897e-04,  1.67123574e-07,  3.45700733e-04,  1.29370100e-04,  0.00000000e+00],
        [-4.89638797e-05,  1.86305659e-08,  1.14226041e-04,  7.00281218e-08,  1.29370100e-04,  5.18511035e-05,  0.00000000e+00],
        [ 0.00000000e+00,  0.00000000e+00,  0.00000000e+00,  0.00000000e+00,  0.00000000e+00,  0.00000000e+00,  0.00000000e+00]
      ], decimal=8)


  def test_LS1_cs1(self):
      "Propagate charge state 1 through LS1."
      PS_X = 0; PS_PX = 1; PS_Y = 2; PS_PY = 3; PS_S = 4; PS_PS = 5

      MeVtoeV = 1e6

      def sqr(a): return a*a;

      P = GLPSParser()

      file_name = 'LS1.lat'

      with open(os.path.join(datadir, file_name), 'rb') as inf:
        M = Machine(inf, extra={'cstate':1})

        S = M.allocState({})
        M.propagate(S, 0, 1)

        self.assertAlmostEqual(S.real_IonZ, 0.14285714285714285, 14)
        self.assertAlmostEqual(S.real_IonEs, 931494320.0, 14)

        self.assertAlmostEqual(S.ref_IonZ, 0.13865546218487396, 14)
        self.assertAlmostEqual(S.ref_IonEs, 931494320.0, 14)
        self.assertAlmostEqual(S.ref_IonEk, 500000.0, 14)

        def checkConsist(self, S, P='real'):
            self.assertEqual(getattr(S, P+'_IonW')    , getattr(S, P+'_IonEs')+getattr(S, P+'_IonEk'))
            self.assertEqual(getattr(S, P+'_gamma')   , getattr(S, P+'_IonW')/getattr(S, P+'_IonEs'))
            self.assertEqual(getattr(S, P+'_beta')    , math.sqrt(1e0-1e0/sqr(getattr(S, P+'_gamma'))))
            self.assertEqual(getattr(S, P+'_bg')      , getattr(S, P+'_beta')*getattr(S, P+'_gamma'))

        checkConsist(self, S, 'real')
        checkConsist(self, S, 'ref')
        self.assertAlmostEqual(S.ref_phis,  0.0, 14)

        assert_aequal(S.moment0_env,
          [7.31266700e-03,   1.47656500e-05,   3.44984900e-03, -7.39768469e-06,   2.29788600e-02,   2.06010000e-03, 1.00000000e+00],
        decimal=8)

        assert_aequal(S.moment1_env, [
          [ 2.7932384438e+00,  8.8728620447e-04, -2.0796518392e-02, -4.6319090287e-05,  7.7964408668e-04, -1.2165264951e-05,  0.0000000000e+00],
          [ 8.8728620447e-04,  4.2768582732e-06,  1.9363966143e-05, -2.6132244852e-08,  5.4267659715e-07, -8.6693752768e-09,  0.0000000000e+00],
          [-2.0796518392e-02,  1.9363966143e-05,  2.7170190301e+00,  1.1114753191e-03, -1.0378823624e-03, -2.4536315081e-06,  0.0000000000e+00],
          [-4.6319090287e-05, -2.6132244852e-08,  1.1114753191e-03,  4.4323845352e-06, -1.0715375602e-06, -7.6510496305e-09,  0.0000000000e+00],
          [ 7.7964408668e-04,  5.4267659715e-07, -1.0378823624e-03, -1.0715375602e-06,  1.1048374735e-03,  1.1686344655e-05,  0.0000000000e+00],
          [-1.2165264951e-05, -8.6693752768e-09, -2.4536315081e-06, -7.6510496305e-09,  1.1686344655e-05,  1.7419755309e-06,  0.0000000000e+00],
          [ 0.0000000000e+00,  0.0000000000e+00,  0.0000000000e+00,  0.0000000000e+00,  0.0000000000e+00,  0.0000000000e+00,  0.0000000000e+00]
        ], decimal=8)

        M.propagate(S, 1, len(M))

        self.assertAlmostEqual(S.real_IonZ, 0.14285714285714285, 12)
        self.assertAlmostEqual(S.real_IonEs, 931494320.0, 12)

        self.assertAlmostEqual(S.ref_IonZ, 0.13865546218487396, 12)
        self.assertAlmostEqual(S.ref_IonEs, 931494320.0, 12)

        checkConsist(self, S, 'real')
        self.assertAlmostEqual(S.real_phis,  1532.222661902655, 14)
        self.assertAlmostEqual(S.real_IonEk,  17087328.436305404, 14)

        checkConsist(self, S, 'ref')
        self.assertAlmostEqual(S.ref_bg,  0.19243502172784563, 14)
        self.assertAlmostEqual(S.ref_phis,  1532.2018533339335, 14)
        self.assertAlmostEqual(S.ref_IonEk,  17090412.218117952, 14)

        assert_aequal(S.moment0_env,
          [-1.9433786002e+00, 3.4278187693e-04, -1.8746179484e+00, -7.8544844580e-04, 2.0808595165e-02, -3.0837818125e-03, 1.0000000000e+00],
          decimal=8)

        assert_aequal(S.moment1_env, [
            [1.6617407565e+00, -1.1561293961e-04, 4.2375967312e-01, 5.2156789327e-06, -3.2654073339e-04, -1.6658632146e-04, 0.0000000000e+00],
            [-1.1561293961e-04, 2.0635559854e-07, 1.7108084373e-06, -3.3149516396e-08, 1.5559377703e-07, 8.0035719827e-08, 0.0000000000e+00],
            [4.2375967312e-01, 1.7108084373e-06, 2.5993321057e+00, -3.1005696288e-05, -4.2384040453e-04, -2.0736703233e-04, 0.0000000000e+00],
            [5.2156789327e-06, -3.3149516396e-08, -3.1005696288e-05, 1.3142616035e-07, 1.4856230934e-07, 7.0098541878e-08, 0.0000000000e+00],
            [-3.2654073339e-04, 1.5559377703e-07, -4.2384040453e-04, 1.4856230934e-07, 7.5949774421e-04, 3.9605913591e-04, 0.0000000000e+00],
            [-1.6658632146e-04, 8.0035719827e-08, -2.0736703233e-04, 7.0098541878e-08, 3.9605913591e-04, 2.0888918966e-04, 0.0000000000e+00],
            [0.0000000000e+00, 0.0000000000e+00, 0.0000000000e+00, 0.0000000000e+00, 0.0000000000e+00, 0.0000000000e+00, 0.0000000000e+00],
        ], decimal=8)


  def test_through_stripper(self):
      "Propagate charge states 0-1 to Charge Stripper, and 0-4 afterwards"
      PS_X = 0; PS_PX = 1; PS_Y = 2; PS_PY = 3; PS_S = 4; PS_PS = 5

      MeVtoeV = 1e6

      def sqr(a): return a*a;

      P = GLPSParser()

      file_name = 'to_strl.lat'

      with open(os.path.join(datadir, file_name), 'rb') as inf:
        M = Machine(inf)

        S = M.allocState({})
        M.propagate(S, 0, 1)

        self.assertAlmostEqual(S.real_IonZ, 0.13865546218487396, 14)
        self.assertAlmostEqual(S.real_IonEs, 931494320.0, 14)
        self.assertAlmostEqual(S.real_IonEk, 500309.99500000, 14)

        self.assertAlmostEqual(S.ref_IonZ, 0.13865546218487396, 14)
        self.assertAlmostEqual(S.ref_IonEs, 931494320.0, 14)
        self.assertAlmostEqual(S.ref_IonEk, 500000.0, 14)

        def checkConsist(self, S, P='real'):
          self.assertEqual(getattr(S, P+'_IonW')    , getattr(S, P+'_IonEs')+getattr(S, P+'_IonEk'))
          self.assertEqual(getattr(S, P+'_gamma')   , getattr(S, P+'_IonW')/getattr(S, P+'_IonEs'))
          self.assertEqual(getattr(S, P+'_beta')    , math.sqrt(1e0-1e0/sqr(getattr(S, P+'_gamma'))))
          self.assertEqual(getattr(S, P+'_bg')      , getattr(S, P+'_beta')*getattr(S, P+'_gamma'))

        checkConsist(self, S, 'real')
        checkConsist(self, S, 'ref')
        self.assertAlmostEqual(S.ref_phis,  0.0, 14)

        assert_aequal(S.moment0_env,
          [3.3444511955e-03,1.2841341839e-05,8.3106826155e-03,-5.0277881002e-07,1.1632697213e-02,1.2028520756e-03,1.0000000000e+00],
        decimal=8)

        assert_aequal(S.moment1_env, [
          [ 2.7784895410e+00,  2.4291087928e-04, -2.8819146986e-03, -1.3098890047e-05,  5.3594209519e-04, -1.7328482801e-05,  0.0000000000e+00],
          [ 2.4291087928e-04,  4.0675140732e-06,  9.1913901431e-06, -2.2427664520e-08,  3.5190231444e-07, -1.1570581361e-10,  0.0000000000e+00],
          [-2.8819146986e-03,  9.1913901431e-06,  2.5433968050e+00,  2.3922914294e-04, -8.7109112521e-04, -2.3018796436e-06,  0.0000000000e+00],
          [-1.3098890047e-05, -2.2427664520e-08,  2.3922914294e-04,  4.6600708351e-06, -8.7385463546e-07,  1.7247530508e-08,  0.0000000000e+00],
          [ 5.3594209519e-04,  3.5190231444e-07, -8.7109112521e-04, -8.7385463546e-07,  1.0267518940e-03,  1.0056757657e-05,  0.0000000000e+00],
          [-1.7328482801e-05, -1.1570581361e-10, -2.3018796436e-06,  1.7247530508e-08,  1.0056757657e-05,  2.6314343481e-06,  0.0000000000e+00],
          [ 0.0000000000e+00,  0.0000000000e+00,  0.0000000000e+00,  0.0000000000e+00,  0.0000000000e+00,  0.0000000000e+00,  0.0000000000e+00]
        ], decimal=8)

        M.propagate(S, 1, len(M))

        self.assertAlmostEqual(S.real_IonZ, 0.31932773109243695, 12)
        self.assertAlmostEqual(S.real_IonEs, 931494320.0, 12)

        self.assertAlmostEqual(S.ref_IonZ, 0.3277310924369748, 12)
        self.assertAlmostEqual(S.ref_IonEs, 931494320.0, 12)

        checkConsist(self, S, 'real')
        self.assertAlmostEqual(S.real_phis,  2023.6417753279932, 14)
        self.assertAlmostEqual(S.real_IonEk,  16829746.76773572, 14)

        checkConsist(self, S, 'ref')
        self.assertAlmostEqual(S.ref_phis,  2023.6512871901991, 14)
        self.assertAlmostEqual(S.ref_IonEk,  16816951.191958785, 14)

        assert_aequal(S.moment0_env,
            [1.0923868437e+00,3.5728723827e-04,1.4420500562e+00,-1.1959772575e-04,-1.7378583261e-02,-2.2999692554e-03,1.0000000000e+00],
          decimal=8)

        assert_aequal(S.moment1_env, [
          [ 4.3986359287e+00, -1.0616512564e-04,  1.4411487761e-02, -3.0296029747e-04, -1.8869234689e-02,  2.0393230824e-03,  0.0000000000e+00],
          [-1.0616512564e-04,  4.9570504271e-07,  3.0901853376e-04,  3.9626728990e-08,  8.2462663320e-07, -9.8476930606e-07,  0.0000000000e+00],
          [ 1.4411487761e-02,  3.0901853376e-04,  6.5727958190e+00,  2.7540054815e-03,  7.7069472866e-03,  5.2455475777e-03,  0.0000000000e+00],
          [-3.0296029747e-04,  3.9626728990e-08,  2.7540054815e-03,  1.8784036222e-06,  9.3079770177e-06,  4.3760834645e-06,  0.0000000000e+00],
          [-1.8869234689e-02,  8.2462663320e-07,  7.7069472866e-03,  9.3079770177e-06,  5.0800488894e-04,  5.1394909319e-04,  0.0000000000e+00],
          [ 2.0393230824e-03, -9.8476930606e-07,  5.2455475777e-03,  4.3760834645e-06,  5.1394909319e-04,  1.3051561248e-03,  0.0000000000e+00],
          [ 0.0000000000e+00,  0.0000000000e+00,  0.0000000000e+00,  0.0000000000e+00,  0.0000000000e+00,  0.0000000000e+00,  0.0000000000e+00]
        ], decimal=8)


  def test_to_end(self):
      "Propagate charge states 0-1 to Charge Stripper, and 0-4 to end of structure with mis-alignments."
      PS_X = 0; PS_PX = 1; PS_Y = 2; PS_PY = 3; PS_S = 4; PS_PS = 5

      MeVtoeV = 1e6

      def sqr(a): return a*a;

      P = GLPSParser()

      file_name = 'to_chg_str_err.lat'

      with open(os.path.join(datadir, file_name), 'rb') as inf:
        M = Machine(inf)

        S = M.allocState({})
        M.propagate(S, 0, 1)

        self.assertAlmostEqual(S.real_IonZ, 0.13865546218487396, 14)
        self.assertAlmostEqual(S.real_IonEs, 931494320.0, 14)
        self.assertAlmostEqual(S.real_IonEk, 500309.99500000, 14)

        self.assertAlmostEqual(S.ref_IonZ, 0.13865546218487396, 14)
        self.assertAlmostEqual(S.ref_IonEs, 931494320.0, 14)
        self.assertAlmostEqual(S.ref_IonEk, 500000.0, 14)

        def checkConsist(self, S, P='real'):
          self.assertEqual(getattr(S, P+'_IonW')    , getattr(S, P+'_IonEs')+getattr(S, P+'_IonEk'))
          self.assertEqual(getattr(S, P+'_gamma')   , getattr(S, P+'_IonW')/getattr(S, P+'_IonEs'))
          self.assertEqual(getattr(S, P+'_beta')    , math.sqrt(1e0-1e0/sqr(getattr(S, P+'_gamma'))))
          self.assertEqual(getattr(S, P+'_bg')      , getattr(S, P+'_beta')*getattr(S, P+'_gamma'))

        checkConsist(self, S, 'real')
        checkConsist(self, S, 'ref')
        self.assertAlmostEqual(S.ref_phis,  0.0, 14)

        assert_aequal(S.moment0_env,
          [3.3444511955e-03,1.2841341839e-05,8.3106826155e-03,-5.0277881002e-07,1.1632697213e-02,1.2028520756e-03,1.0000000000e+00],
        decimal=8)

        assert_aequal(S.moment1_env, [
            [ 2.7784895410e+00,  2.4291087928e-04, -2.8819146986e-03, -1.3098890047e-05,  5.3594209519e-04, -1.7328482801e-05,  0.0000000000e+00],
            [ 2.4291087928e-04,  4.0675140732e-06,  9.1913901431e-06, -2.2427664520e-08,  3.5190231444e-07, -1.1570581361e-10,  0.0000000000e+00],
            [-2.8819146986e-03,  9.1913901431e-06,  2.5433968050e+00,  2.3922914294e-04, -8.7109112521e-04, -2.3018796436e-06,  0.0000000000e+00],
            [-1.3098890047e-05, -2.2427664520e-08,  2.3922914294e-04,  4.6600708351e-06, -8.7385463546e-07,  1.7247530508e-08,  0.0000000000e+00],
            [ 5.3594209519e-04,  3.5190231444e-07, -8.7109112521e-04, -8.7385463546e-07,  1.0267518940e-03,  1.0056757657e-05,  0.0000000000e+00],
            [-1.7328482801e-05, -1.1570581361e-10, -2.3018796436e-06,  1.7247530508e-08,  1.0056757657e-05,  2.6314343481e-06,  0.0000000000e+00],
            [ 0.0000000000e+00,  0.0000000000e+00,  0.0000000000e+00,  0.0000000000e+00,  0.0000000000e+00,  0.0000000000e+00,  0.0000000000e+00]
        ], decimal=8)

        M.propagate(S, 1, len(M))

        self.assertAlmostEqual(S.real_IonZ, 0.319327731092, 12)
        self.assertAlmostEqual(S.real_IonEs, 931494320.0, 12)

        self.assertAlmostEqual(S.ref_IonZ, 0.327731092437, 12)
        self.assertAlmostEqual(S.ref_IonEs, 931494320.0, 12)

        checkConsist(self, S, 'real')
        self.assertAlmostEqual(S.real_phis,  2023.640231990024, 12)
        self.assertAlmostEqual(S.real_IonEk,  16827008.470873594284, 12)

        checkConsist(self, S, 'ref')
        self.assertAlmostEqual(S.ref_phis,  2023.651314326443, 12)
        self.assertAlmostEqual(S.ref_IonEk,  16816951.191958785057, 12)

        assert_aequal(S.moment0_env,
            [6.0625179097e+00,9.3784744367e-04,1.0685277284e+00,-2.2663813468e-03,-1.8084741428e-02,-2.2133294018e-03,1.0000000000e+00],
          decimal=8)

        assert_aequal(S.moment1_env, [
          [ 1.3441952388e+01,  4.1493100647e-03,  3.5581817834e+01,  1.9146148774e-02,  3.7791956156e-03, -7.2604434996e-03,  0.0000000000e+00],
          [ 4.1493100647e-03,  4.5228253867e-06,  5.0067911136e-03,  1.1327665007e-05, -2.6303767247e-06, -1.4406372295e-05,  0.0000000000e+00],
          [ 3.5581817834e+01,  5.0067911136e-03,  1.4559022426e+02,  5.0499060613e-02,  2.0825579323e-02,  6.3711134701e-03,  0.0000000000e+00],
          [ 1.9146148774e-02,  1.1327665007e-05,  5.0499060613e-02,  4.5027125399e-05, -4.5961351604e-06, -3.5715788183e-05,  0.0000000000e+00],
          [ 3.7791956156e-03, -2.6303767247e-06,  2.0825579323e-02, -4.5961351604e-06,  2.2832789853e-04,  4.9093691735e-04,  0.0000000000e+00],
          [-7.2604434996e-03, -1.4406372295e-05,  6.3711134701e-03, -3.5715788183e-05,  4.9093691735e-04,  1.3348204463e-03,  0.0000000000e+00],
          [ 0.0000000000e+00,  0.0000000000e+00,  0.0000000000e+00,  0.0000000000e+00,  0.0000000000e+00,  0.0000000000e+00,  0.0000000000e+00]
        ], decimal=8)
