
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
            'moment0':asfarray([3.3444511955e-03,1.2841341839e-05,8.3106826155e-03,-5.0277881002e-07,1.1632697213e-02,1.2028520756e-03,1.0000000000e+00]),
            'state':asfarray([
                [  2.7784731401e+00,  2.4290292602e-04, -2.8618246244e-03, -1.3070393048e-05,  4.8904782352e-04, -2.0871532360e-05,  0.0000000000e+00,],
                [  2.4290292602e-04,  4.0675102164e-06,  9.2011324292e-06, -2.2413845461e-08,  3.2916186027e-07, -1.8338379704e-09,  0.0000000000e+00,],
                [ -2.8618246244e-03,  9.2011324292e-06,  2.5433721958e+00,  2.3919423578e-04, -8.1364836923e-04,  2.0381500435e-06,  0.0000000000e+00,],
                [ -1.3070393048e-05, -2.2413845461e-08,  2.3919423578e-04,  4.6600213206e-06, -7.9237429129e-07,  2.3403696009e-08,  0.0000000000e+00,],
                [  4.8904782352e-04,  3.2916186027e-07, -8.1364836923e-04, -7.9237429129e-07,  8.9266895435e-04, -7.3744068470e-08,  0.0000000000e+00,],
                [ -2.0871532360e-05, -1.8338379704e-09,  2.0381500435e-06,  2.3403696009e-08, -7.3744068470e-08,  1.8660344733e-06,  0.0000000000e+00,],
                [  0.0000000000e+00,  0.0000000000e+00,  0.0000000000e+00,  0.0000000000e+00,  0.0000000000e+00,  0.0000000000e+00,  0.0000000000e+00,],
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
            'moment0':asfarray([4.2690278079e-03,  1.2841341839e-05,  8.2744825412e-03, -5.0277881002e-07,  7.1994323704e-03,  1.2028520756e-03,  1.0000000000e+00]),
            'state':asfarray([
                [  2.8345371345e+00,  5.3576366160e-04, -3.2566047638e-03, -1.4684189921e-05,  5.9038764541e-04, -2.1003568694e-05,  0.0000000000e+00,],
                [  5.3576366160e-04,  4.0675102164e-06,  7.5873355560e-06, -2.2413845461e-08,  3.3586613225e-07, -1.8338379704e-09,  0.0000000000e+00,],
                [ -3.2566047638e-03,  7.5873355560e-06,  2.6019737163e+00,  5.7471577086e-04, -8.8451430669e-04,  3.7232161561e-06,  0.0000000000e+00,],
                [ -1.4684189921e-05, -2.2413845461e-08,  5.7471577086e-04,  4.6600213206e-06, -8.7910277532e-07,  2.3403696009e-08,  0.0000000000e+00,],
                [  5.9038764541e-04,  3.3586613225e-07, -8.8451430669e-04, -8.7910277532e-07,  9.1889975064e-04, -6.9656039192e-06,  0.0000000000e+00,],
                [ -2.1003568694e-05, -1.8338379704e-09,  3.7232161561e-06,  2.3403696009e-08, -6.9656039192e-06,  1.8660344733e-06,  0.0000000000e+00,],
                [  0.0000000000e+00,  0.0000000000e+00,  0.0000000000e+00,  0.0000000000e+00,  0.0000000000e+00,  0.0000000000e+00,  0.0000000000e+00,],
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
            'moment0':asfarray([9.8156783003e-03,  2.0041768075e-05, -5.2898622977e-03, -1.8638562712e-04, -1.7223152071e-02,  1.7187912211e-03,  1.0000000000e+00]),
            'state':asfarray([
                [  4.5737594920e+00,  6.0669283133e-03, -1.0742727858e-02, -3.1666439717e-05,  9.7062495255e-04,  2.1806021168e-05,  0.0000000000e+00,],
                [  6.0669283133e-03,  1.0279484098e-05, -1.0674677182e-05, -4.6435893047e-08,  1.1295245345e-06,  3.1617597182e-08,  0.0000000000e+00,],
                [ -1.0742727858e-02, -1.0674677182e-05,  4.4608321029e+00,  6.1805540423e-03, -1.2275842312e-03, -4.9248303794e-05,  0.0000000000e+00,],
                [ -3.1666439717e-05, -4.6435893047e-08,  6.1805540423e-03,  1.0960884518e-05, -1.9968111429e-06, -6.6642890830e-08,  0.0000000000e+00,],
                [  9.7062495255e-04,  1.1295245345e-06, -1.2275842312e-03, -1.9968111429e-06,  1.1512454713e-03,  1.7573680559e-05,  0.0000000000e+00,],
                [  2.1806021168e-05,  3.1617597182e-08, -4.9248303794e-05, -6.6642890830e-08,  1.7573680559e-05,  1.7361753083e-06,  0.0000000000e+00,],
                [  0.0000000000e+00,  0.0000000000e+00,  0.0000000000e+00,  0.0000000000e+00,  0.0000000000e+00,  0.0000000000e+00,  0.0000000000e+00,],
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
            'moment0':asfarray([-5.8659453275e-03, -4.0712274747e-05, -5.9179321956e-02, -1.1742389615e-04, -3.9532448827e-02,  1.7187912211e-03,  1.0000000000e+00]),
            'state':asfarray([
                [  8.0532807014e+00, -1.1996731464e-03, -1.2051303773e-02,  7.8258361648e-05,  5.0089780129e-04,  5.0696746897e-06,  0.0000000000e+00,],
                [ -1.1996731464e-03,  1.4386078126e-06,  9.9250124183e-05, -3.4702178884e-09, -3.3745081731e-07,  7.6550916976e-11,  0.0000000000e+00,],
                [ -1.2051303773e-02,  9.9250124183e-05,  8.0924884168e+00, -9.7900498004e-04, -1.1109526590e-03, -7.1207029942e-05,  0.0000000000e+00,],
                [  7.8258361648e-05, -3.4702178884e-09, -9.7900498004e-04,  1.4431711475e-06,  7.4362886116e-09,  1.4338999955e-08,  0.0000000000e+00,],
                [  5.0089780129e-04, -3.3745081731e-07, -1.1109526590e-03,  7.4362886116e-09,  9.8726482203e-04, -5.0096703455e-06,  0.0000000000e+00,],
                [  5.0696746897e-06,  7.6550916976e-11, -7.1207029942e-05,  1.4338999955e-08, -5.0096703455e-06,  1.7361753083e-06,  0.0000000000e+00,],
                [  0.0000000000e+00,  0.0000000000e+00,  0.0000000000e+00,  0.0000000000e+00,  0.0000000000e+00,  0.0000000000e+00,  0.0000000000e+00,],
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
            'moment0':asfarray([-5.8659453275e-03, -4.0712274747e-05, -5.9179321956e-02, -1.1742389615e-04, -3.9532448827e-02,  1.7187912211e-03,  1.0000000000e+00]),
            'state':asfarray([
                [  8.0532807014e+00, -1.1996731464e-03, -1.2051303773e-02,  7.8258361648e-05,  5.0089780129e-04,  5.0696746897e-06,  0.0000000000e+00,],
                [ -1.1996731464e-03,  1.4386078126e-06,  9.9250124183e-05, -3.4702178884e-09, -3.3745081731e-07,  7.6550916976e-11,  0.0000000000e+00,],
                [ -1.2051303773e-02,  9.9250124183e-05,  8.0924884168e+00, -9.7900498004e-04, -1.1109526590e-03, -7.1207029942e-05,  0.0000000000e+00,],
                [  7.8258361648e-05, -3.4702178884e-09, -9.7900498004e-04,  1.4431711475e-06,  7.4362886116e-09,  1.4338999955e-08,  0.0000000000e+00,],
                [  5.0089780129e-04, -3.3745081731e-07, -1.1109526590e-03,  7.4362886116e-09,  9.8726482203e-04, -5.0096703455e-06,  0.0000000000e+00,],
                [  5.0696746897e-06,  7.6550916976e-11, -7.1207029942e-05,  1.4338999955e-08, -5.0096703455e-06,  1.7361753083e-06,  0.0000000000e+00,],
                [  0.0000000000e+00,  0.0000000000e+00,  0.0000000000e+00,  0.0000000000e+00,  0.0000000000e+00,  0.0000000000e+00,  0.0000000000e+00,],
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
            'moment0':asfarray([-9.3631759606e-02, -1.2341993322e-04, -2.9374835471e-01, -1.1199135936e-03, -3.1818121415e-02,  6.0555868280e-04,  1.0000000000e+00]),
            'state':asfarray([
                [  3.5741613005e+00,  3.9046794009e-03,  2.6304937215e-01,  1.8832388456e-04, -3.7280843396e-04, -2.4524514220e-05,  0.0000000000e+00,],
                [  3.9046794009e-03,  5.3876891115e-06,  2.4887009931e-04,  1.4930685430e-07, -5.0692192162e-07, -1.6096969491e-08,  0.0000000000e+00,],
                [  2.6304937215e-01,  2.4887009931e-04,  3.6366600765e+00,  3.2873705063e-03, -2.2761267004e-04, -4.3929837302e-06,  0.0000000000e+00,],
                [  1.8832388456e-04,  1.4930685430e-07,  3.2873705063e-03,  4.0656323075e-06, -1.4318751121e-07,  3.8124461919e-08,  0.0000000000e+00,],
                [ -3.7280843396e-04, -5.0692192162e-07, -2.2761267004e-04, -1.4318751121e-07,  4.6797778782e-04,  2.8620197873e-05,  0.0000000000e+00,],
                [ -2.4524514220e-05, -1.6096969491e-08, -4.3929837302e-06,  3.8124461919e-08,  2.8620197873e-05,  5.0405336231e-06,  0.0000000000e+00,],
                [  0.0000000000e+00,  0.0000000000e+00,  0.0000000000e+00,  0.0000000000e+00,  0.0000000000e+00,  0.0000000000e+00,  0.0000000000e+00,],
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
            'moment0':asfarray([2.2103713766e+00,  3.3787339937e-03,  5.6689616131e-01, -1.9481100199e-04, -1.7814120019e-02, -5.4391793335e-03,  1.0000000000e+00]),
            'state':asfarray([
                [  2.6117833554e+00,  3.9301215251e-03,  1.4132329362e-01, -2.7122411451e-04,  6.8350668826e-05, -1.2074736723e-04,  0.0000000000e+00,],
                [  3.9301215251e-03,  6.0428894076e-06,  2.1668291365e-04, -4.1551664781e-07,  9.3707366679e-08, -1.6000265076e-07,  0.0000000000e+00,],
                [  1.4132329362e-01,  2.1668291365e-04,  1.3069340982e+00, -1.8336813219e-03,  2.5149058968e-05, -1.0935064333e-05,  0.0000000000e+00,],
                [ -2.7122411451e-04, -4.1551664781e-07, -1.8336813219e-03,  2.8944502415e-06, -6.4748711286e-08,  7.5438718878e-08,  0.0000000000e+00,],
                [  6.8350668826e-05,  9.3707366679e-08,  2.5149058968e-05, -6.4748711286e-08,  9.1409950081e-05, -3.5855445481e-05,  0.0000000000e+00,],
                [ -1.2074736723e-04, -1.6000265076e-07, -1.0935064333e-05,  7.5438718878e-08, -3.5855445481e-05,  4.2438739822e-05,  0.0000000000e+00,],
                [  0.0000000000e+00,  0.0000000000e+00,  0.0000000000e+00,  0.0000000000e+00,  0.0000000000e+00,  0.0000000000e+00,  0.0000000000e+00,],
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
            'moment0':asfarray([-1.6829029034e+00,  1.5337535377e-04, -1.6541522987e+00, -7.6642847191e-04,  9.9069891853e-03, -1.8048333267e-03,  1.0000000000e+00]),
            'state':asfarray([
                [  1.7200334179e+00,  3.3695570204e-04,  9.7502901814e-02, -5.7700135591e-05, -4.6850246200e-04, -8.8876731227e-05,  0.0000000000e+00,],
                [  3.3695570204e-04,  3.2879861631e-07, -1.0902596761e-04, -3.9259507748e-08, -8.5183456878e-08, -1.7587644231e-08,  0.0000000000e+00,],
                [  9.7502901814e-02, -1.0902596761e-04,  2.7859177924e+00,  2.8014644925e-04, -4.7392409186e-05, -1.8817840516e-05,  0.0000000000e+00,],
                [ -5.7700135591e-05, -3.9259507748e-08,  2.8014644925e-04,  1.6981943386e-07,  1.5857493370e-07,  7.0064048240e-08,  0.0000000000e+00,],
                [ -4.6850246200e-04, -8.5183456878e-08, -4.7392409186e-05,  1.5857493370e-07,  4.9882862165e-04,  2.5058625210e-04,  0.0000000000e+00,],
                [ -8.8876731227e-05, -1.7587644231e-08, -1.8817840516e-05,  7.0064048240e-08,  2.5058625210e-04,  1.3196776300e-04,  0.0000000000e+00,],
                [  0.0000000000e+00,  0.0000000000e+00,  0.0000000000e+00,  0.0000000000e+00,  0.0000000000e+00,  0.0000000000e+00,  0.0000000000e+00,],
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
            'moment0':asfarray([1.1399399086e+00,3.5332829960e-04,1.6963530125e+00,1.2840348237e-05,-1.6985609790e-02,-1.7053590681e-03,1.0000000000e+00]),
            'state':asfarray([
                [  3.9107054878e+00,  1.4077269633e-04,  1.5946156191e-01, -5.0140752974e-04, -2.1492674747e-02, -2.4476354523e-03,  0.0000000000e+00,],
                [  1.4077269633e-04,  4.0972926582e-07,  2.6001935149e-04,  7.1090223550e-08,  4.7229988322e-07,  9.7388195980e-08,  0.0000000000e+00,],
                [  1.5946156191e-01,  2.6001935149e-04,  3.5878843101e+00,  1.2229326516e-03,  2.1620812170e-03,  7.9439996216e-04,  0.0000000000e+00,],
                [ -5.0140752974e-04,  7.1090223550e-08,  1.2229326516e-03,  1.0587018969e-06,  6.4904502593e-06,  8.9280373156e-07,  0.0000000000e+00,],
                [ -2.1492674747e-02,  4.7229988322e-07,  2.1620812170e-03,  6.4904502593e-06,  5.1204412908e-04,  5.1595713518e-04,  0.0000000000e+00,],
                [ -2.4476354523e-03,  9.7388195980e-08,  7.9439996216e-04,  8.9280373156e-07,  5.1595713518e-04,  1.2775972163e-03,  0.0000000000e+00,],
                [  0.0000000000e+00,  0.0000000000e+00,  0.0000000000e+00,  0.0000000000e+00,  0.0000000000e+00,  0.0000000000e+00,  0.0000000000e+00,],
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

      assert_aequal(S.moment0,
        [-7.88600000e-04, 1.08371000e-05, 1.33734300e-02, 6.67853400e-06, -1.84772900e-04, 3.09995000e-04, 1.00000000e+00],
      decimal=8)

      assert_aequal(S.state, [
        [ 2.76309450e+00, -4.28247337e-04,  1.58178569e-02,  2.15594191e-05,  1.86380506e-04, -2.99394487e-05,  0.00000000e+00],
        [-4.28247337e-04,  3.84946607e-06, -1.38385440e-06, -1.85409878e-08,  1.06777952e-07,  5.28564016e-09,  0.00000000e+00],
        [ 1.58178569e-02, -1.38385440e-06,  2.36251226e+00, -6.69320460e-04, -5.80099939e-04,  6.71651534e-06,  0.00000000e+00],
        [ 2.15594191e-05, -1.85409878e-08, -6.69320460e-04,  4.89711389e-06, -5.01614882e-07,  5.57484222e-08,  0.00000000e+00],
        [ 1.86380506e-04,  1.06777952e-07, -5.80099939e-04, -5.01614882e-07,  6.71687185e-04, -1.23223342e-05,  0.00000000e+00],
        [-2.99394487e-05,  5.28564016e-09,  6.71651534e-06,  5.57484222e-08, -1.23223342e-05,  1.99524669e-06,  0.00000000e+00],
        [ 0.00000000e+00,  0.00000000e+00,  0.00000000e+00,  0.00000000e+00,  0.00000000e+00,  0.00000000e+00,  0.00000000e+00]
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

      assert_aequal(S.moment0,
        [-1.31487198e+00, -6.14116226e-04, -7.31682439e-01, -7.46618429e-04, -2.39819002e-03, -4.72758705e-04,  1.00000000e+00],
        decimal=8)

      assert_aequal(S.state, [
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

        assert_aequal(S.moment0,
          [7.31266700e-03,   1.47656500e-05,   3.44984900e-03, -7.39768469e-06,   2.29788600e-02,   2.06010000e-03, 1.00000000e+00],
        decimal=8)

        assert_aequal(S.state, [
          [2.79323844382261, 0.000887286204468529, -0.0207965183919297, -4.63190902866223e-05, 0.000779644086677272, -1.21652649506989e-05, 0.0],
          [0.000887286204468529, 4.27685827324966e-06, 1.93639661427031e-05, -2.61322448516117e-08, 5.42676597149308e-07, -8.66937527683823e-09, 0.0],
          [-0.0207965183919297, 1.93639661427031e-05, 2.71701903011773, 0.00111147531913367, -0.00103788236242888, -2.45363150807461e-06, 0.0],
          [-4.63190902866223e-05, -2.61322448516117e-08, 0.00111147531913367, 4.43238453519665e-06, -1.07153756021073e-06, -7.65104963052334e-09, 0.0],
          [0.000779644086677272, 5.42676597149308e-07, -0.00103788236242888, -1.07153756021073e-06, 0.00110483747345619, 1.16863446546579e-05, 0.0],
          [-1.21652649506989e-05, -8.66937527683823e-09, -2.45363150807461e-06, -7.65104963052334e-09, 1.16863446546579e-05, 1.74197553089731e-06, 0.0],
          [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
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

        assert_aequal(S.moment0,
          [-1.9433786002e+00, 3.4278187693e-04, -1.8746179484e+00, -7.8544844580e-04, 2.0808595165e-02, -3.0837818125e-03, 1.0000000000e+00],
          decimal=8)

        assert_aequal(S.state, [
            [1.6617407565e+00, -1.1561293961e-04, 4.2375967312e-01, 5.2156789327e-06, -3.2654073339e-04, -1.6658632146e-04, 0.0000000000e+00],
            [-1.1561293961e-04, 2.0635559854e-07, 1.7108084373e-06, -3.3149516396e-08, 1.5559377703e-07, 8.0035719827e-08, 0.0000000000e+00],
            [4.2375967312e-01, 1.7108084373e-06, 2.5993321057e+00, -3.1005696288e-05, -4.2384040453e-04, -2.0736703233e-04, 0.0000000000e+00],
            [5.2156789327e-06, -3.3149516396e-08, -3.1005696288e-05, 1.3142616035e-07, 1.4856230934e-07, 7.0098541878e-08, 0.0000000000e+00],
            [-3.2654073339e-04, 1.5559377703e-07, -4.2384040453e-04, 1.4856230934e-07, 7.5949774421e-04, 3.9605913591e-04, 0.0000000000e+00],
            [-1.6658632146e-04, 8.0035719827e-08, -2.0736703233e-04, 7.0098541878e-08, 3.9605913591e-04, 2.0888918966e-04, 0.0000000000e+00],
            [0.0000000000e+00, 0.0000000000e+00, 0.0000000000e+00, 0.0000000000e+00, 0.0000000000e+00, 0.0000000000e+00, 0.0000000000e+00],
        ], decimal=8)


  def test_to_chg_str_cs0(self):
    "Propagate charge state 0 through up to Charge Stripper with mis-alignments."
    PS_X = 0; PS_PX = 1; PS_Y = 2; PS_PY = 3; PS_S = 4; PS_PS = 5

    MeVtoeV = 1e6

    def sqr(a): return a*a;

    P = GLPSParser()

    file_name = 'to_chg_str_err.lat'

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

      assert_aequal(S.moment0,
        [-7.88600000e-04, 1.08371000e-05, 1.33734300e-02, 6.67853400e-06, -1.84772900e-04, 3.09995000e-04, 1.00000000e+00],
      decimal=8)

      assert_aequal(S.state, [
        [ 2.76309450e+00, -4.28247337e-04,  1.58178569e-02,  2.15594191e-05,  1.86380506e-04, -2.99394487e-05,  0.00000000e+00],
        [-4.28247337e-04,  3.84946607e-06, -1.38385440e-06, -1.85409878e-08,  1.06777952e-07,  5.28564016e-09,  0.00000000e+00],
        [ 1.58178569e-02, -1.38385440e-06,  2.36251226e+00, -6.69320460e-04, -5.80099939e-04,  6.71651534e-06,  0.00000000e+00],
        [ 2.15594191e-05, -1.85409878e-08, -6.69320460e-04,  4.89711389e-06, -5.01614882e-07,  5.57484222e-08,  0.00000000e+00],
        [ 1.86380506e-04,  1.06777952e-07, -5.80099939e-04, -5.01614882e-07,  6.71687185e-04, -1.23223342e-05,  0.00000000e+00],
        [-2.99394487e-05,  5.28564016e-09,  6.71651534e-06,  5.57484222e-08, -1.23223342e-05,  1.99524669e-06,  0.00000000e+00],
        [ 0.00000000e+00,  0.00000000e+00,  0.00000000e+00,  0.00000000e+00,  0.00000000e+00,  0.00000000e+00,  0.00000000e+00]
      ], decimal=8)

      M.propagate(S, 1, len(M))

      self.assertAlmostEqual(S.real_IonZ, 0.13865546218487396, 12)
      self.assertAlmostEqual(S.real_IonEs, 931494320.0, 12)

      self.assertAlmostEqual(S.ref_IonZ, 0.13865546218487396, 12)
      self.assertAlmostEqual(S.ref_IonEs, 931494320.0, 12)

      checkConsist(self, S, 'real')
      self.assertAlmostEqual(S.real_phis,  1602.4717408629808, 14)
      self.assertAlmostEqual(S.real_IonEk,  17089939.459636927, 14)

      checkConsist(self, S, 'ref')
      self.assertAlmostEqual(S.ref_bg,  0.19243502172784563, 14)
      self.assertAlmostEqual(S.ref_phis,  1602.4730487563952, 14)
      self.assertAlmostEqual(S.ref_IonEk,  17090412.218117952, 14)

      assert_aequal(S.moment0,
         [ 2.90434086e+00,  1.19171594e-03, -3.88357134e+00, -1.49883368e-03, -1.58920581e-03, -3.94313801e-04,  1.00000000e+00],
        decimal=8)

      assert_aequal(S.state, [
         [ 8.91033527e-02, -9.67351585e-05,  3.19930426e-02,  1.46443418e-05, -1.15621071e-06, -1.10724626e-04,  0.00000000e+00],
         [-9.67351585e-05,  3.70578399e-06, -1.70784631e-04, -4.34728844e-07,  5.03306187e-08, -7.24623758e-07,  0.00000000e+00],
         [ 3.19930426e-02, -1.70784631e-04,  4.41132936e-01,  6.32587731e-04, -2.65947638e-04,  2.59042891e-04,  0.00000000e+00],
         [ 1.46443418e-05, -4.34728844e-07,  6.32587731e-04,  1.65571086e-06, -4.49804549e-07,  2.96003425e-07,  0.00000000e+00],
         [-1.15621071e-06,  5.03306187e-08, -2.65947638e-04, -4.49804549e-07,  3.57440825e-05,  2.54131690e-05,  0.00000000e+00],
         [-1.10724626e-04, -7.24623758e-07,  2.59042891e-04,  2.96003425e-07,  2.54131690e-05,  5.23153894e-05,  0.00000000e+00],
         [ 0.00000000e+00,  0.00000000e+00,  0.00000000e+00,  0.00000000e+00,  0.00000000e+00,  0.00000000e+00,  0.00000000e+00]
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

        assert_aequal(S.moment0,
          [3.3444511955e-03,1.2841341839e-05,8.3106826155e-03,-5.0277881002e-07,1.1632697213e-02,1.2028520756e-03,1.0000000000e+00],
        decimal=8)

        assert_aequal(S.state, [
            [  2.7784731401e+00,  2.4290292602e-04, -2.8618246244e-03, -1.3070393048e-05,  4.8904782352e-04, -2.0871532360e-05,  0.0000000000e+00,],
            [  2.4290292602e-04,  4.0675102164e-06,  9.2011324292e-06, -2.2413845461e-08,  3.2916186027e-07, -1.8338379704e-09,  0.0000000000e+00,],
            [ -2.8618246244e-03,  9.2011324292e-06,  2.5433721958e+00,  2.3919423578e-04, -8.1364836923e-04,  2.0381500435e-06,  0.0000000000e+00,],
            [ -1.3070393048e-05, -2.2413845461e-08,  2.3919423578e-04,  4.6600213206e-06, -7.9237429129e-07,  2.3403696009e-08,  0.0000000000e+00,],
            [  4.8904782352e-04,  3.2916186027e-07, -8.1364836923e-04, -7.9237429129e-07,  8.9266895435e-04, -7.3744068470e-08,  0.0000000000e+00,],
            [ -2.0871532360e-05, -1.8338379704e-09,  2.0381500435e-06,  2.3403696009e-08, -7.3744068470e-08,  1.8660344733e-06,  0.0000000000e+00,],
            [  0.0000000000e+00,  0.0000000000e+00,  0.0000000000e+00,  0.0000000000e+00,  0.0000000000e+00,  0.0000000000e+00,  0.0000000000e+00,],
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

        assert_aequal(S.moment0,
           [1.1399399086e+00,3.5332829960e-04,1.6963530125e+00,1.2840348237e-05,-1.6985609790e-02,-1.7053590681e-03,1.0000000000e+00],
          decimal=8)

        assert_aequal(S.state, [
            [  3.9107054878e+00,  1.4077269633e-04,  1.5946156191e-01, -5.0140752974e-04, -2.1492674747e-02, -2.4476354523e-03,  0.0000000000e+00,],
            [  1.4077269633e-04,  4.0972926582e-07,  2.6001935149e-04,  7.1090223550e-08,  4.7229988322e-07,  9.7388195980e-08,  0.0000000000e+00,],
            [  1.5946156191e-01,  2.6001935149e-04,  3.5878843101e+00,  1.2229326516e-03,  2.1620812170e-03,  7.9439996216e-04,  0.0000000000e+00,],
            [ -5.0140752974e-04,  7.1090223550e-08,  1.2229326516e-03,  1.0587018969e-06,  6.4904502593e-06,  8.9280373156e-07,  0.0000000000e+00,],
            [ -2.1492674747e-02,  4.7229988322e-07,  2.1620812170e-03,  6.4904502593e-06,  5.1204412908e-04,  5.1595713518e-04,  0.0000000000e+00,],
            [ -2.4476354523e-03,  9.7388195980e-08,  7.9439996216e-04,  8.9280373156e-07,  5.1595713518e-04,  1.2775972163e-03,  0.0000000000e+00,],
            [  0.0000000000e+00,  0.0000000000e+00,  0.0000000000e+00,  0.0000000000e+00,  0.0000000000e+00,  0.0000000000e+00,  0.0000000000e+00,],
        ], decimal=8)
