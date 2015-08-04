from __future__ import print_function

import unittest

from .._internal import dictshow, dict2conf

class testGood(unittest.TestCase):
    data = [
      ({}, ""),
      ({'a':1}, "Name: 'a' type i\n"),
      ({'a':[1,2]}, ">>> a\ni\ni\n<<<\n"),
    ]

    def test_good(self):
        for I, E in self.data:
            try:
                A = dictshow(dict2conf(I))
                self.assertEqual(A, E)
            except:
                print('Error on', I, E)
                raise

    def test_fail(self):
        # argument not a {}
        self.assertRaises(Exception, dict2conf, None)
        self.assertRaises(Exception, dict2conf, set())

        # unknown type
        self.assertRaises(KeyError, dict2conf, {'a':None})
        self.assertRaises(KeyError, dict2conf, {'a':set()})
