#!/usr/bin/env python
import unittest

import lsst.utils.tests

from lsst.skymap.ringsSkyMap import RingsSkyMap
from helper import skyMapTestCase


config = RingsSkyMap.ConfigClass()
config.numRings = 3


class RingsTestCase(skyMapTestCase.SkyMapTestCase):

    def setUp(self):
        s_cls = skyMapTestCase.SkyMapTestCase
        s_cls._NumTracts = 26  # Number of tracts to expect
        s_cls._NeighborAngularSeparation = None  # Expected tract separation
        s_cls._SkyMapClass = RingsSkyMap  # Class of SkyMap to test
        s_cls._SkyMapConfig = config  # Configuration to use
        s_cls._SkyMapName = "rings"  # Name of SkyMap class to test
        s_cls._numNeighbors = None  # Number of neighbours

    def testTractSeparation(self):
        self.skipTest("A particular tract separation is not important for RingsSkyMap")


class MemoryTester(lsst.utils.tests.MemoryTestCase):
    pass


def setup_module(module):
    lsst.utils.tests.init()


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
