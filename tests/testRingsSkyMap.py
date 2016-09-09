#!/usr/bin/env python

import lsst.utils.tests as utilsTests
import unittest

import lsst.afw.geom as afwGeom
from SkyMapTestCase import SkyMapTestCase
from lsst.skymap.ringsSkyMap import RingsSkyMap

config = RingsSkyMap.ConfigClass()
config.numRings = 3


class RingsTestCase(SkyMapTestCase):
    _NumTracts = 26 # Number of tracts to expect
    _NeighborAngularSeparation = None # Expected tract separation
    _SkyMapClass = RingsSkyMap # Class of SkyMap to test
    _SkyMapConfig = config # Configuration to use
    _SkyMapName = "rings" # Name of SkyMap class to test
    _numNeighbors = None # Number of neighbours

    def testTractSeparation(self):
        self.skipTest("A particular tract separation is not important for RingsSkyMap")


def suite():
    """Return a suite containing all the test cases in this module.
    """
    utilsTests.init()

    suites = [
        unittest.makeSuite(RingsTestCase),
        unittest.makeSuite(utilsTests.MemoryTestCase),
    ]

    return unittest.TestSuite(suites)


def run(shouldExit=False):
    """Run the tests"""
    utilsTests.run(suite(), shouldExit)

if __name__ == "__main__":
    run(True)
