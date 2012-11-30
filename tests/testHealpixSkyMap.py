#!/usr/bin/env python

import lsst.utils.tests as utilsTests
import unittest


try:
    import healpy
except:
    import sys
    print >>sys.stderr, "WARNING: not testing HealpixSkyMap because healpy can't be imported."
    sys.exit(0)


import lsst.afw.geom as afwGeom
from SkyMapTestCase import SkyMapTestCase
from lsst.skymap.healpixSkyMap import HealpixSkyMap

config = HealpixSkyMap.ConfigClass()
global nside
nside = 2**config.nSide

class HealpixTestCase(SkyMapTestCase):
    _NumTracts = healpy.nside2npix(nside) # Number of tracts to expect
    _NeighborAngularSeparation = healpy.max_pixrad(nside) * afwGeom.degrees # Expected tract separation
    _SkyMapClass = HealpixSkyMap # Class of SkyMap to test
    _SkyMapName = "healpix" # Name of SkyMap class to test
    _numNeighbors = 1 # Number of neighbours



def suite():
    """Return a suite containing all the test cases in this module.
    """
    utilsTests.init()

    suites = [
        unittest.makeSuite(HealpixTestCase),
        unittest.makeSuite(utilsTests.MemoryTestCase),
    ]

    return unittest.TestSuite(suites)


def run(shouldExit=False):
    """Run the tests"""
    utilsTests.run(suite(), shouldExit)

if __name__ == "__main__":
    run(True)
