#!/usr/bin/env python
from __future__ import print_function
import unittest
import lsst.afw.geom as afwGeom
import lsst.utils.tests
from helper import skyMapTestCase

try:
    import healpy
except:
    healpy = None

from lsst.skymap.healpixSkyMap import HealpixSkyMap


config = HealpixSkyMap.ConfigClass()
global nside
nside = 2**config.log2NSide


class HealpixTestCase(skyMapTestCase.SkyMapTestCase):

    def setUp(self):
        if not healpy:
            self.skipTest("Missing healpy dependency.")
        s_cls = skyMapTestCase.SkyMapTestCase
        s_cls._NumTracts = healpy.nside2npix(nside)  # Number of tracts to expect
        s_cls._NeighborAngularSeparation = healpy.max_pixrad(
            nside) * afwGeom.radians  # Expected tract separation
        s_cls._SkyMapClass = HealpixSkyMap  # Class of SkyMap to test
        s_cls._SkyMapName = "healpix"  # Name of SkyMap class to test
        s_cls._numNeighbors = 1  # Number of neighbours


class MemoryTester(lsst.utils.tests.MemoryTestCase):
    pass


def setup_module(module):
    lsst.utils.tests.init()


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
