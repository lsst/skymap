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


class HealpixTestCase(skyMapTestCase.SkyMapTestCase):

    def setUp(self):
        if not healpy:
            self.skipTest("Missing healpy dependency.")

        self.config = HealpixSkyMap.ConfigClass()
        nside = 2**self.config.log2NSide
        self._NumTracts = healpy.nside2npix(nside)  # Number of tracts to expect
        self._NeighborAngularSeparation = healpy.max_pixrad(
            nside) * afwGeom.radians  # Expected tract separation
        self._SkyMapClass = HealpixSkyMap  # Class of SkyMap to test
        self._SkyMapName = "healpix"  # Name of SkyMap class to test
        self._numNeighbors = 1  # Number of neighbours

    def testSha1Compare(self):
        """Test that HealpixSkyMap's extra state is included in its hash."""
        defaultSkyMap = self.getSkyMap()
        for log2NSide in (3, 4):
            config = self.getConfig()
            config.log2NSide = log2NSide
            skyMap = self.getSkyMap(config=config)
            self.assertNotEqual(skyMap, defaultSkyMap)
        for nest in (True,):
            config = self.getConfig()
            config.nest = nest
            skyMap = self.getSkyMap(config=config)
            self.assertNotEqual(skyMap, defaultSkyMap)

    def tearDown(self):
        if hasattr(self, "config"):
            del self.config


class MemoryTester(lsst.utils.tests.MemoryTestCase):
    pass


def setup_module(module):
    lsst.utils.tests.init()


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
