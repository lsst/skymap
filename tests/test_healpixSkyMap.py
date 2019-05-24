import unittest
import lsst.geom as geom
import lsst.utils.tests
from helper import skyMapTestCase

try:
    import healpy
except Exception:
    healpy = None

from lsst.skymap.healpixSkyMap import HealpixSkyMap


class HealpixTestCase(skyMapTestCase.SkyMapTestCase):

    def setUp(self):
        if not healpy:
            self.skipTest("Missing healpy dependency.")

        config = HealpixSkyMap.ConfigClass()
        nside = 2**config.log2NSide
        self.setAttributes(
            SkyMapClass=HealpixSkyMap,
            name="healpix",
            config=config,
            numTracts=healpy.nside2npix(nside),
            numNeighbors=1,
            neighborAngularSeparation=healpy.max_pixrad(nside) * geom.radians,
        )

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
