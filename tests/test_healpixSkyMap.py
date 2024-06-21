import unittest
import lsst.geom as geom
import lsst.utils.tests
from helper import skyMapTestCase

import hpgeom

from lsst.skymap.healpixSkyMap import HealpixSkyMap


# TODO: Remove with DM-44799
class HealpixTestCase(skyMapTestCase.SkyMapTestCase):

    def setUp(self):
        with self.assertWarns(FutureWarning):
            config = HealpixSkyMap.ConfigClass()
        nside = 2**config.log2NSide
        self.setAttributes(
            SkyMapClass=HealpixSkyMap,
            name="healpix",
            config=config,
            numTracts=hpgeom.nside_to_npixel(nside),
            numNeighbors=1,
            neighborAngularSeparation=hpgeom.max_pixel_radius(nside) * geom.degrees,
        )

    def getSkyMap(self, config=None):
        with self.assertWarns(FutureWarning):
            skymap = super().getSkyMap(config=config)
        return skymap

    def getConfig(self):
        with self.assertWarns(FutureWarning):
            config = super().getConfig()
        return config

    # def testRegistry(self):
    #     with self.assertWarns(FutureWarning):
    #         super().testRegistry()

    def testBasicAttributes(self):
        with self.assertWarns(FutureWarning):
            super().testBasicAttributes()

    def testPickle(self):
        with self.assertWarns(FutureWarning):
            super().testPickle()

    def testTractSeparation(self):
        with self.assertWarns(FutureWarning):
            super().testTractSeparation()

    def testFindPatchList(self):
        with self.assertWarns(FutureWarning):
            super().testFindPatchList()

    def testFindTractPatchList(self):
        with self.assertWarns(FutureWarning):
            super().testFindTractPatchList()

    def testTractContains(self):
        with self.assertWarns(FutureWarning):
            super().testTractContains()

    def testTractInfoGetPolygon(self):
        with self.assertWarns(FutureWarning):
            super().testTractInfoGetPolygon()

    def testPatchInfoGetPolygon(self):
        with self.assertWarns(FutureWarning):
            super().testPatchInfoGetPolygon()

    def testTractInfoGetRegion(self):
        with self.assertWarns(FutureWarning):
            super().testTractInfoGetRegion()

    def testDm14809(self):
        with self.assertWarns(FutureWarning):
            super().testDm14809()

    def testNumbering(self):
        with self.assertWarns(FutureWarning):
            super().testNumbering()

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
