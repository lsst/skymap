import unittest
import numpy as np

import lsst.utils.tests
import lsst.geom

from lsst.skymap.discreteSkyMap import DiscreteSkyMap


class TractInfoTestCase(lsst.utils.tests.TestCase):
    def setUp(self):
        # These are from the PS1 Medium-Deep fields
        coords = [(10.6750, 41.2667),  # M31
                  (36.2074, -04.5833),  # XMM-LSS
                  ]
        config = DiscreteSkyMap.ConfigClass()
        config.raList = [c[0] for c in coords]
        config.decList = [c[1] for c in coords]
        config.radiusList = [2] * len(coords)
        config.validate()

        self.skymap = DiscreteSkyMap(config=config)

    def testProperties(self):
        """Test accessing via properties."""
        tractInfo = self.skymap[0]

        self.assertEqual(tractInfo.bbox, tractInfo.getBBox())
        self.assertEqual(tractInfo.ctr_coord, tractInfo.getCtrCoord())
        self.assertEqual(tractInfo.tract_id, tractInfo.getId())
        self.assertEqual(tractInfo.num_patches, tractInfo.getNumPatches())
        self.assertEqual(tractInfo.patch_border, tractInfo.getPatchBorder())
        self.assertEqual(tractInfo.patch_inner_dimensions, tractInfo.getPatchInnerDimensions())
        self.assertEqual(tractInfo.tract_overlap, tractInfo.getTractOverlap())
        self.assertEqual(tractInfo.vertex_list, tractInfo.getVertexList())
        self.assertEqual(tractInfo.inner_sky_polygon, tractInfo.getInnerSkyPolygon())
        self.assertEqual(tractInfo.outer_sky_polygon, tractInfo.getOuterSkyPolygon())
        self.assertEqual(tractInfo.wcs, tractInfo.getWcs())

    def testContains(self):
        """Simple tests for the contains method, including bad inputs."""
        tractInfo = self.skymap[0]
        coord = tractInfo.ctr_coord

        self.assertTrue(tractInfo.contains(coord))

        badCoord = lsst.geom.SpherePoint(np.nan*lsst.geom.degrees, 0.0*lsst.geom.degrees)
        self.assertFalse(tractInfo.contains(badCoord))

        badCoord = lsst.geom.SpherePoint(0.0*lsst.geom.degrees, np.nan*lsst.geom.degrees)
        self.assertFalse(tractInfo.contains(badCoord))


class MemoryTester(lsst.utils.tests.MemoryTestCase):
    pass


def setup_module(module):
    lsst.utils.tests.init()


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
