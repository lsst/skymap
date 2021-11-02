import unittest

import lsst.utils.tests
import lsst.geom
import lsst.sphgeom

from lsst.skymap.discreteSkyMap import DiscreteSkyMap
from lsst.skymap import Index2D


class PatchInfoTestCase(lsst.utils.tests.TestCase):
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
        patchInfo = tractInfo[0]

        self.assertEqual(patchInfo.index, patchInfo.getIndex())
        self.assertEqual(patchInfo.inner_bbox, patchInfo.getInnerBBox())
        self.assertEqual(patchInfo.outer_bbox, patchInfo.getOuterBBox())
        self.assertEqual(patchInfo.num_cells, patchInfo.getNumCells())
        self.assertEqual(patchInfo.cell_border, patchInfo.getCellBorder())
        self.assertEqual(patchInfo.cell_inner_dimensions, patchInfo.getCellInnerDimensions())
        self.assertEqual(patchInfo.wcs, patchInfo.getWcs())
        self.assertEqual(patchInfo.wcs, tractInfo.wcs)
        self.assertEqual(patchInfo.sequential_index, patchInfo.getSequentialIndex())
        self.assertEqual(patchInfo.inner_sky_polygon,
                         patchInfo.getInnerSkyPolygon(tractWcs=tractInfo.wcs))
        self.assertEqual(patchInfo.outer_sky_polygon,
                         patchInfo.getOuterSkyPolygon(tractWcs=tractInfo.wcs))

    def testIndexes(self):
        """Test accessing by index via tuple and Index2D"""
        tractInfo = self.skymap[0]

        index = Index2D(x=2, y=5)
        sequential_index = tractInfo.getSequentialPatchIndexFromPair(index)
        patchInfo = tractInfo[index]

        self.assertEqual(tractInfo.getSequentialPatchIndexFromPair((2, 5)), sequential_index)
        self.assertEqual(tractInfo.getPatchInfo(index), patchInfo)
        self.assertEqual(tractInfo.getPatchInfo((2, 5)), patchInfo)
        self.assertEqual(tractInfo[index], patchInfo)
        self.assertEqual(tractInfo[sequential_index], patchInfo)
        self.assertEqual(tractInfo[(2, 5)], patchInfo)

        for patchInfo in tractInfo:
            self.assertEqual(tractInfo[patchInfo.index], patchInfo)

    def testSequentialIndex(self):
        """Test patch sequential indices."""
        tractInfo = self.skymap[0]

        for patchInfo in tractInfo:
            self.assertEqual(patchInfo.sequential_index,
                             tractInfo.getSequentialPatchIndex(patchInfo))

    def testSkyPolygons(self):
        """Test sky polygons."""
        tractInfo = self.skymap[0]
        patchInfo = tractInfo[0]

        wcs = tractInfo.wcs

        inner_bbox = patchInfo.inner_bbox
        points = wcs.pixelToSky(lsst.geom.Box2D(inner_bbox).getCorners())
        poly = lsst.sphgeom.ConvexPolygon([p.getVector() for p in points])

        self.assertEqual(patchInfo.inner_sky_polygon, poly)

    def testNoCells(self):
        """Test retrieving cells for a legacy tract builder."""
        tractInfo = self.skymap[0]
        patchInfo = tractInfo[0]

        self.assertRaises(IndexError, patchInfo.__getitem__, 0)
        self.assertRaises(IndexError, patchInfo.__getitem__, 10)

        self.assertRaises(IndexError, patchInfo.getCellIndexPair, 0)
        self.assertRaises(IndexError, patchInfo.getCellIndexPair, 10)


class MemoryTester(lsst.utils.tests.MemoryTestCase):
    pass


def setup_module(module):
    lsst.utils.tests.init()


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
