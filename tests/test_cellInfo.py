import unittest

import lsst.utils.tests
import lsst.geom

from lsst.skymap.discreteSkyMap import DiscreteSkyMap
from lsst.skymap import Index2D


class CellInfoTestCase(lsst.utils.tests.TestCase):
    def setUp(self):
        # These are from the PS1 Medium-Deep fields
        coords = [(10.6750, 41.2667),  # M31
                  (36.2074, -04.5833),  # XMM-LSS
                  ]
        config = DiscreteSkyMap.ConfigClass()
        config.raList = [c[0] for c in coords]
        config.decList = [c[1] for c in coords]
        config.radiusList = [2] * len(coords)
        config.tractBuilder = "cells"
        config.validate()

        self.skymap = DiscreteSkyMap(config=config)

    def testProperties(self):
        """Test accessing via properties."""
        tractInfo = self.skymap[0]
        patchInfo = tractInfo[0]
        cellInfo = patchInfo[0]

        self.assertEqual(cellInfo.index, cellInfo.getIndex())
        self.assertEqual(cellInfo.inner_bbox, cellInfo.getInnerBBox())
        self.assertEqual(cellInfo.outer_bbox, cellInfo.getOuterBBox())
        self.assertEqual(cellInfo.wcs, cellInfo.getWcs())
        self.assertEqual(cellInfo.wcs, tractInfo.wcs)
        self.assertEqual(cellInfo.sequential_index, cellInfo.getSequentialIndex())
        self.assertEqual(cellInfo.inner_sky_polygon,
                         cellInfo.getInnerSkyPolygon(tractWcs=tractInfo.wcs))
        self.assertEqual(cellInfo.outer_sky_polygon,
                         cellInfo.getOuterSkyPolygon(tractWcs=tractInfo.wcs))

    def testIndexes(self):
        """Test accessing by index via tuple and Index2D"""
        tractInfo = self.skymap[0]
        patchInfo = tractInfo[0]

        index = Index2D(x=5, y=2)
        sequential_index = patchInfo.getSequentialCellIndexFromPair(index)
        cellInfo = patchInfo[index]

        self.assertEqual(patchInfo.getSequentialCellIndexFromPair((5, 2)), sequential_index)
        self.assertEqual(patchInfo.getCellInfo(index), cellInfo)
        self.assertEqual(patchInfo.getCellInfo((5, 2)), cellInfo)
        self.assertEqual(patchInfo[index], cellInfo)
        self.assertEqual(patchInfo[sequential_index], cellInfo)
        self.assertEqual(patchInfo[(5, 2)], cellInfo)

        for cellInfo in patchInfo:
            self.assertEqual(patchInfo[cellInfo.index], cellInfo)

    def testSequentialIndex(self):
        """Test cell sequential indices."""
        tractInfo = self.skymap[0]
        patchInfo = tractInfo[0]

        for cellInfo in patchInfo:
            self.assertEqual(cellInfo.sequential_index,
                             patchInfo.getSequentialCellIndex(cellInfo))

    def testSkyPolygons(self):
        """Test retrieving sky polygons."""
        tractInfo = self.skymap[0]
        patchInfo = tractInfo[0]
        cellInfo = patchInfo[0]

        wcs = tractInfo.wcs

        innerBBox = cellInfo.inner_bbox
        points = wcs.pixelToSky(lsst.geom.Box2D(innerBBox).getCorners())
        poly = lsst.sphgeom.ConvexPolygon([p.getVector() for p in points])

        self.assertEqual(cellInfo.inner_sky_polygon, poly)


class MemoryTester(lsst.utils.tests.MemoryTestCase):
    pass


def setup_module(module):
    lsst.utils.tests.init()


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
