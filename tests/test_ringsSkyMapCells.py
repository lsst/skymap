import unittest
# import math
# import numpy as np

import lsst.utils.tests
import lsst.geom

from lsst.skymap.ringsSkyMap import RingsSkyMap
from helper import skyMapTestCase


class RingsCellsTestCase(skyMapTestCase.SkyMapTestCase):
    def setUp(self):
        config = RingsSkyMap.ConfigClass()
        config.numRings = 3
        config.patchBuilder = "cells"

        self.setAttributes(
            SkyMapClass=RingsSkyMap,
            name="rings",
            config=config,
            numTracts=26,
            neighborAngularSeparation=None,  # no uniform tract separation
            numNeighbors=None,    # ignored because neighborAngularSeparation=None
        )

    def testCellInfo(self):
        """Test retrieving cell info."""
        skymap = self.getSkyMap()

        tractInfo = skymap[13]

        cellConfig = self.config.patchBuilder['cells']

        numCellsPerPatch = cellConfig.numCellsPerPatchInner + 2
        patchInnerSize = (cellConfig.numCellsPerPatchInner*cellConfig.cellInnerDimensions[0],
                          cellConfig.numCellsPerPatchInner*cellConfig.cellInnerDimensions[1])
        patchOuterSize = (numCellsPerPatch*cellConfig.cellInnerDimensions[0],
                          numCellsPerPatch*cellConfig.cellInnerDimensions[1])

        for patchInfo in tractInfo:
            # Check that all the patchInfos have the correct size.
            self.assertEqual(len(patchInfo), numCellsPerPatch**2)
            innerBBox = patchInfo.getInnerBBox()
            outerBBox = patchInfo.getOuterBBox()
            self.assertEqual(innerBBox.getDimensions()[0], patchInnerSize[0])
            self.assertEqual(innerBBox.getDimensions()[1], patchInnerSize[1])
            self.assertEqual(outerBBox.getDimensions()[0], patchOuterSize[0])
            self.assertEqual(outerBBox.getDimensions()[1], patchOuterSize[1])

        # Look at a patch in first corner, middle, and last corner.
        numPatches = tractInfo.getNumPatches()
        names = ['ll', 'center', 'ur']
        patchInfos = [tractInfo[0],
                      tractInfo[(numPatches[0]//2, numPatches[1]//2)],
                      tractInfo[(numPatches[0] - 1, numPatches[1] - 1)]]
        for name, patchInfo in zip(names, patchInfos):
            # Check that the corner of the edge patches are what we expect.
            # (the dimensions of the patches were checked above)
            index = patchInfo.getIndex()
            innerCorner = (index[0]*patchInnerSize[0], index[1]*patchInnerSize[1])

            self.assertEqual(patchInfo.getInnerBBox().getBeginX(), innerCorner[0])
            self.assertEqual(patchInfo.getInnerBBox().getBeginY(), innerCorner[1])

            self.assertEqual(patchInfo.getOuterBBox().getBeginX(),
                             innerCorner[0] - cellConfig.cellInnerDimensions[0])
            self.assertEqual(patchInfo.getOuterBBox().getBeginY(),
                             innerCorner[1] - cellConfig.cellInnerDimensions[1])

            # And confirm that the lower left corner of the lower left patch
            # goes negative.
            if name == 'll':
                self.assertLess(patchInfo.getOuterBBox().getBeginX(), 0)
                self.assertLess(patchInfo.getOuterBBox().getBeginY(), 0)

            patchInnerBBox = patchInfo.getInnerBBox()

            # Check that all the cells are where we expect.
            for cellInfo in patchInfo:
                index = cellInfo.getIndex()

                innerBBox = cellInfo.getInnerBBox()
                outerBBox = cellInfo.getOuterBBox()

                innerMin = ((index[0] - 1)*cellConfig.cellInnerDimensions[0] + patchInnerBBox.getBeginX(),
                            (index[1] - 1)*cellConfig.cellInnerDimensions[1] + patchInnerBBox.getBeginY())
                outerMin = (innerMin[0] - cellConfig.cellBorder,
                            innerMin[1] - cellConfig.cellBorder)
                self.assertEqual(innerBBox.getBeginX(), innerMin[0])
                self.assertEqual(innerBBox.getBeginY(), innerMin[1])
                self.assertEqual(innerBBox.getWidth(), cellConfig.cellInnerDimensions[0])
                self.assertEqual(innerBBox.getHeight(), cellConfig.cellInnerDimensions[1])
                self.assertEqual(outerBBox.getBeginX(), outerMin[0])
                self.assertEqual(outerBBox.getBeginY(), outerMin[1])
                self.assertEqual(outerBBox.getWidth(),
                                 cellConfig.cellInnerDimensions[0] + 2*cellConfig.cellBorder)
                self.assertEqual(outerBBox.getHeight(),
                                 cellConfig.cellInnerDimensions[1] + 2*cellConfig.cellBorder)
                if index[0] == 0:
                    self.assertLess(innerBBox.getBeginX(), patchInnerBBox.getBeginX())
                    self.assertEqual(innerBBox.getEndX(), patchInnerBBox.getBeginX())
                if index[1] == 0:
                    self.assertLess(innerBBox.getBeginY(), patchInnerBBox.getBeginY())
                    self.assertEqual(innerBBox.getEndY(), patchInnerBBox.getBeginY())
                if index[0] == patchInfo.getNumCells()[0] - 1:
                    self.assertEqual(innerBBox.getBeginX(),
                                     patchInfo.getInnerBBox().getEndX())
                    self.assertGreater(innerBBox.getEndX(),
                                       patchInfo.getInnerBBox().getEndX())
                if index[1] == patchInfo.getNumCells()[1] - 1:
                    self.assertEqual(innerBBox.getBeginY(),
                                     patchInfo.getInnerBBox().getEndY())
                    self.assertGreater(innerBBox.getEndY(),
                                       patchInfo.getInnerBBox().getEndY())


class MemoryTester(lsst.utils.tests.MemoryTestCase):
    pass


def setup_module(module):
    lsst.utils.tests.init()


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
