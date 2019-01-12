# This file is part of skymap.
#
# Developed for the LSST Data Management System.
# This product includes software developed by the LSST Project
# (http://www.lsst.org).
# See the COPYRIGHT file at the top-level directory of this distribution
# for details of code ownership.
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

import unittest

import lsst.utils.tests

try:
    from lsst.daf.butler import DataId, DataIdPackerDimensions, DimensionGraph, DimensionSet
    HAVE_DAF_BUTLER = True
except ImportError:
    HAVE_DAF_BUTLER = False

from lsst.skymap.packers import SkyMapDataIdPacker


@unittest.skipUnless(HAVE_DAF_BUTLER, "daf_butler not setup")
class SkyMapDataIdPackerTestCase(lsst.utils.tests.TestCase):

    def setUp(self):
        self.universe = DimensionGraph.fromConfig()
        self.given = DimensionSet(universe=self.universe, elements=["SkyMap"])
        self.parameters = dict(skymap="unimportant", tractMax=5, patchNxMax=3, patchNyMax=3)

    def testWithoutFilter(self):
        covered = DimensionSet(universe=self.universe, elements=["Tract", "Patch"])
        dimensions = DataIdPackerDimensions(given=self.given, required=self.given.union(covered))
        dataId = DataId(skymap=self.parameters["skymap"], tract=2, patch=6, universe=self.universe)
        packer = SkyMapDataIdPacker(dimensions, **self.parameters)
        packedId = packer.pack(dataId)
        self.assertLessEqual(packedId.bit_length(), packer.maxBits)
        self.assertEqual(packer.unpack(packedId), dataId)

    def testWithFilter(self):
        covered = DimensionSet(universe=self.universe, elements=["Tract", "Patch", "AbstractFilter"])
        dimensions = DataIdPackerDimensions(given=self.given, required=self.given.union(covered))
        dataId = DataId(skymap=self.parameters["skymap"], tract=2, patch=6, abstract_filter="g",
                        universe=self.universe)
        packer = SkyMapDataIdPacker(dimensions, **self.parameters)
        packedId = packer.pack(dataId)
        self.assertLessEqual(packedId.bit_length(), packer.maxBits)
        self.assertEqual(packer.unpack(packedId), dataId)


class MemoryTester(lsst.utils.tests.MemoryTestCase):
    pass


def setup_module(module):
    lsst.utils.tests.init()


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
