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
    from lsst.daf.butler import ExpandedDataCoordinate, DimensionUniverse, DimensionGraph, DataCoordinate
    HAVE_DAF_BUTLER = True
except ImportError:
    HAVE_DAF_BUTLER = False

from lsst.skymap.packers import SkyMapDimensionPacker


@unittest.skipUnless(HAVE_DAF_BUTLER, "daf_butler not setup")
class SkyMapDimensionPackerTestCase(lsst.utils.tests.TestCase):

    def setUp(self):
        self.universe = DimensionUniverse()
        self.fixed = ExpandedDataCoordinate(
            DimensionGraph(universe=self.universe, names=["skymap"]),
            values=("unimportant",),
            records={
                "skymap": self.universe["skymap"].RecordClass.fromDict({
                    "name": "unimportant",
                    "tract_max": 5,
                    "patch_nx_max": 3,
                    "patch_ny_max": 3,
                })
            })

    def testWithoutFilter(self):
        dimensions = DimensionGraph(universe=self.universe, names=["tract", "patch"])
        dataId = DataCoordinate.standardize(
            skymap=self.fixed["skymap"],
            tract=2,
            patch=6,
            universe=self.universe
        )
        packer = SkyMapDimensionPacker(self.fixed, dimensions)
        packedId = packer.pack(dataId)
        self.assertLessEqual(packedId.bit_length(), packer.maxBits)
        self.assertEqual(packer.unpack(packedId), dataId)

    def testWithFilter(self):
        dimensions = DimensionGraph(universe=self.universe, names=["tract", "patch", "abstract_filter"])
        dataId = DataCoordinate.standardize(
            skymap=self.fixed["skymap"],
            tract=2,
            patch=6,
            abstract_filter="g",
            universe=self.universe
        )
        packer = SkyMapDimensionPacker(self.fixed, dimensions)
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
