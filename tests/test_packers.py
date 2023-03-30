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
    from lsst.daf.butler import DimensionUniverse, DimensionGraph, DataCoordinate
    HAVE_DAF_BUTLER = True
except ImportError:
    HAVE_DAF_BUTLER = False

from lsst.skymap.packers import SkyMapDimensionPacker


@unittest.skipUnless(HAVE_DAF_BUTLER, "daf_butler not setup")
class SkyMapDimensionPackerTestCase(lsst.utils.tests.TestCase):

    def setUp(self):
        self.universe = DimensionUniverse()
        self.fixed = DataCoordinate.fromFullValues(
            DimensionGraph(universe=self.universe, names=["skymap"]),
            values=("unimportant",),
        ).expanded(
            records={
                "skymap": self.universe["skymap"].RecordClass(
                    name="unimportant",
                    tract_max=5,
                    patch_nx_max=3,
                    patch_ny_max=3,
                )
            }
        )

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
        dimensions = DimensionGraph(universe=self.universe, names=["tract", "patch", "band"])
        dataId = DataCoordinate.standardize(
            skymap=self.fixed["skymap"],
            tract=2,
            patch=6,
            band="g",
            universe=self.universe
        )
        packer = SkyMapDimensionPacker(self.fixed, dimensions)
        packedId = packer.pack(dataId)
        self.assertLessEqual(packedId.bit_length(), packer.maxBits)
        self.assertEqual(packer.unpack(packedId), dataId)

    def test_bad_dimensions(self):
        with self.assertRaises(ValueError):
            SkyMapDimensionPacker(
                self.fixed,
                DimensionGraph(universe=self.universe, names=["tract", "patch", "visit"]),
            )
        with self.assertRaises(ValueError):
            SkyMapDimensionPacker(
                self.fixed,
                DimensionGraph(universe=self.universe, names=["tract", "patch", "detector"]),
            )

    def test_from_config(self):
        data_id = DataCoordinate.standardize(
            skymap=self.fixed["skymap"],
            tract=2,
            patch=6,
            band="g",
            universe=self.universe
        )
        config = SkyMapDimensionPacker.ConfigClass()
        config.n_tracts = 5
        config.n_patches = 9
        config.n_bands = 3
        config.bands = {"r": 0, "g": 1}
        packer = SkyMapDimensionPacker.from_config(data_id, config=config)
        packed_id = packer.pack(data_id)
        self.assertLessEqual(packed_id.bit_length(), packer.maxBits)
        self.assertEqual(packer.unpack(packed_id), data_id)

    def test_from_config_no_bands(self):
        data_id = DataCoordinate.standardize(
            skymap=self.fixed["skymap"],
            tract=2,
            patch=6,
            universe=self.universe
        )
        config = SkyMapDimensionPacker.ConfigClass()
        config.n_tracts = 5
        config.n_patches = 9
        packer = SkyMapDimensionPacker.from_config(data_id, config=config)
        packed_id = packer.pack(data_id)
        self.assertLessEqual(packed_id.bit_length(), packer.maxBits)
        self.assertEqual(packer.unpack(packed_id), data_id)

    def test_from_config_default_bands(self):
        data_id = DataCoordinate.standardize(
            skymap=self.fixed["skymap"],
            tract=2,
            patch=6,
            band="g",
            universe=self.universe
        )
        config = SkyMapDimensionPacker.ConfigClass()
        config.n_tracts = 5
        config.n_patches = 9
        config.n_bands = None
        packer = SkyMapDimensionPacker.from_config(data_id, config=config)
        packed_id = packer.pack(data_id)
        self.assertLessEqual(packed_id.bit_length(), packer.maxBits)
        self.assertEqual(packer.unpack(packed_id), data_id)


class MemoryTester(lsst.utils.tests.MemoryTestCase):
    pass


def setup_module(module):
    lsst.utils.tests.init()


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
