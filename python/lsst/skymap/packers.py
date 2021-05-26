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

__all__ = ("SkyMapDimensionPacker",)

from lsst.daf.butler import DimensionPacker, DimensionGraph, DataCoordinate


class SkyMapDimensionPacker(DimensionPacker):
    """A `DimensionPacker` for tract, patch and optionally band,
    given a SkyMap.

    Parameters
    ----------
    fixed : `lsst.daf.butler.DataCoordinate`
        Expanded data ID that must include at least the skymap dimension.
    dimensions : `lsst.daf.butler.DimensionGraph`
        The dimensions of data IDs packed by this instance.  Must include
        skymap, tract, and patch, and may include band.
    """

    SUPPORTED_FILTERS = (
        [None]
        + list("ugrizyUBGVRIZYJHK")  # split string into single chars
        + [f"N{d}" for d in (387, 515, 656, 816, 1010)]  # HSC narrow-bands
    )
    """band names supported by this packer.

    New filters should be added to the end of the list to maximize
    compatibility with existing IDs.
    """

    @classmethod
    def getIntFromFilter(cls, name):
        """Return an integer that represents the band with the given
        name.
        """
        try:
            return cls.SUPPORTED_FILTERS.index(name)
        except ValueError:
            raise NotImplementedError(f"band '{name}' not supported by this ID packer.")

    @classmethod
    def getFilterNameFromInt(cls, num):
        """Return an band name from its integer representation.
        """
        return cls.SUPPORTED_FILTERS[num]

    @classmethod
    def getMaxIntForFilters(cls):
        return len(cls.SUPPORTED_FILTERS)

    @classmethod
    def configure(cls, dimensions):
        # Docstring inherited from DataIdPacker.configure
        assert dimensions.given == ["skymap"]
        assert dimensions.required.issuperset(["tract", "patch"])
        metadata = {"skymap": ["tract_max", "patch_nx_max", "patch_ny_max"]}
        kwds = {}
        return metadata, kwds

    def __init__(self, fixed: DataCoordinate, dimensions: DimensionGraph):
        super().__init__(fixed, dimensions)
        record = fixed.records["skymap"]
        self._skyMapName = record.name
        self._patchMax = record.patch_nx_max * record.patch_ny_max
        self._tractPatchMax = self._patchMax*record.tract_max
        if "band" in dimensions:
            self._filterMax = self.getMaxIntForFilters()
        else:
            self._filterMax = None

    @property
    def maxBits(self) -> int:
        # Docstring inherited from DataIdPacker.maxBits
        packedMax = self._tractPatchMax
        if self._filterMax is not None:
            packedMax *= self._filterMax
        return packedMax.bit_length()

    def _pack(self, dataId: DataCoordinate) -> int:
        # Docstring inherited from DataIdPacker.pack
        packed = dataId["patch"] + self._patchMax*dataId["tract"]
        if self._filterMax is not None:
            packed += self.getIntFromFilter(dataId["band"])*self._tractPatchMax
        return packed

    def unpack(self, packedId: int) -> DataCoordinate:
        # Docstring inherited from DataIdPacker.unpack
        d = {"skymap": self._skyMapName}
        if self._filterMax is not None:
            d["band"] = self.getFilterNameFromInt(packedId // self._tractPatchMax)
            packedId %= self._tractPatchMax
        d["tract"] = packedId // self._patchMax
        d["patch"] = packedId % self._patchMax
        return DataCoordinate.standardize(d, graph=self.dimensions)
