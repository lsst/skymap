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

__all__ = ("SkyMapDataIdPacker",)

from lsst.daf.butler import DataIdPacker, DataId


class SkyMapDataIdPacker(DataIdPacker):
    """A `DataIdPacker` for Tract, Patch and optionally AbstractFilter, given
    a SkyMap.

    Parameters
    ----------
    dimensions : `DataIdPackerDimensions`
        Struct containing dimensions related to this `DataIdPacker`.  Must
        have SkyMap as the only given dimension, Tract, Patch, and possibly
        AbstractFilter as the covered dimensions, and all of these as required
        dimensions.
    skymap : `str`
        SkyMap name from `Registry`.
    tractMax : `int`
        Maximum (exclusive) tract index for this skymap.
    patchNxMax : `int`
        Maximum (exclusive) patch index in the x direction.
    patchNyMax : `int`
        Maximum (exclusive) patch index in the y direction.
    """

    SUPPORTED_FILTERS = [None] + list("ugrizyUBGVRIZYJHK")  # split string into single chars
    """AbstractFilter names supported by this packer.

    New filters should be added to the end of the list to maximize
    compatibility with existing IDs.
    """

    @classmethod
    def getIntFromFilter(cls, name):
        """Return an integer that represents the AbstractFilter with the given
        name.
        """
        try:
            return cls.SUPPORTED_FILTERS.index(name)
        except ValueError:
            raise NotImplementedError(f"AbstractFilter '{name}'' not supported by this ID packer.")

    @classmethod
    def getFilterNameFromInt(cls, num):
        """Return an AbstractFilter name from its integer representation.
        """
        return cls.SUPPORTED_FILTERS[num]

    @classmethod
    def getMaxIntForFilters(cls):
        return len(cls.SUPPORTED_FILTERS)

    @classmethod
    def configure(cls, dimensions):
        # Docstring inherited from DataIdPacker.configure
        assert dimensions.given == ["SkyMap"]
        assert dimensions.required.issuperset(["Tract", "Patch"])
        metadata = {"SkyMap": ["tract_max", "patch_nx_max", "patch_ny_max"]}
        kwds = {}
        return metadata, kwds

    def __init__(self, dimensions, skymap, tractMax, patchNxMax, patchNyMax):
        self._skyMapName = skymap
        self._patchMax = patchNxMax*patchNyMax
        self._tractPatchMax = self._patchMax*tractMax
        if "AbstractFilter" in dimensions.required:
            self._filterMax = self.getMaxIntForFilters()
        else:
            self._filterMax = None

    @property
    def maxBits(self):
        # Docstring inherited from DataIdPacker.maxBits
        packedMax = self._tractPatchMax
        if self._filterMax is not None:
            packedMax *= self._filterMax
        return packedMax.bit_length()

    def _pack(self, dataId):
        # Docstring inherited from DataIdPacker.pack
        packed = dataId["patch"] + self._patchMax*dataId["tract"]
        if self._filterMax is not None:
            packed += self.getIntFromFilter(dataId["abstract_filter"])*self._tractPatchMax
        return packed

    def unpack(self, packedId):
        # Docstring inherited from DataIdPacker.unpack
        d = {"skymap": self._skyMapName}
        if self._filterMax is not None:
            d["abstract_filter"] = self.getFilterNameFromInt(packedId // self._tractPatchMax)
            packedId %= self._tractPatchMax
        d["tract"] = packedId // self._patchMax
        d["patch"] = packedId % self._patchMax
        return DataId(d, dimensions=self.dimensions.required)
