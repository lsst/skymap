# LSST Data Management System
# Copyright 2008, 2009, 2010 LSST Corporation.
#
# This product includes software developed by the
# LSST Project (http://www.lsst.org/).
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
# You should have received a copy of the LSST License Statement and
# the GNU General Public License along with this program.  If not,
# see <http://www.lsstcorp.org/LegalNotices/>.
#
import abc
import numbers

import lsst.pex.config as pexConfig
import lsst.pex.exceptions
import lsst.geom as geom
from .patchInfo import PatchInfo

__all__ = ["patchBuilderRegistry", "BasePatchBuilder",
           "OldPatchBuilder", "CellPatchBuilder"]


class BasePatchBuilderConfig(pexConfig.Config):
    """Configuration that is to be shared amongst all patch builders."""
    pass


class BasePatchBuilder(metaclass=abc.ABCMeta):
    """Base class for patch builders.

    Parameters
    ----------
    config : `lsst.pexConfig.Config`
        Input for configuring the algorithm
    """
    def __init__(self, config):
        self.config = config

    @abc.abstractmethod
    def setupPatches(self, minBBox, wcs):
        """Set up the patches of a particular size.

        We grow the bounding box to hold an exact multiple of
        the desired size (patchInnerDimensions), while keeping
        the center roughly the same.  We return the final
        bounding box, and the number of patches in each dimension
        (as an Extent2I).

        Parameters
        ----------
        minBBox : `lsst.geom.Box2I`
            Minimum bounding box for tract
        wcs : `lsst.afw.geom.SkyWcs`
            Wcs object

        Returns
        -------
        bbox : `lsst.geom.Box2I
            final bounding box, number of patches
        numPatches : `lsst.geom.Extent2I`
        """
        raise NotImplementedError("Must be implemented by a subclass.")

    @abc.abstractmethod
    def findPatch(self, tractId, wcs, coord):
        """Find the patch containing the specified coord.

        Parameters
        ----------
        tractId : `int`
            Tract ID number (used for error logging).
        wcs : `lsst.afw.image.SkyWcs`
            The WCS of the associated tract.
        coord : `lsst.geom.SpherePoint`
            ICRS sky coordinate to search for.

        Returns
        -------
        result : `lsst.skymap.PatchInfo`
            PatchInfo of patch whose inner bbox contains the specified coord

        Raises
        ------
        LookupError
            If coord is not in tract or we cannot determine the
            pixel coordinate (which likely means the coord is off the tract).
        """
        raise NotImplementedError("Must be implemented by a subclass.")

    @abc.abstractmethod
    def findPatchList(self, wcs, coordList):
        """Find patches containing the specified list of coords.

        Parameters
        ----------
        wcs : `lsst.afw.image.SkyWcs`
            The WCS of the associated tract.
        coordList : `list` of `lsst.geom.SpherePoint`
            ICRS sky coordinates to search for.

        Returns
        -------
        result : `list` of `lsst.skymap.PatchInfo`
            List of PatchInfo for patches that contain, or may contain, the
            specified region. The list will be empty if there is no overlap.

        Notes
        -----
        **Warning:**

        - This may give incorrect answers on regions that are larger than a
          tract.

        - This uses a naive algorithm that may find some patches that do not
          overlap the region (especially if the region is not a rectangle
          aligned along patch x,y).
        """
        raise NotImplementedError("Must be implemented by a subclass.")

    def getPatchBorder(self):
        return self._patchBorder

    @abc.abstractmethod
    def getPatchInfo(self, index):
        """Return information for the specified patch.

        Parameters
        ----------
        index : `tuple` of `int`
            Index of patch, as a pair of ints;
            or a sequential index as returned by getSequentialPatchIndex;
            negative values are not supported.

        Returns
        -------
        result : `lsst.skymap.PatchInfo`
            The patch info for that index.

        Raises
        ------
        IndexError
            If index is out of range.
        """
        raise NotImplementedError("Must be implemented by a subclass")

    def getPatchInnerDimensions(self):
        """Get dimensions of inner region of the patches (all are the same)
        """
        return self._patchInnerDimensions

    def getPatchIndexPair(self, sequentialIndex):
        """Convert sequential index into patch index (x,y) pair.

        Parameters
        ----------
        sequentialIndex : `int`

        Returns
        -------
        x, y : `int`, `int`
        """
        nx, ny = self._numPatches
        x = sequentialIndex % nx
        y = (sequentialIndex - x) // nx
        return (x, y)


class OldPatchBuilderConfig(BasePatchBuilderConfig):
    patchInnerDimensions = pexConfig.ListField(
        doc="dimensions of inner region of patches (x,y pixels)",
        dtype=int,
        length=2,
        default=(4000, 4000),
    )
    patchBorder = pexConfig.Field(
        doc="border between patch inner and outer bbox (pixels)",
        dtype=int,
        default=100,
    )

    def validate(self):
        if len(self.patchInnerDimensions) != 2:
            raise ValueError("patchInnerDimensions must be 2 ints.")


class OldPatchBuilder(BasePatchBuilder):
    ConfigClass = OldPatchBuilderConfig

    def __init__(self, config):
        super().__init__(config)

        self._patchInnerDimensions = geom.Extent2I(*(val
                                                     for val in config.patchInnerDimensions))
        self._patchBorder = config.patchBorder
        self._initialized = False

    def setupPatches(self, minBBox, wcs):
        bbox = geom.Box2I(minBBox)
        bboxMin = bbox.getMin()
        bboxDim = bbox.getDimensions()
        numPatches = geom.Extent2I(0, 0)
        for i, innerDim in enumerate(self._patchInnerDimensions):
            num = (bboxDim[i] + innerDim - 1) // innerDim  # round up
            deltaDim = (innerDim*num) - bboxDim[i]
            if deltaDim > 0:
                bboxDim[i] = innerDim * num
                bboxMin[i] -= deltaDim // 2
            numPatches[i] = num
        bbox = geom.Box2I(bboxMin, bboxDim)
        self._numPatches = numPatches
        # The final tract BBox starts at zero.
        self._tractBBox = geom.Box2I(geom.Point2I(0, 0), bbox.getDimensions())
        self._initialized = True

        return bbox, numPatches

    def getPatchInfo(self, index):
        # This should always be initialized
        assert self._initialized
        if isinstance(index, numbers.Number):
            index = self.getPatchIndexPair(index)
        if (not 0 <= index[0] < self._numPatches[0]) \
                or (not 0 <= index[1] < self._numPatches[1]):
            raise IndexError("Patch index %s is not in range [0-%d, 0-%d]" %
                             (index, self._numPatches[0] - 1, self._numPatches[1] - 1))
        innerMin = geom.Point2I(*[index[i] * self._patchInnerDimensions[i] for i in range(2)])
        innerBBox = geom.Box2I(innerMin, self._patchInnerDimensions)
        if not self._tractBBox.contains(innerBBox):
            raise RuntimeError(
                "Bug: patch index %s valid but inner bbox=%s not contained in tract bbox=%s" %
                (index, innerBBox, self._tractBBox))
        outerBBox = geom.Box2I(innerBBox)
        outerBBox.grow(self.getPatchBorder())
        outerBBox.clip(self._tractBBox)
        return PatchInfo(
            index=index,
            innerBBox=innerBBox,
            outerBBox=outerBBox,
        )

    def findPatch(self, tractId, wcs, coord):
        try:
            pixel = wcs.skyToPixel(coord)
        except (lsst.pex.exceptions.DomainError, lsst.pex.exceptions.RuntimeError):
            # Point must be way off the tract
            raise LookupError("Unable to determine pixel position for coordinate %s" % (coord,))
        pixelInd = geom.Point2I(pixel)
        if not self._tractBBox.contains(pixelInd):
            raise LookupError("coord %s is not in tract %s" % (coord, tractId))
        patchInd = tuple(int(pixelInd[i]/self._patchInnerDimensions[i]) for i in range(2))
        return self.getPatchInfo(patchInd)

    def findPatchList(self, wcs, coordList):
        box2D = geom.Box2D()
        for coord in coordList:
            try:
                pixelPos = wcs.skyToPixel(coord)
            except (lsst.pex.exceptions.DomainError, lsst.pex.exceptions.RuntimeError):
                # the point is so far off the tract that its pixel position cannot be computed
                continue
            box2D.include(pixelPos)
        bbox = geom.Box2I(box2D)
        bbox.grow(self.getPatchBorder())
        bbox.clip(self._tractBBox)
        if bbox.isEmpty():
            return ()

        llPatchInd = tuple(int(bbox.getMin()[i]/self._patchInnerDimensions[i]) for i in range(2))
        urPatchInd = tuple(int(bbox.getMax()[i]/self._patchInnerDimensions[i]) for i in range(2))
        return tuple(self.getPatchInfo((xInd, yInd))
                     for xInd in range(llPatchInd[0], urPatchInd[0]+1)
                     for yInd in range(llPatchInd[1], urPatchInd[1]+1))


class CellPatchBuilderConfig(BasePatchBuilderConfig):
    cellInnerDimensions = pexConfig.ListField(
        doc="dimensions of inner region of cells (x,y pixels)",
        dtype=int,
        length=2,
        default=(150, 150),
    )
    cellBorder = pexConfig.Field(
        doc="Border between cell inner and outer bbox (pixels)",
        dtype=int,
        default=50,
    )
    numCellsPerPatchInner = pexConfig.Field(
        doc=("Number of cells per inner patch.  There will be a buffer of "
             "one cell all the way around the boundary."),
        dtype=int,
        default=20,
    )

    def validate(self):
        if len(self.cellInnerDimensions) != 2:
            raise ValueError("cellInnerDimensions must be 2 ints.")

        if self.cellInnerDimensions[0] != self.cellInnerDimensions[1]:
            raise ValueError("cellInnerDimensions must be equal (for square cells).")


class CellPatchBuilder(BasePatchBuilder):
    ConfigClass = CellPatchBuilderConfig

    def __init__(self, config):
        super().__init__(config)

        self._cellInnerDimensions = geom.Extent2I(*(val
                                                    for val in config.cellInnerDimensions))
        self._cellBorder = config.cellBorder
        self._numCellsPerPatchInner = config.numCellsPerPatchInner
        self._patchInnerDimensions = geom.Extent2I(*(val*self._numCellsPerPatchInner
                                                     for val in config.cellInnerDimensions))
        self._patchBorder = config.cellInnerDimensions[0]
        self._initialized = False

    def setupPatches(self, minBBox, wcs):
        # This, so far, is the same as the Old version.
        bbox = geom.Box2I(minBBox)
        bboxMin = bbox.getMin()
        bboxDim = bbox.getDimensions()
        numPatches = geom.Extent2I(0, 0)
        for i, innerDim in enumerate(self._patchInnerDimensions):
            num = (bboxDim[i] + innerDim - 1) // innerDim  # round up
            deltaDim = (innerDim*num) - bboxDim[i]
            if deltaDim > 0:
                bboxDim[i] = innerDim * num
                bboxMin[i] -= deltaDim // 2
            numPatches[i] = num
        bbox = geom.Box2I(bboxMin, bboxDim)
        self._numPatches = numPatches
        # The final tract BBox starts at zero.
        self._tractBBox = geom.Box2I(geom.Point2I(0, 0), bbox.getDimensions())
        self._initialized = True

        return bbox, numPatches

    def getPatchInfo(self, index):
        # This should always be initialized
        assert self._initialized
        if isinstance(index, numbers.Number):
            index = self.getPatchIndexPair(index)
        if (not 0 <= index[0] < self._numPatches[0]) \
           or (not 0 <= index[1] < self._numPatches[1]):
            raise IndexError("Patch index %s is not in range [0-%d, 0-%d]" %
                             (index, self._numPatches[0] - 1, self._numPatches[1] - 1))
        innerMin = geom.Point2I(*[index[i] * self._patchInnerDimensions[i] for i in range(2)])
        innerBBox = geom.Box2I(innerMin, self._patchInnerDimensions)
        if not self._tractBBox.contains(innerBBox):
            raise RuntimeError(
                "Bug: patch index %s valid but inner bbox=%s not contained in tract bbox=%s" %
                (index, innerBBox, self._tractBBox))
        outerBBox = geom.Box2I(innerBBox)
        outerBBox.grow(self.getPatchBorder())
        # We do not clip the patch
        return PatchInfo(
            index=index,
            innerBBox=innerBBox,
            outerBBox=outerBBox,
        )

    def findPatch(self, tractId, wcs, coord):
        # This is the same as above
        try:
            pixel = wcs.skyToPixel(coord)
        except (lsst.pex.exceptions.DomainError, lsst.pex.exceptions.RuntimeError):
            # Point must be way off the tract
            raise LookupError("Unable to determine pixel position for coordinate %s" % (coord,))
        pixelInd = geom.Point2I(pixel)
        if not self._tractBBox.contains(pixelInd):
            raise LookupError("coord %s is not in tract %s" % (coord, tractId))
        # This should probably be the same as above because we only
        # care about the INNER dimensions.
        patchInd = tuple(int(pixelInd[i]/self._patchInnerDimensions[i]) for i in range(2))
        return self.getPatchInfo(patchInd)

    def findPatchList(self, wcs, coordList):
        # Same as above
        box2D = geom.Box2D()
        for coord in coordList:
            try:
                pixelPos = wcs.skyToPixel(coord)
            except (lsst.pex.exceptions.DomainError, lsst.pex.exceptions.RuntimeError):
                # the point is so far off the tract that its pixel position cannot be computed
                continue
            box2D.include(pixelPos)
        bbox = geom.Box2I(box2D)
        bbox.grow(self.getPatchBorder())
        bbox.clip(self._tractBBox)
        if bbox.isEmpty():
            return ()

        llPatchInd = tuple(int(bbox.getMin()[i]/self._patchInnerDimensions[i]) for i in range(2))
        urPatchInd = tuple(int(bbox.getMax()[i]/self._patchInnerDimensions[i]) for i in range(2))
        return tuple(self.getPatchInfo((xInd, yInd))
                     for xInd in range(llPatchInd[0], urPatchInd[0]+1)
                     for yInd in range(llPatchInd[1], urPatchInd[1]+1))


patchBuilderRegistry = pexConfig.makeRegistry(
    doc="A registry of Patch Builders (subclasses of PatchBuilder)",
)

patchBuilderRegistry.register("old", OldPatchBuilder)
patchBuilderRegistry.register("cells", CellPatchBuilder)
