#
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

__all__ = ["TractInfo"]

import numpy as np
from deprecated.sphinx import deprecated

import lsst.pex.exceptions
import lsst.geom as geom
from lsst.sphgeom import ConvexPolygon, Box

from .detail import makeSkyPolygonFromBBox, Index2D


class TractInfo:
    """Information about a tract in a SkyMap sky pixelization

    Parameters
    ----------
    id : `int`
        tract ID
    tractBuilder : Subclass of `lsst.skymap.BaseTractBuilder`
        Object used to compute patch geometry.
    ctrCoord : `lsst.geom.SpherePoint`
        ICRS sky coordinate of center of inner region of tract; also used as
        the CRVAL for the WCS.
    vertexCoordList : `list` of `lsst.geom.SpherePoint`
        Vertices that define the boundaries of the inner region.
    tractOverlap : `lsst.geom.Angle`
        Minimum overlap between adjacent sky tracts; this defines the minimum
        distance the tract extends beyond the inner region in all directions.
    wcs : `lsst.afw.image.SkyWcs`
        WCS for tract. The reference pixel will be shifted as required so that
        the lower left-hand pixel (index 0,0) has pixel position 0.0, 0.0.
    innerBoxCorners : `list` [`lsst.sphgeom.LonLat`], optional
        If set then the ``inner_sky_region`` will be a `lsst.sphgeom.Box` with
        these corners as opposed to a `lsst.sphgeom.ConvexPolygon` built from
        the ``vertex_list``.

    Notes
    -----
    The tract is subdivided into rectangular patches. Each patch has the
    following properties:

    - An inner region defined by an inner bounding box. The inner regions of
      the patches exactly tile the tract, and all inner regions have the same
      dimensions. The tract is made larger as required to make this work.

    - An outer region defined by an outer bounding box. The outer region
      extends beyond the inner region by patchBorder pixels in all directions,
      except there is no border at the edges of the tract.
      Thus patches overlap each other but never extend off the tract.
      If you do not want any overlap between adjacent patches then set
      patchBorder to 0.

    - An index that consists of a pair of integers:

      * 0 <= x index < numPatches[0]

      * 0 <= y index < numPatches[1]

      Patch 0,0 is at the minimum corner of the tract bounding box.

    - It is not enforced that ctrCoord is the center of vertexCoordList, but
      SkyMap relies on it.
    """
    def __init__(self, id, tractBuilder, ctrCoord, vertexCoordList, tractOverlap, wcs, innerBoxCorners=None):
        self._id = id
        self._ctrCoord = ctrCoord
        self._vertexCoordList = tuple(vertexCoordList)
        self._tractOverlap = tractOverlap
        self._tractBuilder = tractBuilder
        self._innerBoxCorners = innerBoxCorners

        minBBox = self._minimumBoundingBox(wcs)
        initialBBox, self._numPatches = self._tractBuilder.setupPatches(minBBox, wcs)
        self._bbox, self._wcs = self._finalOrientation(initialBBox, wcs)

    def _minimumBoundingBox(self, wcs):
        """Calculate the minimum bounding box for the tract, given the WCS.

        The bounding box is created in the frame of the supplied WCS,
        so that it's OK if the coordinates are negative.

        We compute the bounding box that holds all the vertices and the
        desired overlap.
        """
        minBBoxD = geom.Box2D()
        halfOverlap = self._tractOverlap / 2.0
        for vertexCoord in self._vertexCoordList:
            if self._tractOverlap == 0:
                minBBoxD.include(wcs.skyToPixel(vertexCoord))
            else:
                numAngles = 24
                angleIncr = geom.Angle(360.0, geom.degrees) / float(numAngles)
                for i in range(numAngles):
                    offAngle = angleIncr * i
                    offCoord = vertexCoord.offset(offAngle, halfOverlap)
                    pixPos = wcs.skyToPixel(offCoord)
                    minBBoxD.include(pixPos)
        return minBBoxD

    def _finalOrientation(self, bbox, wcs):
        """Determine the final orientation

        We offset everything so the lower-left corner is at 0,0
        and compute the final Wcs.

        Parameters
        ----------
        bbox : `lsst.geom.Box2I`
            Current bounding box.
        wcs : `lsst.afw.geom.SkyWcs
            Current Wcs.

        Returns
        -------
        finalBBox : `lsst.geom.Box2I`
            Revised bounding box.
        wcs : `lsst.afw.geom.SkyWcs`
            Revised Wcs.
        """
        finalBBox = geom.Box2I(geom.Point2I(0, 0), bbox.getDimensions())
        # shift the WCS by the same amount as the bbox; extra code is required
        # because simply subtracting makes an Extent2I
        pixPosOffset = geom.Extent2D(finalBBox.getMinX() - bbox.getMinX(),
                                     finalBBox.getMinY() - bbox.getMinY())
        wcs = wcs.copyAtShiftedPixelOrigin(pixPosOffset)
        return finalBBox, wcs

    def getSequentialPatchIndex(self, patchInfo):
        """Return a single integer that uniquely identifies
        the given patch within this tract.

        Parameters
        ----------
        patchInfo : `lsst.skymap.PatchInfo`

        Returns
        -------
        sequentialIndex : `int`
        """
        return self._tractBuilder.getSequentialPatchIndex(patchInfo)

    def getSequentialPatchIndexFromPair(self, index):
        """Return a single integer that uniquely identifies
        the patch index within the tract.

        Parameters
        ----------
        index : `lsst.skymap.Index2D`

        Returns
        -------
        sequentialIndex : `int`
        """
        return self._tractBuilder.getSequentialPatchIndexFromPair(index)

    def getPatchIndexPair(self, sequentialIndex):
        """Convert sequential index into patch index (x,y) pair.

        Parameters
        ----------
        sequentialIndex : `int`

        Returns
        -------
        x, y : `lsst.skymap.Index2D`
        """
        return self._tractBuilder.getPatchIndexPair(sequentialIndex)

    def findPatch(self, coord):
        """Find the patch containing the specified coord.

        Parameters
        ----------
        coord : `lsst.geom.SpherePoint`
            ICRS sky coordinate to search for.

        Returns
        -------
        result : `lsst.skymap.PatchInfo`
            PatchInfo of patch whose inner bbox contains the specified coord

        Raises
        ------
        LookupError
            Raised if coord is not in tract or we cannot determine the
            pixel coordinate (which likely means the coord is off the tract).
        """
        try:
            pixel = self.wcs.skyToPixel(coord)
        except (lsst.pex.exceptions.DomainError, lsst.pex.exceptions.RuntimeError):
            # Point must be way off the tract
            raise LookupError("Unable to determine pixel position for coordinate %s" % (coord,))
        pixelInd = geom.Point2I(pixel)
        if not self._bbox.contains(pixelInd):
            raise LookupError("coord %s is not in tract %s" % (coord, self.tract_id))
        # This should probably be the same as above because we only
        # care about the INNER dimensions.
        patchInd = tuple(int(pixelInd[i]/self.patch_inner_dimensions[i]) for i in range(2))
        return self.getPatchInfo(patchInd)

    def findPatchList(self, coordList):
        """Find patches containing the specified list of coords.

        Parameters
        ----------
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
        box2D = geom.Box2D()
        for coord in coordList:
            try:
                pixelPos = self.wcs.skyToPixel(coord)
            except (lsst.pex.exceptions.DomainError, lsst.pex.exceptions.RuntimeError):
                # The point is so far off the tract that its pixel position
                # cannot be computed.
                continue
            box2D.include(pixelPos)
        bbox = geom.Box2I(box2D)
        bbox.grow(self.getPatchBorder())
        bbox.clip(self._bbox)
        if bbox.isEmpty():
            return ()

        llPatchInd = tuple(int(bbox.getMin()[i]/self.patch_inner_dimensions[i]) for i in range(2))
        urPatchInd = tuple(int(bbox.getMax()[i]/self.patch_inner_dimensions[i]) for i in range(2))
        return tuple(self.getPatchInfo((xInd, yInd))
                     for xInd in range(llPatchInd[0], urPatchInd[0]+1)
                     for yInd in range(llPatchInd[1], urPatchInd[1]+1))

    def getBBox(self):
        """Get bounding box of tract (as an geom.Box2I)
        """
        return geom.Box2I(self._bbox)

    bbox = property(getBBox)

    def getCtrCoord(self):
        """Get ICRS sky coordinate of center of tract
        (as an lsst.geom.SpherePoint)
        """
        return self._ctrCoord

    ctr_coord = property(getCtrCoord)

    def getId(self):
        """Get ID of tract
        """
        return self._id

    tract_id = property(getId)

    def getNumPatches(self):
        """Get the number of patches in x, y.

        Returns
        -------
        result : `lsst.skymap.Index2D`
            The number of patches in x, y
        """
        return self._numPatches

    num_patches = property(getNumPatches)

    def getPatchBorder(self):
        return self._tractBuilder.getPatchBorder()

    patch_border = property(getPatchBorder)

    def getPatchInfo(self, index):
        """Return information for the specified patch.

        Parameters
        ----------
        index : `typing.NamedTuple` ['x': `int`, 'y': `int`]
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
            Raised if index is out of range.
        """
        return self._tractBuilder.getPatchInfo(index, self._wcs)

    def getPatchInnerDimensions(self):
        """Get dimensions of inner region of the patches (all are the same)
        """
        return self._tractBuilder.getPatchInnerDimensions()

    patch_inner_dimensions = property(getPatchInnerDimensions)

    def getTractOverlap(self):
        """Get minimum overlap of adjacent sky tracts.
        """
        return self._tractOverlap

    tract_overlap = property(getTractOverlap)

    def getVertexList(self):
        """Get list of ICRS sky coordinates of vertices that define the
        boundary of the inner region.

        Notes
        -----
        **warning:** this is not a deep copy.
        """
        return self._vertexCoordList

    vertex_list = property(getVertexList)

    # TODO: Remove with DM-44799
    @deprecated(reason="getInnerSkyPolygon()/inner_sky_polygon has been deprecated in favor of "
                       "inner_sky_region, and will be removed after v28.",
                category=FutureWarning, version=28)
    def getInnerSkyPolygon(self):
        """Get inner on-sky region as a sphgeom.ConvexPolygon.
        """
        skyUnitVectors = [sp.getVector() for sp in self.getVertexList()]
        return ConvexPolygon.convexHull(skyUnitVectors)

    inner_sky_polygon = property(getInnerSkyPolygon)

    def getInnerSkyRegion(self):
        """Get inner on-sky region.
        """
        if self._innerBoxCorners:
            return Box(point1=self._innerBoxCorners[0], point2=self._innerBoxCorners[1])
        else:
            skyUnitVectors = [sp.getVector() for sp in self.getVertexList()]
            return ConvexPolygon.convexHull(skyUnitVectors)

    inner_sky_region = property(getInnerSkyRegion)

    def getOuterSkyPolygon(self):
        """Get outer on-sky region as a sphgeom.ConvexPolygon
        """
        return makeSkyPolygonFromBBox(bbox=self.getBBox(), wcs=self.getWcs())

    outer_sky_polygon = property(getOuterSkyPolygon)

    def getWcs(self):
        """Get WCS of tract.

        Returns
        -------
        wcs : `lsst.afw.geom.SkyWcs`
            The WCS of this tract
        """
        return self._wcs

    wcs = property(getWcs)

    def __str__(self):
        return "TractInfo(id=%s)" % (self._id,)

    def __repr__(self):
        return "TractInfo(id=%s, ctrCoord=%s)" % (self._id, self._ctrCoord.getVector())

    def __iter__(self):
        xNum, yNum = self.getNumPatches()
        for y in range(yNum):
            for x in range(xNum):
                yield self.getPatchInfo(Index2D(x=x, y=y))

    def __len__(self):
        xNum, yNum = self.getNumPatches()
        return xNum*yNum

    def __getitem__(self, index):
        return self.getPatchInfo(index)

    def contains(self, coord):
        """Does this tract contain the coordinate?"""
        try:
            pixel = self.getWcs().skyToPixel(coord)
        except (lsst.pex.exceptions.DomainError, lsst.pex.exceptions.RuntimeError):
            # Point must be way off the tract
            return False
        if not np.isfinite(pixel.getX()) or not np.isfinite(pixel.getY()):
            # Point is definitely off the tract
            return False
        return self.getBBox().contains(geom.Point2I(pixel))


class ExplicitTractInfo(TractInfo):
    """Information for a tract specified explicitly.

    A tract is placed at the explicitly defined coordinates, with the nominated
    radius.  The tracts are square (i.e., the radius is really a half-size).

    Parameters
    ----------
    id : : `int`
        tract ID
    tractBuilder : Subclass of `lsst.skymap.BaseTractBuilder`
        Object used to compute patch geometry.
    ctrCoord : `lsst.geom.SpherePoint`
        ICRS sky coordinate of center of inner region of tract; also used as
        the CRVAL for the WCS.
    radius : `lsst.geom.Angle`
        Radius of the tract.
    tractOverlap : `lsst.geom.Angle`
        Minimum overlap between adjacent sky tracts; this defines the minimum
        distance the tract extends beyond the inner region in all directions.
    wcs : `lsst.afw.image.SkyWcs`
        WCS for tract. The reference pixel will be shifted as required so that
        the lower left-hand pixel (index 0,0) has pixel position 0.0, 0.0.
    innerBoxCorners : `list` [`lsst.sphgeom.LonLat`], optional
        If set then the ``inner_sky_region`` will be a `lsst.sphgeom.Box` with
        these corners as oppsed to a `lsst.sphgeom.ConvexPolygon` built from
        the ``vertex_list``.
    """
    def __init__(self, id, tractBuilder, ctrCoord, radius, tractOverlap, wcs, innerBoxCorners=None):
        # We don't want TractInfo setting the bbox on the basis of vertices,
        # but on the radius.
        vertexList = []
        self._radius = radius
        super(ExplicitTractInfo, self).__init__(
            id,
            tractBuilder,
            ctrCoord,
            vertexList,
            tractOverlap,
            wcs,
            innerBoxCorners=innerBoxCorners,
        )
        # Shrink the box slightly to make sure the vertices are in the tract
        bboxD = geom.BoxD(self.getBBox())
        bboxD.grow(-0.001)
        finalWcs = self.getWcs()
        self._vertexCoordList = finalWcs.pixelToSky(bboxD.getCorners())

    def _minimumBoundingBox(self, wcs):
        """Calculate the minimum bounding box for the tract, given the WCS, and
        the nominated radius.
        """
        bbox = geom.Box2D()
        for i in range(4):
            cornerCoord = self._ctrCoord.offset(i*90*geom.degrees, self._radius + self._tractOverlap)
            pixPos = wcs.skyToPixel(cornerCoord)
            bbox.include(pixPos)
        return bbox
