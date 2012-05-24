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
from lsst.pex.exceptions import LsstCppException
import lsst.afw.coord as afwCoord
import lsst.afw.geom as afwGeom
import lsst.afw.image as afwImage
from .patchInfo import PatchInfo

__all__ = ["TractInfo"]

class TractInfo(object):
    """Information about a tract in a SkyMap sky pixelization
    
    The tract is subdivided into rectangular patches. Each patch has the following properties:
    - An inner region defined by an inner bounding. The inner regions of the patches exactly tile the tract,
      and all inner regions have the same dimensions. The tract is made larger as required to make this work.
    - An outer region defined by an outer bounding box. The outer region extends beyond the inner region
      by patchBorder pixels in all directions, except there is no border at the edges of the tract.
      Thus patches overlap each other but never extend off the tract. If you do not want any overlap
      between adjacent patches then set patchBorder to 0.
    - An index that consists of a pair of integers:
        0 <= x index < numPatches[0]
        0 <= y index < numPatches[1]
      Patch 0,0 is at the minimum corner of the tract bounding box.
    """
    def __init__(self, id, patchInnerDimensions, patchBorder, ctrCoord, vertexCoordList, tractOverlap,
        wcsFactory):
        """Construct a TractInfo

        @param[in] id: tract ID
        @param[in] patchInnerDimensions: dimensions of inner region of patches (x,y pixels)
        @param[in] patchBorder: overlap between adjacent patches (in pixels, one int)
        @param[in] ctrCoord: sky coordinate of center of inner region of tract, as an afwCoord.Coord;
            also used as the CRVAL for the WCS.
        @param[in] vertexCoordList: list of sky coordinates (afwCoord.Coord)
            of vertices that define the boundaries of the inner region
        @param[in] tractOverlap: minimum overlap between adjacent sky tracts; an afwGeom.Angle;
            this defines the minimum distance the tract extends beyond the inner region in all directions
        @param[in] wcsFactory: a skymap.detail.WcsFactory object
        
        @warning
        - It is not enforced that ctrCoord is the center of vertexCoordList, but SkyMap relies on it
        - vertexCoordList will likely become a geom SphericalConvexPolygon someday.
        """
        self._id = id
        try:
            assert len(patchInnerDimensions) == 2
            self._patchInnerDimensions = afwGeom.Extent2I(*(int(val) for val in patchInnerDimensions))
        except:
            raise TypeError("patchInnerDimensions=%s; must be two ints" % (patchInnerDimensions,))
        self._patchBorder = int(patchBorder)
        self._ctrCoord = ctrCoord
        self._vertexCoordList = tuple(coord.clone() for coord in vertexCoordList)
        self._tractOverlap = tractOverlap
        
        # We don't know how big the tract will be, yet, so start by computing everything
        # as if the tract center was at pixel position 0, 0; then shift all pixel positions
        # so that the tract's bbox starts from 0,0
        initialCRPixPos = afwGeom.Point2D(0.0, 0.0)
        initialWcs = wcsFactory.makeWcs(crPixPos=initialCRPixPos, crValCoord=self._ctrCoord)

        # compute minimum bounding box that will hold all corners and tractOverlap
        minBBoxD = afwGeom.Box2D()
        halfOverlap = self._tractOverlap / 2.0
        for vertexCoord in self._vertexCoordList:
            vertexDeg = vertexCoord.getPosition(afwGeom.degrees)
            if self._tractOverlap == 0:
                minBBoxD.include(initialWcs.skyToPixel(vertexCoord))
            else:
                numAngles = 24
                angleIncr = afwGeom.Angle(360.0, afwGeom.degrees) / float(numAngles)
                for i in range(numAngles):
                    offAngle = angleIncr * i
                    offCoord = vertexCoord.clone()
                    offCoord.offset(offAngle, halfOverlap)
                    pixPos = initialWcs.skyToPixel(offCoord)
                    minBBoxD.include(pixPos)
        initialBBox = afwGeom.Box2I(minBBoxD)

        # grow initialBBox to hold exactly a multiple of patchInnerDimensions in x,y,
        # while keeping the center roughly the same
        initialBBoxMin = initialBBox.getMin()
        bboxDim = initialBBox.getDimensions()
        self._numPatches = afwGeom.Extent2I(0, 0)
        for i, innerDim in enumerate(self._patchInnerDimensions):
            numPatches =  (bboxDim[i] + innerDim - 1) // innerDim # round up
            deltaDim = (innerDim * numPatches) - bboxDim[i]
            if deltaDim > 0:
                bboxDim[i] = innerDim * numPatches
                initialBBoxMin[i] -= deltaDim // 2
            self._numPatches[i] = numPatches
        initialBBox = afwGeom.Box2I(initialBBoxMin, bboxDim)

        # compute final bbox the same size but with LL corner = 0,0; use that to compute final WCS
        self._bbox = afwGeom.Box2I(afwGeom.Point2I(0, 0), initialBBox.getDimensions())
        # crpix for final Wcs is shifted by the same amount as bbox
        pixPosOffset = afwGeom.Extent2D(self._bbox.getMinX() - initialBBox.getMinX(),
                                        self._bbox.getMinY() - initialBBox.getMinY())
        crPixPos = initialCRPixPos + pixPosOffset
        self._wcs = wcsFactory.makeWcs(crPixPos=crPixPos, crValCoord=self._ctrCoord)

    def findPatch(self, coord):
        """Find the patch containing the specified coord

        @param[in] coord: sky coordinate (afwCoord.Coord)
        @return PatchInfo of patch whose inner bbox contains the specified coord
        
        @raise LookupError if coord is not in tract
        
        @note This routine will be more efficient if coord is ICRS.
        """
        pixelInd = afwGeom.Point2I(self.getWcs().skyToPixel(coord.toIcrs()))
        if not self.getBBox().contains(pixelInd):
            raise LookupError("coord %s is not in tract %s" % (coord, self.getId()))
        patchInd = tuple(int(pixelInd[i]/self._patchInnerDimensions[i]) for i in range(2))
        return self.getPatchInfo(patchInd)
    
    def findPatchList(self, coordList):
        """Find patches containing the specified list of coords
        
        @param[in] coordList: list of sky coordinates (afwCoord.Coord)
        @return list of PatchInfo for patches that contain, or may contain, the specified region.
            The list will be empty if there is no overlap.
        
        @warning:
        * This may give incorrect answers on regions that are larger than a tract
        * This uses a naive algorithm that may find some patches that do not overlap the region
            (especially if the region is not a rectangle aligned along patch x,y).
        """
        box2D = afwGeom.Box2D()
        for coord in coordList:
            try:
                pixelPos = self.getWcs().skyToPixel(coord.toIcrs())
            except LsstCppException:
                # the point is so far off the tract that its pixel position cannot be computed
                continue
            box2D.include(pixelPos)
        bbox = afwGeom.Box2I(box2D)
        bbox.grow(self.getPatchBorder())
        bbox.clip(self.getBBox())
        if bbox.isEmpty():
            return ()

        llPatchInd = tuple(int(bbox.getMin()[i]/self._patchInnerDimensions[i]) for i in range(2))
        urPatchInd = tuple(int(bbox.getMax()[i]/self._patchInnerDimensions[i]) for i in range(2))
        return tuple(self.getPatchInfo((xInd, yInd))
            for xInd in range(llPatchInd[0], urPatchInd[0]+1)
            for yInd in range(llPatchInd[1], urPatchInd[1]+1))

    def getBBox(self):
        """Get bounding box of tract (as an afwGeom.Box2I)
        """
        return afwGeom.Box2I(self._bbox)
    
    def getCtrCoord(self):
        """Get sky coordinate of center of tract (as an afwCoord.Coord)
        """
        return self._ctrCoord

    def getId(self):
        """Get ID of tract
        """
        return self._id
    
    def getNumPatches(self):
        """Get the number of patches in x, y
        
        @return the number of patches in x, y
        """
        return self._numPatches

    def getPatchBorder(self):
        """Get batch border
        
        @return patch border (pixels)
        """
        return self._patchBorder
    
    def getPatchInfo(self, index):
        """Return information for the specified patch
        
        @param[in] index: index of patch, as a pair of ints
        @return patch info, an instance of PatchInfo

        @raise IndexError if index is out of range
        """
        if (not 0 <= index[0] < self._numPatches[0]) \
            or (not 0 <= index[1] < self._numPatches[1]):
            raise IndexError("Patch index %s is not in range [0-%d, 0-%d]" % \
                (index, self._numPatches[0]-1, self._numPatches[1]-1))
        innerMin = afwGeom.Point2I(*[index[i] * self._patchInnerDimensions[i] for i in range(2)])
        innerBBox = afwGeom.Box2I(innerMin, self._patchInnerDimensions)
        if not self._bbox.contains(innerBBox):
            raise RuntimeError(
                "Bug: patch index %s valid but inner bbox=%s not contained in tract bbox=%s" % \
                (index, innerBBox, self._bbox))
        outerBBox = afwGeom.Box2I(innerBBox)
        outerBBox.grow(self.getPatchBorder())
        outerBBox.clip(self._bbox)
        return PatchInfo(
            index = index,
            innerBBox = innerBBox,
            outerBBox = outerBBox,
        )
    
    def getPatchInnerDimensions(self):
        """Get dimensions of inner region of the patches (all are the same)
        
        @return dimensions of inner region of the patches (as an afwGeom Extent2I)
        """
        return self._patchInnerDimensions
    
    def getTractOverlap(self):
        """Get minimum overlap of adjacent sky tracts
        
        @return minimum overlap between adjacent sky tracts, as an afwGeom Angle
        """
        return self._tractOverlap

    def getVertexList(self):
        """Get list of sky coordinates of vertices that define the boundary of the inner region
        
        @warning: this is not a deep copy
        @warning vertexCoordList will likely become a geom SphericalConvexPolygon someday.
        """
        return self._vertexCoordList
    
    def getWcs(self):
        """Get WCS of tract
        
        @warning: this is not a deep copy
        """
        return self._wcs

    def __str__(self):
        return "TractInfo(id=%s)" % (self._id,)
    
    def __repr__(self):
        return "TractInfo(id=%s, ctrCoord=%s)" % (self._id, self._ctrCoord.getVector())
