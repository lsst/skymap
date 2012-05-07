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
import math
import lsst.afw.coord as afwCoord
import lsst.afw.geom as afwGeom
import lsst.afw.image as afwImage
from .patchInfo import PatchInfo

__all__ = ["TractInfo"]

class TractInfo(object):
    """Information about a tract in a SkyMap sky pixelization
    
    The tract is subdivided into rectangular patches. Each patch has the following properties:
    - An inner region defined by an inner bounding. The inner regions of the patches exactly tile the tract,
      and all inner regions have the same dimensions (the tract is made larger as required to make this true).
    - An outer region defined by an outer bounding box. The outer region extends beyond the inner region
      by patchBorder pixels in all directions the border is 0 at the edges of the tract. Thus patches
      overlap each other but never extend off the tract. Note that patchBorder may be 0.
    - An index that consists of a pair of integers:
        0 <= x index < numPatches[0]
        0 <= y index < numPatches[1]
      Patch 0,0 is at the minimum corner of the tract bounding box.
    """
    def __init__(self, id, numPatches, patchBorder, ctrCoord, vertexCoordList, tractOverlap, wcsFactory):
        """Construct a TractInfo

        @param[in] id: tract ID
        @param[in] numPatches: number of patches in a tract along the (x, y) direction
        @param[in] patchBorder: overlap between adjacent patches, in pixels; an int
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
#         print "TractInfo(id=%s, ctrCoord=%s, tractOverlap=%0.1f)" % \
#             (id, ctrCoord.getPosition(afwGeom.degrees), tractOverlap)
        self._id = id
        self._numPatches = numPatches
        self._patchBorder = patchBorder
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

        # grow initialBBox to hold exactly a multiple of numPatches, while keeping the center roughly the same
        initialBBoxMin = initialBBox.getMin()
        bboxDim = initialBBox.getDimensions()
        self._patchInnerDim = afwGeom.Extent2I(0, 0)
        for i, numPatches in enumerate(self._numPatches):
            self._patchInnerDim[i] = bboxDim[i] // numPatches
            bboxDim[i] = self._patchInnerDim[i] * numPatches
        initialBBox = afwGeom.Box2I(initialBBoxMin, bboxDim)

        # compute final bbox the same size but with LL corner = 0,0; use that to compute final WCS
        self._bbox = afwGeom.Box2I(afwGeom.Point2I(0, 0), initialBBox.getDimensions())
        # crpix for final Wcs is shifted by the same amount as bbox
        pixPosOffset = afwGeom.Extent2D(self._bbox.getMinX() - initialBBox.getMinX(),
                                        self._bbox.getMinY() - initialBBox.getMinY())
        crPixPos = initialCRPixPos + pixPosOffset
        self._wcs = wcsFactory.makeWcs(crPixPos=crPixPos, crValCoord=self._ctrCoord)
    
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
    
    def getPatchInnerDim(self):
        """Get dimensions of inner region of the patches (all are the same)
        
        @return dimensions of inner region of the patches (as an afwGeom Extent2I)
        """
        return self._patchInnerDim
    
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
        innerMin = afwGeom.Point2I(*[index[i] * self._patchInnerDim[i] for i in range(2)])
        innerBBox = afwGeom.Box2I(innerMin, self._patchInnerDim)
        if not self._bbox.contains(innerBBox):
            raise RuntimeError(
                "Bug: patch index %s valid but inner bbox=%s not contained in tract bbox=%s" % \
                (index, innerBBox, self._bbox))
        outerBBox = afwGeom.Box2I(innerBBox)
        outerBBox.grow(self._patchBorder)
        outerBBox.clip(self._bbox)
        return PatchInfo(
            index = index,
            innerBBox = innerBBox,
            outerBBox = outerBBox,
        )

    def findPatch(self, coord):
        """Find the patch containing the specified coord

        @param[in] coord: sky coordinate (afwCoord.Coord)
        @return PatchInfo of patch whose inner bbox contains the specified coord
        
        @raise RuntimeError if coord is not in tract
        
        @note This routine will be more efficient if coord is ICRS.
        """
        pixelInd = afwGeom.Point2I(self.getWcs().skyToPixel(coord.toIcrs()))
        if not self.getBBox().contains(pixelInd):
            raise RuntimeError("coord %s is not in tract %s" % (coord, self.getId()))
        patchInd = tuple(int(pixelInd[i]/self._patchInnerDim[i]) for i in range(2))
        return self.getPatchInfo(patchInd)
    
    def getTractOverlap(self):
        """Get minimum overlap of adjacent sky tracts
        
        @return minimum overlap between adjacent sky tracts, as an afwGeom Angle
        """
        return self._tractOverlap
    
    def getWcs(self):
        """Get WCS of tract
        
        @warning: this is not a deep copy
        """
        return self._wcs

    def getVertexList(self):
        """Get list of sky coordinates of vertices that define the boundary of the inner region
        
        @warning: this is not a deep copy
        @warning vertexCoordList will likely become a geom SphericalConvexPolygon someday.
        """
        return self._vertexCoordList

    def __str__(self):
        return "TractInfo(id=%s)" % (self._id,)
    
    def __repr__(self):
        return "TractInfo(id=%s, ctrCoord=%s)" % (self._id, self._ctrCoord.getVector())
