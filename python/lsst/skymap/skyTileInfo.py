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

_RadPerDeg = math.pi / 180.0

class SkyTileInfo(object):
    """Information about a sky tile in a SkyMap sky pixelization
    
    @todo Provide a way returning a geometry.SphericalConvexPolygon;
    one question is whether the geometry is ready; it certainly doesn't work with afwCoord yet.
    """
    def __init__(self, id, ctrCoord, vertexCoordList, overlap, wcsFactory):
        """Construct a SkyTileInfo

        @param[in] id: sky tile ID
        @param[in] ctrCoord: sky coordinate of center of inner region of tile, as an afwCoord.Coord;
            also used as the CRVAL for the WCS.
        @param[in] vertexCoordList: list of sky coordinates (afwCoord.Coord)
            of vertices that define the boundaries of the inner region
        @param[in] overlap: minimum overlap between adjacent sky tiles (rad);
            this defines the minimum distance the sky tile extends beyond the inner region in all directions
        @param[in] wcsFactory: a skymap.detail.WcsFactory object
        
        @warning
        - It is not enforced that ctrCoord is the center of vertexCoordList, but SkyMap relies on it
        - vertexCoordList will likely become a geom SphericalConvexPolygon someday.
        """
#         print "SkyTileInfo(id=%s, ctrCoord=%s, overlap=%0.1f)" % \
#             (id, ctrCoord.getPosition(afwCoord.DEGREES), overlap)
        self._id = id
        self._ctrCoord = ctrCoord
        self._vertexCoordList = tuple(coord.clone() for coord in vertexCoordList)
        self._overlap = float(overlap)
        
        DebugMinId = 0 # print extra information if id < DebugMinId

        # We don't know how big the tile will be, yet, so start by computing everything
        # as if the tile center was at pixel position 0, 0; then shift all pixel positions
        # so that the tile's bbox starts from 0,0
        initialCRPixPos = afwGeom.Point2D(0.0, 0.0)
        initialWcs = wcsFactory.makeWcs(crPixPos=initialCRPixPos, crValCoord=self._ctrCoord)

        # compute minimum bounding box that will hold all corners and overlap
        minBBoxD = afwGeom.Box2D()
        if id < DebugMinId:
            print "center position =", self._ctrCoord.getPosition(afwCoord.DEGREES)
        for vertexCoord in self._vertexCoordList:
            vertexDeg = vertexCoord.getPosition(afwCoord.DEGREES)
            if self._overlap == 0:
                minBBoxD.include(initialWcs.skyToPixel(vertexCoord))
            else:
                halfOverlap = self._overlap / 2.0
                numAngles = 24
                angleIncr = _RadPerDeg * 360.0 / float(numAngles)
                for i in range(numAngles):
                    offAngle = angleIncr * i
                    offCoord = vertexCoord.clone()
                    offCoord.offset(offAngle, halfOverlap)
                    pixPos = initialWcs.skyToPixel(offCoord)
                    if id < DebugMinId:
                        print "vertexDeg=%s, offDeg=%s, pixPos=%s" % \
                            (vertexDeg, offCoord.getPosition(afwCoord.DEGREES), pixPos)
                    minBBoxD.include(pixPos)
        initialBBox = afwGeom.Box2I(minBBoxD)

        # compute final bbox the same size but with LL corner = 0,0; use that to compute final WCS
        self._bbox = afwGeom.Box2I(afwGeom.Point2I(0, 0), initialBBox.getDimensions())
        # crpix for final Wcs is shifted by the same amount as bbox
        pixPosOffset = afwGeom.Extent2D(self._bbox.getMinX() - initialBBox.getMinX(),
                                        self._bbox.getMinY() - initialBBox.getMinY())
        crPixPos = initialCRPixPos + pixPosOffset
        if id < DebugMinId:
            print "initialBBox=%s; bbox=%s; pixPosOffset=%s; crPixPos=%s" % \
                (initialBBox, self._bbox, pixPosOffset, crPixPos)
        self._wcs = wcsFactory.makeWcs(crPixPos=crPixPos, crValCoord=self._ctrCoord)
    
    def getBBox(self):
        """Get bounding box of sky tile (as an afwGeom.Box2I)
        """
        return afwGeom.Box2I(self._bbox)
    
    def getCtrCoord(self):
        """Get sky coordinate of center of sky tile (as an afwCoord.Coord)
        """
        return self._ctrCoord.clone()

    def getId(self):
        """Get ID of sky tile
        """
        return self._id
    
    def getOverlap(self):
        """Get overlap with adjacent sky tiles (rad)
        """
        return self._overlap
    
    def getWcs(self):
        """Get WCS of sky tile
        
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
        return "SkyTileInfo(id=%s)" % (self._id,)
    
    def __repr__(self):
        return "SkyTileInfo(id=%s, ctrCoord=%s)" % (self._id, self._ctrCoord.getVector())
