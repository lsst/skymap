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
    """Information about a sky tile
    
    @todo: provide a way returning a geometry.SphericalConvexPolygon;
    one question is whether the geometry is ready; it certainly doesn't work with afwCoord yet.
    """
    def __init__(self, id, crValCoord, vertexCoordList, overlap, wcsFactory):
        """Construct a SkyTileInfo

        @param[in] id: sky tile ID
        @param[in] crValCoord: sky coordinate of center of WCS (CRVal), an afwCoord.Coord
        @param[in] vertexCoordList: list of sky coordinates of vertices that define edge of optimal area
        @param[in] overlap: minimum overlap between adjacent sky tiles (rad)
        @param[in] wcsFactory: a skymap.detail.WcsFactory object
        """
#         print "SkyTileInfo(id=%s, crValCoord=%s, overlap=%0.1f)" % \
#             (id, crValCoord.getPosition(afwCoord.DEGREES), overlap)
        self._id = id
        self._vertexCoordList = tuple(coord.clone() for coord in vertexCoordList)
        self._overlap = float(overlap)
        
        DebugMinId = 0 # print extra information if id < DebugMinId

        # We don't know how big the tile will be, yet, so start by computing everything
        # as if the tile center was at pixel position 0, 0; then shift all pixel positions
        # so that the tile's bbox starts from 0,0
        initialCRPixPos = afwGeom.Point2D(0.0, 0.0)
        initialWcs = wcsFactory.makeWcs(crPixPos=initialCRPixPos, crValCoord=crValCoord)

        # compute minimum bounding box that will hold all corners and overlap
        minBBoxD = afwGeom.Box2D()
        if id < DebugMinId:
            print "center position =", crValCoord.getPosition(afwCoord.DEGREES)
        for vertexCoord in self._vertexCoordList:
            vertexDeg = vertexCoord.getPosition(afwCoord.DEGREES)
            if self._overlap == 0:
                minBBoxD.include(initialWcs.skyToPixel(vertexCoord))
            else:
                for i in range(4):
                    offAngle = 90.0 * i
                    offCoord = vertexCoord.clone()
                    offCoord.offset(offAngle * _RadPerDeg, self._overlap / 2.0)
                    pixPos = initialWcs.skyToPixel(offCoord)
                    if id < DebugMinId:
                        print "vertexDeg=%s, offDeg=%s, pixPos=%s" % \
                            (vertexDeg, offCoord.getPosition(afwCoord.DEGREES), pixPos)
                    minBBoxD.include(pixPos)
        initialBBox = afwGeom.Box2I(minBBoxD)

        # compute final bbox the same size but with LL corner = 0,0; use that to compute final WCS
        self._bbox = afwGeom.Box2I(afwGeom.Point2I(0, 0), initialBBox.getDimensions())
        # crpix for final Wcs is shifted by the same amount as bbox
        crPixPos = initialCRPixPos + afwGeom.Extent2D(self._bbox.getMin() - initialBBox.getMin())
        if id < DebugMinId:
            print "initialBBox=%s; bbox=%s; crPixPos=%s" % (initialBBox, self._bbox, crPixPos)
        self._wcs = wcsFactory.makeWcs(crPixPos=crPixPos, crValCoord=crValCoord)
    
    def getBBox(self):
        """Get bounding box of sky tile (as an afwGeom.Box2I)
        """
        return afwGeom.Box2I(self._bbox)

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
        """Get list of sky coordinates of vertices that define edge of optimal area
        
        @warning: this is not a deep copy
        """
        return self._vertexCoordList
