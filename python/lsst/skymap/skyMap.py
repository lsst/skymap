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
"""
@todo
- Consider tweaking pixel scale so the average scale is as specified, rather than the scale at the center
- The sky tiles could be pentagonal (or some approximation), not rectangular
  and still preserve the desired minimum overlap. This would cut down a bit on the number
  of pixels to store per tile, but would make the code more complex.
"""
import math
import numpy

import lsst.afw.coord as afwCoord
import lsst.afw.geom as afwGeom
import detail
import skyTileInfo

_RadPerDeg = math.pi / 180.0
_TinyFloat = numpy.finfo(float).tiny

# Default overlap is 3.5 degrees
_DefaultOverlap = 3.5 * _RadPerDeg

# LSST plate scale is 50 um/arcsec
# LSST pixel size is 10 um
# Default sky pixel scale is 1/sqrt(2) of image pixel scale
# Required units are radians/pixel
_DefaultPlateScale = (10.0 / (50.0 * 3600.0 * math.sqrt(2.0))) * _RadPerDeg

def _coordFromVec(vec, defRA=None):
    """Convert an ICRS cartesian vector to an ICRS Coord
    
    @param[in] vec: an ICRS catesian vector as a sequence of three floats
    @param[in] defRA: the RA to use if the vector is too near a pole; ignored if not near a pole
    
    @throw RuntimeError if vec too near a pole and defRA is None
    """
    if abs(vec[0]) < _TinyFloat and abs(vec[1]) < _TinyFloat:
        if defRA == None:
            raise RuntimeError("At pole and defRA==None")
        if vec[2] > 0:
            dec = 90.0
        else:
            dec = -90.0
        return afwCoord.makeCoord(afwCoord.ICRS, afwGeom.Point2D(defRA, dec), afwCoord.DEGREES)
    return afwCoord.makeCoord(afwCoord.ICRS, afwGeom.Point3D(*vec))
        

class SkyMap(object):
    """Metadata about a sky map pixelization.
        
    SkyMap divides the sky into 12 overlapping SkyTiles arranged as the faces of a dodecahedron.
    Each sky tile is a rectangular Exposure using the specified WCS projection and nominal pixel scale.
    
    @note
    - The inner region of each sky tile is defined to be the region closer to the center of that tile
      than to the center of any other tile.
    - The native coordinate system is ICRS
    """
    def __init__(self,
        overlap = _DefaultOverlap,
        pixelScale = _DefaultPlateScale,
        projection = "STG",
        withFacesOnPoles = False,
    ):
        """Construct a SkyMap

        @param[in] overlap: minimum overlap between adjacent sky tiles (rad)
        @param[in] pixelScale: nominal pixel scale (rad/pixel)
        @param[in] projection: one of the FITS WCS projection codes, such as:
          - STG: stereographic projection
          - MOL: Molleweide's projection
          - TAN: tangent-plane projection
        @param[in] withFacesOnPoles: if True center a face on each pole, else put a vertex on each pole
        """
        self._overlap = float(overlap)
        self._pixelScale = float(pixelScale)
        self._projection = str(projection)
        self._dodecahedron = detail.Dodecahedron(withFacesOnPoles)
        self._skyTileInfoList = []
        self._wcsFactory = detail.WcsFactory(self._pixelScale, self._projection)

        for id in range(12):
            faceVec = self._dodecahedron.getFace(id)
            faceCoord = _coordFromVec(faceVec)
            faceRA = faceCoord.getPosition(afwCoord.DEGREES)[0]
            vertexVecList = self._dodecahedron.getVertices(id)
            
            self._skyTileInfoList.append(skyTileInfo.SkyTileInfo(
                id = id,
                ctrCoord = faceCoord,
                vertexCoordList = [_coordFromVec(vec, defRA=faceRA) for vec in vertexVecList],
                overlap = self._overlap,
                wcsFactory = self._wcsFactory,
            ))
            
    def getOverlap(self):
        """Get the minimum overlap between adjacent sky tiles (rad)
        """
        return self._overlap
    
    def getPixelScale(self):
        """Get the nominal pixel scale (rad/pixel)
        """
        return self._pixelScale
    
    def getProjection(self):
        """Get the projection as a FITS WCS code
        """
        return self._projection

    def getSkyTileId(self, coord):
        """Return the ID of the sky tile whose inner region includes the coord.
        
        @param[in] coord: sky coordinate (afwCoord.Coord)
        
        @note
        - This routine will be more efficient if coord is ICRS.
        - If coord is equidistant between two sky tile centers then one of the two tiles
          is arbitrarily chosen (in an explicitly undefined manner).
        """
        return self._dodecahedron.getFaceInd(coord.toIcrs().getVector())
    
    def getNumSkyTiles(self):
        """Get the number of sky tiles
        """
        return len(self._skyTileInfoList)

    def getSkyTileInfo(self, id):
        """Get information about a sky tile
        
        @param[in] id: sky tile ID
        """
        return self._skyTileInfoList[id]
