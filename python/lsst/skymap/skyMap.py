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
- Tweak pixel scale so the average scale is as specified, rather than the scale at the center
- Bug fix: overlap=0 results in nonsense
- The outer boundary of the sky tile should be pentagonal, not rectangular; deal with it
  and in particular use this knowledge to compute the number of pixels per tile
"""
import math
import lsst.afw.coord as afwCoord
import lsst.afw.geom as afwGeom
import detail
import skyTileInfo

# LSST plate scale is 50 um/arcsec
# LSST pixel size is 10 um
# Default sky pixel scale is 1/sqrt(2) of image pixel scale
_DefaultPlateScale = 10.0 / (50.0 * 3600.0 * math.sqrt(2.0))

_RadPerDeg = math.pi / 180.0

class SkyMap(object):
    """Information about sky tiles
    """
    def __init__(self,
        overlap = 3.5 * _RadPerDeg,
        pixelScale = _DefaultPlateScale,
        projection = "STG",
        withFacesOnPoles = False,
    ):
        """Inputs:
        - overlap: minimum overlap between adjacent sky tiles (rad)
        - pixelScale: nominal pixel scale in degrees/pixel
        - projection: one of the FITS WCS projection codes, such as:
          - STG: stereographic projection
          - MOL: Molleweide's projection
          - TAG: tangent-plane projection
        - withFacesOnPoles: if True center a face on each pole, else put a vertex on each pole
        """
        self._overlap = float(overlap)
        self._pixelScale = float(pixelScale)
        self._projection = str(projection)
        self._dodecahedron = detail.Dodecahedron(withFacesOnPoles)
        self._skyTileInfoList = []
        self._wcsFactory = detail.WcsFactory(self._pixelScale, self._projection)

        def coordFromVec(vec, defRA=None):
            if abs(vec[0]) < 1e-14 and abs(vec[1]) < 1e-14:
                if defRA == None:
                    raise RuntimeError("At pole and defRA==None")
                if vec[2] > 0:
                    dec = 90.0
                else:
                    dec = -90.0
                return afwCoord.makeCoord(afwCoord.ICRS, afwGeom.makePointD(defRA, dec), afwCoord.DEGREES)
            return afwCoord.makeCoord(afwCoord.ICRS, afwGeom.makePointD(*vec))
        
        for ind in range(12):
            faceVec = self._dodecahedron.getFace(ind)
            faceCoord = coordFromVec(faceVec)
            faceRA = faceCoord.getPosition(afwCoord.DEGREES)[0]
            vertexVecList = self._dodecahedron.getVertices(ind)
            
            self._skyTileInfoList.append(skyTileInfo.SkyTileInfo(
                ind = ind,
                ctrCoord = faceCoord,
                vertexCoordList = [coordFromVec(vec, defRA=faceRA) for vec in vertexVecList],
                overlap = self._overlap,
                wcsFactory = self._wcsFactory,
            ))
            
    def getOverlap(self):
        """Get the minimum overlap between adjacent sky tiles (rad)
        """
        return self._overlap
    
    def getPixelScale(self):
        """Get the pixel scale in degrees/pixel
        """
        return self._pixelScale
    
    def getProjection(self):
        """Get the projection as a FITS WCS code
        """
        return self._projection

    def getSkyTileInd(self, coord):
        """Return the index of the sky tile whose optimal area includes the coord.
        
        If coord is on a boundary between two sky tiles then one of the two tiles will be returned
        without warning; it is explicitly not defined how the choice is made.
        """
        return self._dodecahedron.getFaceInd(coord.getVector())

    def getSkyTileInfo(self, ind):
        """Get information about a sky tile
        """
        return self._skyTileInfoList[ind]
