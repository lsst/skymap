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
- The sky faces could be pentagonal (or some approximation), not rectangular
  and still preserve the desired minimum overlap. This would cut down a bit on the number
  of pixels to store per face, but would make the code more complex.
"""
import math
import numpy

import lsst.afw.coord as afwCoord
import lsst.afw.geom as afwGeom
from . import detail
from .skyFaceInfo import SkyFaceInfo

_TinyFloat = numpy.finfo(float).tiny

# Default overlap is 3.5 degrees
_DefaultOverlap = afwGeom.Angle(3.5, afwGeom.degrees)

# LSST plate scale is 50 um/arcsec
# LSST pixel size is 10 um
# Default sky pixel scale is 1/sqrt(2) of image pixel scale
_DefaultPlateScale = afwGeom.Angle(10.0 / (50.0 * math.sqrt(2.0)), afwGeom.arcseconds)

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
        
    SkyMap divides the sky into 12 overlapping SkyFaces arranged as the faces of a dodecahedron.
    Each sky face is a rectangular Exposure using the specified WCS projection and nominal pixel scale.
    
    @note
    - The inner region of each sky face is defined to be the region closer to the center of that face
      than to the center of any other face.
    - The native coordinate system is ICRS
    """
    def __init__(self,
        overlap = _DefaultOverlap,
        pixelScale = _DefaultPlateScale,
        projection = "STG",
        withFacesOnPoles = False,
    ):
        """Construct a SkyMap

        @param[in] overlap: minimum overlap between adjacent sky faces; an afwGeom.Angle
        @param[in] pixelScale: nominal pixel scale (angle on sky/pixel); an afwGeom.Angle
        @param[in] projection: one of the FITS WCS projection codes, such as:
          - STG: stereographic projection
          - MOL: Molleweide's projection
          - TAN: tangent-plane projection
        @param[in] withFacesOnPoles: if True center a face on each pole, else put a vertex on each pole
        """
        self._overlap = overlap
        self._pixelScale = pixelScale
        self._projection = str(projection)
        self._dodecahedron = detail.Dodecahedron(withFacesOnPoles)
        self._skyFaceInfoList = []
        self._wcsFactory = detail.WcsFactory(self._pixelScale, self._projection)

        for id in range(12):
            faceVec = self._dodecahedron.getFace(id)
            faceCoord = _coordFromVec(faceVec)
            faceRA = faceCoord.getPosition(afwGeom.degrees)[0]
            vertexVecList = self._dodecahedron.getVertices(id)
            
            self._skyFaceInfoList.append(SkyFaceInfo(
                id = id,
                ctrCoord = faceCoord,
                vertexCoordList = [_coordFromVec(vec, defRA=faceRA) for vec in vertexVecList],
                overlap = self._overlap,
                wcsFactory = self._wcsFactory,
            ))
            
    def getOverlap(self):
        """Get the minimum overlap between adjacent sky faces; an afwGeom.Angle
        """
        return self._overlap
    
    def getPixelScale(self):
        """Get the nominal pixel scale (angle on sky/pixel); an afwGeom.Angle
        """
        return self._pixelScale
    
    def getProjection(self):
        """Get the projection as a FITS WCS code
        """
        return self._projection

    def getSkyFaceId(self, coord):
        """Return the ID of the sky face whose inner region includes the coord.
        
        @param[in] coord: sky coordinate (afwCoord.Coord)
        
        @note
        - This routine will be more efficient if coord is ICRS.
        - If coord is equidistant between two sky face centers then one of the two faces
          is arbitrarily chosen (in an explicitly undefined manner).
        """
        return self._dodecahedron.getFaceInd(coord.toIcrs().getVector())
    
    def getNumSkyFaces(self):
        """Get the number of sky faces
        """
        return len(self._skyFaceInfoList)

    def getSkyFaceInfo(self, id):
        """Get information about a sky face
        
        @param[in] id: sky face ID
        """
        return self._skyFaceInfoList[id]
