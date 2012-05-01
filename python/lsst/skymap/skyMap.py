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
- The sky tracts could be pentagonal (or some approximation), not rectangular
  and still preserve the desired minimum overlap. This would cut down a bit on the number
  of pixels to store per tract, but would make the code more complex.
"""
import math
import numpy

import lsst.afw.coord as afwCoord
import lsst.afw.geom as afwGeom
from . import detail
from .skyTractInfo import SkyTractInfo

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
    @param[in] defRA: the RA to use if the vector is too near a pole (an afwGeom Angle);
                ignored if not near a pole
    
    @throw RuntimeError if vec too near a pole and defRA is None
    """
    if abs(vec[0]) < _TinyFloat and abs(vec[1]) < _TinyFloat:
        if defRA is None:
            raise RuntimeError("At pole and defRA==None")
        if vec[2] > 0:
            dec = 90.0
        else:
            dec = -90.0
        return afwCoord.makeCoord(afwCoord.ICRS, defRA, afwGeom.Angle(dec, afwGeom.degrees))
    return afwCoord.makeCoord(afwCoord.ICRS, afwGeom.Point3D(*vec))
        

class SkyMap(object):
    """Metadata about a sky map pixelization.
        
    SkyMap divides the sky into 12 overlapping SkyTracts arranged as the tracts of a dodecahedron.
    Each sky tract is a rectangular Exposure using the specified WCS projection and nominal pixel scale.
    
    @note
    - The inner region of each sky tract is defined to be the region closer to the center of that tract
      than to the center of any other tract.
    - The native coordinate system is ICRS
    """
    def __init__(self,
        overlap = _DefaultOverlap,
        pixelScale = _DefaultPlateScale,
        projection = "STG",
        withTractsOnPoles = False,
    ):
        """Construct a SkyMap

        @param[in] overlap: minimum overlap between adjacent sky tracts; an afwGeom.Angle
        @param[in] pixelScale: nominal pixel scale (angle on sky/pixel); an afwGeom.Angle
        @param[in] projection: one of the FITS WCS projection codes, such as:
          - STG: stereographic projection
          - MOL: Molleweide's projection
          - TAN: tangent-plane projection
        @param[in] withTractsOnPoles: if True center a tract on each pole, else put a vertex on each pole
        """
        self._overlap = overlap
        self._pixelScale = pixelScale
        self._projection = str(projection)
        self._dodecahedron = detail.Dodecahedron(withFacesOnPoles = withTractsOnPoles)
        self._skyTractInfoList = []
        self._wcsFactory = detail.WcsFactory(self._pixelScale, self._projection)

        for id in range(12):
            tractVec = self._dodecahedron.getFace(id)
            tractCoord = _coordFromVec(tractVec, defRA=afwGeom.Angle(0))
            tractRA = tractCoord.getLongitude()
            vertexVecList = self._dodecahedron.getVertices(id)
            
            self._skyTractInfoList.append(SkyTractInfo(
                id = id,
                ctrCoord = tractCoord,
                vertexCoordList = [_coordFromVec(vec, defRA=tractRA) for vec in vertexVecList],
                overlap = self._overlap,
                wcsFactory = self._wcsFactory,
            ))
            
    def getOverlap(self):
        """Get the minimum overlap between adjacent sky tracts; an afwGeom.Angle
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

    def getSkyTractId(self, coord):
        """Return the ID of the sky tract whose inner region includes the coord.
        
        @param[in] coord: sky coordinate (afwCoord.Coord)
        
        @note
        - This routine will be more efficient if coord is ICRS.
        - If coord is equidistant between two sky tract centers then one of the two tracts
          is arbitrarily chosen (in an explicitly undefined manner).
        """
        return self._dodecahedron.getFaceInd(coord.toIcrs().getVector())
    
    def getNumSkyTracts(self):
        """Get the number of sky tracts
        """
        return len(self._skyTractInfoList)

    def getSkyTractInfo(self, id):
        """Get information about a sky tract
        
        @param[in] id: sky tract ID
        """
        return self._skyTractInfoList[id]
