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
"""
from . import detail

__all__ = ["BaseSkyMap"]

class BaseSkyMap(object):
    """A collection of overlapping Sky Tracts that map part or all of the sky.
    
    Each Sky Tract may be further divided into subregions.
        
    @note
    - The native coordinate system is ICRS.
    - The inner region of each sky tract is defined to be the region closer to the center of that tract
      than to the center of any other tract.
    """
    def __init__(self,
        numPatches,
        overlap,
        pixelScale,
        projection,
    ):
        """Construct a BaseSkyMap

        @param[in] numPatches: number of patches in a tract along the (x, y) direction
        @param[in] overlap: minimum overlap between adjacent sky tracts; an afwGeom.Angle
        @param[in] pixelScale: nominal pixel scale (angle on sky/pixel); an afwGeom.Angle
        @param[in] projection: one of the FITS WCS projection codes, such as:
          - STG: stereographic projection
          - MOL: Molleweide's projection
          - TAN: tangent-plane projection
        """
        try:
            assert len(numPatches) == 2
            self._numPatches = tuple(int(val) for val in numPatches)
        except Exception:
            raise RuntimeError("numPatches = %r must contain 2 ints" % (numPatches,))
        self._overlap = overlap
        self._pixelScale = pixelScale
        self._projection = str(projection)
        self._skyTractInfoList = []
        self._wcsFactory = detail.WcsFactory(self._pixelScale, self._projection)
        
        self._makeTracts()
    
    def _makeTracts(self):
        """Construct the list of tracts
        """
        raise NotImplementedError()
    
    def getNumPatches(self):
        """Return the number of patches within a tract along x, y
        """
        return self._numPatches

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
    
    def findTract(self, coord):
        """Find the tract whose center is nearest the specified coord.
        
        @param[in] coord: sky coordinate (afwCoord.Coord)
        @return SkyTractInfo of tract whose center is nearest the specified coord
        
        @warning:
        - if tracts do not cover the whole sky then the returned tract may not include the coord
        
        @note
        - This routine will be more efficient if coord is ICRS.
        - If coord is equidistant between multiple sky tract centers then one is arbitrarily chosen.
        - The default implementation is not very efficient; subclasses may wish to override.
        """
        icrsCoord = coord.toIcrs()
        distTractInfoList = []
        for tractInfo in self._skyTractInfoList:
            angSep = icrsCoord.angularSeparation(tractInfo.getCtrCoord()).asDegrees()
            distTractInfoList.append((angSep, tract))
        distTractInfoList.sort()
        return distTractInfoList[1]
    
    def findTractAndPatchIndex(self, coord):
        """Find tract whose center is nearest the specified coord and the index within that tract
        
        @param[in] coord: sky coordinate (afwCoord.Coord)
        @raise RuntimeError if coord is not on any tract. This is only possible if the tracts
        do not cover the entire sky.
        """
        icrsCoord = coord.toIcrs()
        tractInfo = self.findNearestTract(icrsCoord)
        return (tractInfo, tractInfo.getPatchIndex(icrsCoord))
    
    def __getitem__(self, ind):
        return self._skyTractInfoList[ind]
    
    def __iter__(self):
        return iter(self._skyTractInfoList)
    
    def __len__(self):
        return len(self._skyTractInfoList)
