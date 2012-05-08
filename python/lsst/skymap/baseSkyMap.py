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
import lsst.afw.geom as afwGeom
from . import detail

__all__ = ["BaseSkyMap"]

class BaseSkyMap(object):
    """A collection of overlapping Tracts that map part or all of the sky.
    
    See TractInfo for more information.
    
    BaseSkyMap is an abstract base class. Subclasses must do the following:
    @li define __init__ and have it construct the TractInfo objects and put them in _tractInfoList
    @li define __getstate__ and __setstate__ to allow pickling (the butler saves sky maps using pickle);
        see DodecaSkyMap for an example of how to do this.
    """
    def __init__(self,
        patchInnerDimensions,
        patchBorder,
        tractOverlap,
        pixelScale,
        projection,
    ):
        """Construct a BaseSkyMap
        
        @warning: inheriting classes must set self._tractInfoList

        @param[in] patchInnerDimensions: dimensions of inner region of patches (x,y pixels)
        @param[in] patchBorder: border between patch inner and outer bbox (pixels); an int
        @param[in] tractOverlap: minimum overlap between adjacent sky tracts, on the sky; an afwGeom.Angle
        @param[in] pixelScale: nominal pixel scale (angle on sky/pixel); an afwGeom.Angle
        @param[in] projection: one of the FITS WCS projection codes, such as:
          - STG: stereographic projection
          - MOL: Molleweide's projection
          - TAN: tangent-plane projection
        """
        try:
            assert len(patchInnerDimensions) == 2
            self._patchInnerDimensions = afwGeom.Extent2I(*(int(val) for val in patchInnerDimensions))
        except Exception:
            raise RuntimeError("patchInnerDimensions = %r must contain 2 ints" % (patchInnerDimensions,))
        self._patchBorder = int(patchBorder)
        self._tractOverlap = tractOverlap
        self._pixelScale = pixelScale
        self._projection = str(projection)
        self._tractInfoList = []
        self._wcsFactory = detail.WcsFactory(self._pixelScale, self._projection)

    def getPatchBorder(self):
        """Get the border between the inner and outer bbox of patches (pixels)
        """
        return self._patchBorder
    
    def getPatchInnerDimensions(self):
        """Get dimensions of inner region of the patches (all are the same)
        
        @return dimensions of inner region of the patches (as an afwGeom Extent2I)
        """
        return self._patchInnerDimensions
    
    def getPixelScale(self):
        """Get the nominal pixel scale (angle on sky/pixel); an afwGeom.Angle
        """
        return self._pixelScale
    
    def getProjection(self):
        """Get the projection as a FITS WCS code
        """
        return self._projection

    def getTractOverlap(self):
        """Get the minimum overlap between adjacent sky tracts, on the sky; an afwGeom.Angle
        """
        return self._tractOverlap
    
    def findTract(self, coord):
        """Find the tract whose center is nearest the specified coord.
        
        @param[in] coord: sky coordinate (afwCoord.Coord)
        @return TractInfo of tract whose center is nearest the specified coord
        
        @warning:
        - if tracts do not cover the whole sky then the returned tract may not include the coord
        
        @note
        - This routine will be more efficient if coord is ICRS.
        - If coord is equidistant between multiple sky tract centers then one is arbitrarily chosen.
        - The default implementation is not very efficient; subclasses may wish to override.
        """
        icrsCoord = coord.toIcrs()
        distTractInfoList = []
        for tractInfo in self._tractInfoList:
            angSep = icrsCoord.angularSeparation(tractInfo.getCtrCoord()).asDegrees()
            distTractInfoList.append((angSep, tractInfo))
        distTractInfoList.sort()
        return distTractInfoList[0][1]
    
    def __getitem__(self, ind):
        return self._tractInfoList[ind]
    
    def __iter__(self):
        return iter(self._tractInfoList)
    
    def __len__(self):
        return len(self._tractInfoList)
