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
import lsst.pex.config as pexConfig
import lsst.afw.geom as afwGeom
from . import detail

__all__ = ["BaseSkyMap"]

class BaseSkyMapConfig(pexConfig.Config):
    patchInnerDimensions = pexConfig.ListField(
        doc = "dimensions of inner region of patches (x,y pixels)",
        dtype = int,
        length = 2,
        default = (4000, 4000),
    )
    patchBorder = pexConfig.Field(
        doc = "border between patch inner and outer bbox (pixels)",
        dtype = int,
        default = 100,
    )
    tractOverlap = pexConfig.Field(
        doc = "minimum overlap between adjacent sky tracts, on the sky (deg)",
        dtype = float,
        default = 1.0,
    )
    pixelScale = pexConfig.Field(
        doc = "nominal pixel scale (arcsec/pixel)",
        dtype = float,
        default = 0.333
    )
    projection = pexConfig.Field(
        doc = """one of the FITS WCS projection codes, such as:
          - STG: stereographic projection
          - MOL: Molleweide's projection
          - TAN: tangent-plane projection
        """,
        dtype = str,
        default = "STG",
    )
        

class BaseSkyMap(object):
    """A collection of overlapping Tracts that map part or all of the sky.
    
    See TractInfo for more information.
    
    BaseSkyMap is an abstract base class. Subclasses must do the following:
    @li define __init__ and have it construct the TractInfo objects and put them in _tractInfoList
    @li define __getstate__ and __setstate__ to allow pickling (the butler saves sky maps using pickle);
        see DodecaSkyMap for an example of how to do this. (Most of that code could be moved
        into this base class, but that would make it harder to handle older versions of pickle data.)
    """
    ConfigClass = BaseSkyMapConfig
    def __init__(self, config=None):
        """Construct a BaseSkyMap
        
        @param[in] config: an instance of self.ConfigClass; if None the default config is used
        """
        if config is None:
            config = self.ConfigClass()
        self.config = config
        self._tractInfoList = []
        self._wcsFactory = detail.WcsFactory(
            pixelScale = afwGeom.Angle(self.config.pixelScale, afwGeom.arcseconds),
            projection = self.config.projection,
        )
    
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
