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
import lsst.afw.coord as afwCoord
import lsst.afw.geom as afwGeom
from .baseSkyMap import BaseSkyMap
from .tractInfo import TractInfo

# Default overlap is 50", which is approximately the overlap SDSS uses between images
_DefaultTractOverlap = afwGeom.Angle(50, afwGeom.arcseconds)

# SDSS plate scale is 0.396"/pixel
_DefaultPlateScale = afwGeom.Angle(0.396, afwGeom.arcseconds)

# SDSS stripe 82 covers a declination range of -1.25 to 1.25 degrees
_DefaultDecRange = (afwGeom.Angle(-1.25, afwGeom.degrees), afwGeom.Angle(1.25, afwGeom.degrees))

_DefaultProjection = "CEA"

_DefaultNumTracts = 2

# yields patches similar in size to the existing FPC files
_DefaultPatchInnerDimensions = (2048, 2048)

# match SDSS normal images; roughly 50"
_DefaultPatchBorder = 128


class EquatSkyMap(BaseSkyMap):
    """Equatorial sky map pixelization, e.g. for SDSS stripe 82 image data.
        
    EquatSkyMap represents an equatorial band of sky divided along declination into overlapping tracts.
    """
    def __init__(self,
        patchInnerDimensions = _DefaultPatchInnerDimensions,
        patchBorder = _DefaultPatchBorder,
        tractOverlap = _DefaultTractOverlap,
        pixelScale = _DefaultPlateScale,
        projection = _DefaultProjection,
        decRange = _DefaultDecRange,
        numTracts = _DefaultNumTracts,
    ):
        """Construct a EquatSkyMap

        @param[in] patchInnerDimensions: dimensions of inner region of patches (x,y pixels)
        @param[in] patchBorder: border between patch inner and outer bbox (pixels); an int
        @param[in] tractOverlap: minimum overlap between adjacent sky tracts; an afwGeom.Angle
        @param[in] pixelScale: nominal pixel scale (angle on sky/pixel); an afwGeom.Angle
        @param[in] projection: one of the FITS WCS projection codes
        @param[in] numTracts: number of tracts along RA (there is only one along Dec)
        @param[in] decRange: range of declination (a pair of afwGeom.Angle)
        
        @warning: some projections, such as "TAN", require at least 3 tracts.
        If you use too few then construction will fail.
        """
        self._version = (1, 0) # for pickle

        try:
            assert(len(decRange) == 2)
            junkRad = [ang.asRadians() for ang in decRange] # test elements are Angles
            self._decRange = tuple(decRange)
        except Exception:
            raise RuntimeError("decRange = %s; must be a pair of afwGeom Angles" % (decRange,))

        BaseSkyMap.__init__(self,
            patchInnerDimensions = patchInnerDimensions,
            patchBorder = patchBorder,
            tractOverlap = tractOverlap,
            pixelScale = pixelScale,
            projection = projection,
        )
    
        midDec = (self._decRange[0] + self._decRange[1]) / 2.0
        tractWidthRA = afwGeom.Angle(360.0 / numTracts, afwGeom.degrees)
        for id in range(numTracts):
            begRA = tractWidthRA * id
            endRA = begRA + tractWidthRA
            vertexCoordList = (
                afwCoord.IcrsCoord(begRA, self._decRange[0]),
                afwCoord.IcrsCoord(endRA, self._decRange[0]),
                afwCoord.IcrsCoord(endRA, self._decRange[1]),
                afwCoord.IcrsCoord(begRA, self._decRange[1]),
            )

            midRA = begRA + tractWidthRA / 2.0
            ctrCoord = afwCoord.IcrsCoord(midRA, midDec)
                
            self._tractInfoList.append(TractInfo(
                id = id,
                patchInnerDimensions = self.getPatchInnerDimensions(),
                patchBorder = self.getPatchBorder(),
                ctrCoord = ctrCoord,
                vertexCoordList = vertexCoordList,
                tractOverlap = self.getTractOverlap(),
                wcsFactory = self._wcsFactory,
            ))

    
    def __getstate__(self):
        """Support pickle
        
        @note: angle arguments are persisted in radians
        """
        return dict(
            version = self._version,
            patchInnerDimensions = tuple(self.getPatchInnerDimensions()),
            patchBorder = self.getPatchBorder(),
            tractOverlap = self.getTractOverlap().asRadians(),
            pixelScale = self.getPixelScale().asRadians(),
            projection = self.getProjection(),
            decRange = [ang.asRadians() for ang in self.getDecRange()],
            numTracts = len(self),
        )
    
    def __setstate__(self, argDict):
        """Support unpickle
        """
        version = argDict.pop("version")
        if version >= (2, 0):
            raise runtimeError("Version = %s >= (2,0); cannot unpickle" % (version,))
        for angleArg in ("tractOverlap", "pixelScale"):
            argDict[angleArg] = afwGeom.Angle(argDict[angleArg], afwGeom.radians)
        argDict["decRange"] = tuple(afwGeom.Angle(ang, afwGeom.radians) for ang in argDict["decRange"])
        self.__init__(**argDict)

    def getDecRange(self):
        """Return the declination range
        
        @return declination range as a pair of afwGeom.Angles
        """
        return self._decRange

    def getVersion(self):
        """Return version (e.g. for pickle)
        
        @return version as a pair of integers
        """
        return self._version
