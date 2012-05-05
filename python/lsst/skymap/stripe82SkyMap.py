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

# yields patches similar in size to the existing FPC files
_DefaultNumPatches = (200, 10)

# match SDSS normal images; roughly 50"
_DefaultPatchBorder = 128

_DefaultNumTracts = 4


class Stripe82SkyMap(BaseSkyMap):
    """Cylindrical sky map pixelization for SDSS stripe 82 data.
        
    Stripe82SkyMap divides the sky into 4 overlapping equatorial tracts
    """
    def __init__(self,
        numPatches = _DefaultNumPatches,
        patchBorder = _DefaultPatchBorder,
        tractOverlap = _DefaultTractOverlap,
        pixelScale = _DefaultPlateScale,
        projection = _DefaultProjection,
        decRange = _DefaultDecRange,
        numTracts = _DefaultNumTracts,
    ):
        """Construct a Stripe82SkyMap

        @param[in] numPatches: number of patches in a tract along the (x=RA, y=Dec) direction
        @param[in] patchBorder: border between patch inner and outer bbox (pixels); an int
        @param[in] tractOverlap: minimum overlap between adjacent sky tracts; an afwGeom.Angle
        @param[in] pixelScale: nominal pixel scale (angle on sky/pixel); an afwGeom.Angle
        @param[in] projection: one of the FITS WCS projection codes
        @param[in] numTracts: number of tracts along RA (there is only one along Dec)
        @param[in] decRange: range of declination (a pair of afwGeom.Angle)
        """
        self._version = (1, 0) # for pickle
        self._numTracts = numTracts
        try:
            assert(len(decRange) == 2)
            junkRad = [ang.asRadians() for ang in decRange] # test elements are Angles
            self._decRange = tuple(decRange)
        except Exception:
            raise RuntimeError("decRange = %s; must be a pair of afwGeom Angles" % (decRange,))

        BaseSkyMap.__init__(self,
            numPatches = numPatches,
            patchBorder = patchBorder,
            tractOverlap = tractOverlap,
            pixelScale = pixelScale,
            projection = projection,
        )
    
    def __getstate__(self):
        """Support pickle
        
        @note: angle arguments are persisted in radians
        """
        return dict(
            version = self._version,
            numPatches = self.getNumPatches(),
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
    
    def _makeTracts(self):
        """Construct _skyTractInfoList
        """
        midDec = (self._decRange[1] - self._decRange[0]) / 2.0
        tractWidth = afwGeom.Angle(360.0 / self._numTracts, afwGeom.degrees)
        for id in range(4):
            begRA = tractWidth * id
            endRA = begRA + tractWidth
            vertexCoordList = (
                afwCoord.IcrsCoord(begRA, self._decRange[0]),
                afwCoord.IcrsCoord(endRA, self._decRange[0]),
                afwCoord.IcrsCoord(endRA, self._decRange[1]),
                afwCoord.IcrsCoord(begRA, self._decRange[1]),
            )

            midRA = begRA + tractWidth / 2.0
            ctrCoord = afwCoord.IcrsCoord(midRA, midDec)
                
            self._skyTractInfoList.append(TractInfo(
                id = id,
                numPatches = self.getNumPatches(),
                patchBorder = self.getPatchBorder(),
                ctrCoord = ctrCoord,
                vertexCoordList = vertexCoordList,
                tractOverlap = self.getTractOverlap(),
                wcsFactory = self._wcsFactory,
            ))

    def getDecRange(self):
        return self._decRange
