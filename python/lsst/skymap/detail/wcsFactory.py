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

import lsst.daf.base as dafBase
import lsst.afw.coord as afwCoord
import lsst.afw.geom as afwGeom
import lsst.afw.image as afwImage

class WcsFactory(object):
    """A factory for creating Wcs objects for the sky tiles.
    """
    def __init__(self, pixelScale, projection, rotation=0*afwGeom.radians):
        """Make a WcsFactory
        
        @param[in] pixelScale: desired scale, as sky/pixel, an afwGeom.Angle
        @param[in] projection: FITS-standard 3-letter name of projection, e.g.:
            TAN (tangent), STG (stereographic), MOL (Mollweide's), AIT (Hammer-Aitoff)
            see Representations of celestial coordinates in FITS (Calabretta and Greisen, 2002)
        @param[in] rotation:   Rotation relative to cardinal, as an lsst.afw.geom.Angle
        """
        if len(projection) != 3:
            raise RuntimeError("projection=%r; must have length 3" % (projection,))
        self._pixelScaleDeg = pixelScale.asDegrees()
        self._projection = str(projection)
        self._rotation = rotation
        cosTerm = self._pixelScaleDeg * math.cos(rotation.asRadians())
        sinTerm = self._pixelScaleDeg * math.sin(rotation.asRadians())
        self._cdMatrix = {"CD1_1": -cosTerm,
                          "CD2_1": sinTerm,
                          "CD1_2": sinTerm,
                          "CD2_2": cosTerm,
                          }
        self._ctypes = [("%-5s%3s" % (("RA", "DEC")[i], self._projection)).replace(" ", "-")
            for i in range(2)]

    def makeWcs(self, crPixPos, crValCoord, **kargs):
        """Make a Wcs
        
        @param[in] crPixPos: crPix for WCS, using the LSST standard; an afwGeom.Point2D or pair of floats
        @param[in] crValCoord: crVal for WCS (afwCoord.Coord)
        **kargs: FITS keyword arguments for WCS
        """
        ps = dafBase.PropertySet()
        crPixFits = [ind + 1.0 for ind in crPixPos] # convert pix position to FITS standard
        crValDeg = crValCoord.getPosition(afwGeom.degrees)
        for i in range(2):
            ip1 = i + 1
            ps.add("CTYPE%1d" % (ip1,), self._ctypes[i])
            ps.add("CRPIX%1d" % (ip1,), crPixFits[i])
            ps.add("CRVAL%1d" % (ip1,), crValDeg[i])
        ps.add("RADECSYS", "ICRS")
        ps.add("EQUINOX", 2000)
        for k,v in self._cdMatrix.items():
            ps.add(k, v)
        for key, value in kargs.iteritems():
            ps.add(key, value)
        return afwImage.makeWcs(ps)
