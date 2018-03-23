from builtins import object
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
import lsst.afw.geom as afwGeom


class WcsFactory(object):
    """A factory for creating Wcs objects for the sky tiles.
    """

    def __init__(self, pixelScale, projection, rotation=0*afwGeom.radians, flipX=False):
        """Make a WcsFactory

        @param[in] pixelScale: desired scale, as sky/pixel, an afwGeom.Angle
        @param[in] projection: FITS-standard 3-letter name of projection, e.g.:
            TAN (tangent), STG (stereographic), MOL (Mollweide's), AIT (Hammer-Aitoff)
            see Representations of celestial coordinates in FITS (Calabretta and Greisen, 2002)
        @param[in] rotation:   Rotation relative to cardinal, as an lsst.afw.geom.Angle
        @param[in] flipX: Flip the X axis?
        """
        if len(projection) != 3:
            raise RuntimeError("projection=%r; must have length 3" % (projection,))
        self._projection = projection
        self._cdMatrix = afwGeom.makeCdMatrix(scale=pixelScale, orientation=rotation, flipX=flipX)

    def makeWcs(self, crPixPos, crValCoord):
        """Make a Wcs

        @param[in] crPixPos: crPix for WCS, using the LSST standard; an afwGeom.Point2D or pair of floats
        @param[in] crValCoord: ICRS crVal for WCS (lsst.afw.geom.SpherePoint)
        """
        return afwGeom.makeSkyWcs(crpix=crPixPos, crval=crValCoord,
                                  cdMatrix=self._cdMatrix, projection=self._projection)
