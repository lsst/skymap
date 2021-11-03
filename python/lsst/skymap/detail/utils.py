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

__all__ = ["coordFromVec", "makeSkyPolygonFromBBox", "Index2D"]

from typing import NamedTuple
import numpy

import lsst.sphgeom
import lsst.geom as geom

_TinyFloat = numpy.finfo(float).tiny


def coordFromVec(vec, defRA=None):
    """Convert an ICRS cartesian vector to an ICRS lsst.geom.SpherePoint

    Parameters
    ----------
    vec : `list` of `float`
        An ICRS catesian vector.
    defRA : `lsst.geom.Angle` or None
        The RA to use if the vector is too near a pole;
        ignored if not near a pole.

    Raises
    ------
    RuntimeError
        If vec too near a pole and defRA is None.
    """
    if abs(vec[0]) < _TinyFloat and abs(vec[1]) < _TinyFloat:
        if defRA is None:
            raise RuntimeError("At pole and defRA==None")
        if vec[2] > 0:
            decDeg = 90.0
        else:
            decDeg = -90.0
        return geom.SpherePoint(defRA, decDeg*geom.degrees)
    return geom.SpherePoint(lsst.sphgeom.Vector3d(*vec))


def makeSkyPolygonFromBBox(bbox, wcs):
    """Make an on-sky polygon from a bbox and a SkyWcs

    Parameters
    ----------
    bbox : `lsst.geom.Box2I` or `lsst.geom.Box2D`
        Bounding box of region, in pixel coordinates
    wcs : `lsst.afw.geom.SkyWcs`
        Celestial WCS

    Returns
    -------
    polygon : `lsst.sphgeom.ConvexPolygon`
        On-sky region
    """
    pixelPoints = geom.Box2D(bbox).getCorners()
    skyPoints = wcs.pixelToSky(pixelPoints)
    return lsst.sphgeom.ConvexPolygon([sp.getVector() for sp in skyPoints])


class Index2D(NamedTuple):
    """Two dimensional index for patches and cells.

    This class contains the x and y values of the location of a patch
    within a tract, or a cell within a patch.

    Parameters
    ----------
    x : `int`
    y : `int`
    """
    x: int
    y: int
