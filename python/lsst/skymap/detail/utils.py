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
import numpy

import lsst.afw.geom as afwGeom

__all__ = ["coordFromVec"]

_TinyFloat = numpy.finfo(float).tiny


def coordFromVec(vec, defRA=None):
    """Convert an ICRS cartesian vector to an ICRS lsst.afw.geom.SpherePoint

    @param[in] vec: an ICRS catesian vector as a sequence of three floats
    @param[in] defRA: the RA to use if the vector is too near a pole (an afwGeom Angle);
                ignored if not near a pole

    @throw RuntimeError if vec too near a pole and defRA is None
    """
    if abs(vec[0]) < _TinyFloat and abs(vec[1]) < _TinyFloat:
        if defRA is None:
            raise RuntimeError("At pole and defRA==None")
        if vec[2] > 0:
            decDeg = 90.0
        else:
            decDeg = -90.0
        return afwGeom.SpherePoint(defRA, decDeg*afwGeom.degrees)
    return afwGeom.SpherePoint(afwGeom.Point3D(*vec))
