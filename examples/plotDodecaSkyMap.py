#!/usr/bin/env python
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
"""Display sky map geometry as a 3d plot
"""
import math
import numpy

from mpl_toolkits.mplot3d import Axes3D  # noqa F401 used by fig.gca
import matplotlib.pyplot as plt

import lsst.geom as geom
from lsst.skymap import DodecaSkyMap

skyMap = DodecaSkyMap()

fig = plt.figure()
ax = fig.gca(projection='3d')

for tractInfo in skyMap:
    # display inner region
    vertexList = list(tractInfo.getVertexList())
    vertexList.append(vertexList[0])  # to close region
    innerPoints = [tuple(coord.getVector()) for coord in vertexList]
    inX, inY, inZ = zip(*innerPoints)
    lineList = ax.plot(inX, inY, inZ, label="Inner tractInfo %s" % (tractInfo.getId(),))
    color = lineList[0].get_color()

    # display center
    centerPoint = numpy.mean(innerPoints[0:-1], axis=0)
    ax.plot([centerPoint[0]], [centerPoint[1]], [centerPoint[2]], ".", color=color)

    # display outer edge; scale to be approximately in the same plane as the inner region
    wcs = tractInfo.getWcs()
    bbox = tractInfo.getBBox()
    outerPixPos = [
        bbox.getMin(),
        geom.Point2I(bbox.getMaxX(), bbox.getMinY()),
        bbox.getMax(),
        geom.Point2I(bbox.getMinX(), bbox.getMaxY()),
        bbox.getMin(),
    ]
    outerPoints = [numpy.array(wcs.pixelToSky(p[0], p[1]).getVector()) for p in outerPixPos]
    innerCtrRadSq = numpy.sum(centerPoint**2)
    outerCtr = numpy.mean(outerPoints[0:-1], axis=0)
    outerCtrRadSq = numpy.sum(outerCtr**2)
    scale = math.sqrt(innerCtrRadSq / outerCtrRadSq)
    outerPoints = [v * scale for v in outerPoints]
    outX, outY, outZ = zip(*outerPoints)
    ax.plot(outX, outY, outZ, color=color)

plt.show()
