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
import pickle
import argparse
from itertools import cycle

from mpl_toolkits.mplot3d import Axes3D  # noqa F401 used by fig.gca
import matplotlib.pyplot as plt

import lsst.afw.geom as afwGeom


def reportSkyMapInfo(skyMap):
    paramDict = skyMap.config.toDict()
    paramNameList = sorted(paramDict)
    print("Sky Map parameters:")
    for paramName in paramNameList:
        param = paramDict[paramName]
        print("skyMap.config.%s = %s" % (paramName, param))

    print("\nSkyMap has %d tracts:" % (len(skyMap)))
    for tractInfo in skyMap:
        wcs = tractInfo.getWcs()
        posBox = afwGeom.Box2D(tractInfo.getBBox())
        pixelPosList = (
            posBox.getMin(),
            afwGeom.Point2D(posBox.getMaxX(), posBox.getMinY()),
            posBox.getMax(),
            afwGeom.Point2D(posBox.getMinX(), posBox.getMaxY()),
        )
        skyPosList = [wcs.pixelToSky(pos).getPosition(afwGeom.degrees) for pos in pixelPosList]
        posStrList = ["(%0.3f, %0.3f)" % tuple(skyPos) for skyPos in skyPosList]
        print("tract %s has corners %s (RA, Dec deg) and %s x %s patches" %
              (tractInfo.getId(), ", ".join(posStrList),
               tractInfo.getNumPatches()[0], tractInfo.getNumPatches()[1]))


def plotSkyMap3d(skyMap):
    fig = plt.figure()
    ax = fig.gca(projection="3d")
    ax.set_axis_off()

    # make sure a complete 1x1x1 cube is shown -- what I really want is to constrain the aspect ratio
    # but that is not yet supported for 3D plots
    for direction in (-1, 1):
        for point in numpy.diag(direction * numpy.array([1, 1, 1])):
            ax.plot([point[0]], [point[1]], [point[2]], 'w')

    for tractInfo in skyMap:
        # display outer edge; scale to be approximately in the same plane as the inner region
        wcs = tractInfo.getWcs()
        posBox = afwGeom.Box2D(tractInfo.getBBox())
        xRange = posBox.getMinX(), posBox.getMaxX()
        yRange = posBox.getMinY(), posBox.getMaxY()

        numX = min(50, max(1, ((xRange[1] - xRange[0]) // 100)))
        numY = min(50, max(1, ((yRange[1] - yRange[0]) // 100)))

        outerPixPosList = \
            [(x1, yRange[0]) for x1 in numpy.linspace(xRange[0], xRange[1], num=numX, endpoint=False)] \
            + [(xRange[1], y1) for y1 in numpy.linspace(yRange[0], yRange[1], num=numY, endpoint=False)] \
            + [(x2, yRange[1]) for x2 in numpy.linspace(xRange[1], xRange[0], num=numX, endpoint=False)] \
            + [(xRange[0], y2) for y2 in numpy.linspace(yRange[1], yRange[0], num=numY, endpoint=False)]
        outerPixPosList.append(outerPixPosList[0])

        outerPoints = [numpy.array(wcs.pixelToSky(p[0], p[1]).getVector()) for p in outerPixPosList]
        outX, outY, outZ = zip(*outerPoints)
        ax.plot(outX, outY, outZ)

    plt.show()


class DefaultProjector:
    """Default projector for plotting a SkyMap in 2D

    Coordinates may be projected and optionally recentered. The recentering
    helps to avoid lines crisscrossing the plot when plotting boundaries.

    Class variables xLabel and yLabel provide labels for the plots.
    """
    xLabel = "RA (radians)"
    yLabel = "sin(Dec)"

    def __init__(self, coord):
        """Constructor

        Provided coordinate will be stored as the center of the area
        of interest, for subsequent recentering.
        """
        self._ra0 = coord.getLongitude().asRadians()
        self._dec0 = coord.getLatitude().asRadians()

    @staticmethod
    def project(coord):
        """Project the provided coordinates"""
        return coord.getLongitude().asRadians(), math.sin(coord.getLatitude().asRadians())

    def recenter(self, x, y):
        """Recenter the projected coordinates.

        The recentering seeks to keep boundaries on the same side of
        the plot as the center.
        """
        if x < 0.5*math.pi and self._ra0 > 1.5*math.pi:
            x += 2*math.pi
        elif x > 1.5*math.pi and self._ra0 < 0.5*math.pi:
            x -= 2*math.pi
        return x, y

    def projectWithRecenter(self, coord):
        """Project the coordinates and recenter"""
        x, y = self.project(coord)
        return self.recenter(x, y)


class PoleProjector(DefaultProjector):
    """A projection at the north pole"""
    xLabel = "x"
    yLabel = "y"

    def project(self, coord):
        ra, dec = coord.getLongitude().asRadians(), coord.getLatitude().asRadians()
        r = 0.5*math.pi - dec
        theta = ra
        return r*math.cos(theta), r*math.sin(theta)

    def recenter(self, x, y):
        """No recentering required."""
        return x, y


def makePlotter(Projector=DefaultProjector):
    """Make a function that will plot a SkyMap in 2D

    The Projector is used to project the center of each tract and its boundaries
    onto the plot.
    """

    def plotSkyMap2d(skyMap):
        plt.figure()
        plt.clf()
        axes = plt.axes()

        for tract, color in zip(skyMap, cycle("rgbcmyk")):
            center = tract.getCtrCoord()
            wcs = tract.getWcs()
            box = tract.getBBox()
            xMin, xMax, yMin, yMax = box.getMinX(), box.getMaxX(), box.getMinY(), box.getMaxY()
            num = 50
            proj = Projector(center)
            x, y = proj.project(center)
            axes.text(x, y, str(tract.getId()), color=color, ha="center", va="center")
            xList = numpy.linspace(xMin, xMax, num=num, endpoint=True)
            yList = numpy.linspace(yMin, yMax, num=num, endpoint=True)
            for xs, ys in ((xList, yMin*numpy.ones(num)),
                           (xMax*numpy.ones(num), yList),
                           (xList, yMax*numpy.ones(num)),
                           (xMin*numpy.ones(num), yList),
                           ):
                coords = [wcs.pixelToSky(afwGeom.Point2D(x1, y1)) for x1, y1 in zip(xs, ys)]
                bounds = [proj.projectWithRecenter(c) for c in coords]
                axes.plot([b[0] for b in bounds], [b[1] for b in bounds], color + '-')

        plt.xlabel(Projector.xLabel)
        plt.ylabel(Projector.yLabel)
        plt.grid(True)
        plt.show()
    return plotSkyMap2d


if __name__ == "__main__":
    plotStyles = {"3d": plotSkyMap3d,
                  "2d": makePlotter(),
                  "pole": makePlotter(PoleProjector),
                  }
    parser = argparse.ArgumentParser()
    parser.add_argument("skymap", nargs=1, help="Path to skymap pickle")
    parser.add_argument("--style", choices=list(plotStyles.keys()), default="3d", help="Plot style to use")
    args = parser.parse_args()

    with open(args.skymap[0], "rb") as f:
        skyMap = pickle.load(f)
        reportSkyMapInfo(skyMap)
        plotStyles[args.style](skyMap)
