#!/usr/bin/env python
import lsst.afw.coord as afwCoord
import lsst.afw.geom as afwGeom
import lsst.afw.image as afwImage
import lsst.skymap

skyMap = lsst.skymap.SkyMap()
for i in range(12):
    skyTileInfo = skyMap.getSkyTileInfo(i)
    ctrCoord = skyTileInfo.getCtrCoord()
    ctrPos = ctrCoord.getPosition(afwCoord.DEGREES)
    dimensions = skyTileInfo.getDimensions()
    numPix = dimensions[0] * dimensions[1]
    wcs = skyTileInfo.getWcs()
    leftCoord = wcs.pixelToSky(0.0, dimensions[1]/2.0)
    rightCoord = wcs.pixelToSky(dimensions[0]-1, dimensions[1]/2.0)
    topCoord = wcs.pixelToSky(dimensions[0]/2.0, dimensions[1] - 1)
    bottomCoord = wcs.pixelToSky(dimensions[0]/2.0, 0)
    xSpan = leftCoord.angularSeparation(rightCoord, afwCoord.DEGREES)
    ySpan = bottomCoord.angularSeparation(topCoord, afwCoord.DEGREES)
    print "sky tile %2d center RA/Dec = %6.1f, %6.1f; dimensions = %.1e x %.1e pixels (%0.1e total); span = %.1f x %.1f deg" % \
        (i, ctrPos[0], ctrPos[1], dimensions[0], dimensions[1], numPix, xSpan, ySpan)
