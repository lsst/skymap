"""
To do:
- Tweak pixel scale so the average scale is as specified, rather than the scale at the center
"""
import math
import dodecahedron
import lsst.daf.base as dafBase
import lsst.afw.coord as afwCoord
import lsst.afw.geom as afwGeom
import lsst.afw.image as afwImage

# LSST plate scale is 50 um/arcsec
# LSST pixel size is 10 um
# Default sky pixel scale is 1/sqrt(2) of image pixel scale
_DefaultPlateScale = 10.0 / (50.0 * 3600.0 * math.sqrt(2.0))

class SkyTileInfo(object):
    def __init__(self, ind, ctrCoord, vertexCoordList, overlap, wcsFactory):
        """
        
        Inputs:
        - ind: index of sky tile
        - ctrCoord: Coord of center of sky tile
        - vertexCoordList: list of Coords of vertices that define edge of optimal area
        - overlap: minimum overlap with optimal area of adjacent sky tiles
        - wcsFactory: a WcsFactory object
        """
        
#         print "SkyTileInfo(ind=%s, ctrCoord=%s, overlap=%0.1f)" % \
#             (ind, ctrCoord.getPosition(afwCoord.DEGREES), overlap)
        self._ind = ind
        self._ctrCoord = ctrCoord
        self._vertexCoordList = vertexCoordList
        self._overlap = float(overlap)
        wcsFactory
        
        DebugMinInd = 0

        # start by computing everything relative to the center of the lower left pixel;
        # compute a WCS and use it to determine the required size of the tile
        basicWcs = wcsFactory.makeWcs(ctrInd=(0, 0), ctrCoord=self._ctrCoord)
        minPos = [0, 0]
        maxPos = [0, 0]
        
        # for now ignore overlap and just compute position at each vertex
        ctrSphPos = self._ctrCoord.getPosition(afwCoord.DEGREES)
        overlapDeg = self._overlap / RadPerDeg
        if ind < DebugMinInd:
            print "center position =", self._ctrCoord.getPosition(afwCoord.DEGREES)
#        overlapDeg = 0
        for vertexCoord in self._vertexCoordList:
            vertexDeg = vertexCoord.getPosition(afwCoord.DEGREES)
            if overlapDeg == 0:
                numToTest = 1
            else:
                numToTest = 4
            for i in range(numToTest):
                offAngle = 90.0 * i
                offCoord = vertexCoord.clone()
                offCoord.offset(offAngle * RadPerDeg, overlapDeg * RadPerDeg)
                pixPos = basicWcs.skyToPixel(offCoord)
                if ind < DebugMinInd:
                    print "vertexPos=%s, offPos=%s, xy=%s" % (vertexDeg, offCoord.getPosition(afwCoord.DEGREES), pixPos)
                for j, pixPosVal in enumerate(pixPos):
                    minPos[j] = min(minPos[j], pixPosVal)
                    maxPos[j] = max(maxPos[j], pixPosVal)
        minInd = [afwImage.positionToIndex(val) for val in minPos]
        maxInd = [afwImage.positionToIndex(val) for val in maxPos]
        if ind < DebugMinInd:
            print "minInd=%s, maxInd=%s" % (minInd, maxInd)
        
        self._dimensions = tuple(1 + maxInd[i] - minInd[i] for i in range(2))
        
        # now offset things so the lower left pixel is 0,0 and compute a new WCS using that
        self._ctrPixInd = tuple(-ind for ind in minInd)
        
        self._wcs = wcsFactory.makeWcs(ctrInd=self._ctrPixInd, ctrCoord=self._ctrCoord)

    def getCtrCoord(self):
        return self._ctrCoord

    def getCtrIndex(self):
        return self._ctrPixInd
    
    def getDimensions(self):
        return self._dimensions

    def getIndex(self):
        """Return sky tile index"""
        return self._ind
    
    def getWcs(self):
        return self._wcs


class WcsFactory(object):
    """A factory for creating Wcs objects for the sky tiles.
    """
    def __init__(self, pixelScale, projection):
        if len(projection) != 3:
            raise RuntimeError("projection=%r; must have length 3" % (projection,))
        self._pixelScale = float(pixelScale)
        self._projection = str(projection)
        self._ctypes = [("%-5s%3s" % (("RA", "DEC")[i], self._projection)).replace(" ", "-")
            for i in range(2)]

    def makeWcs(self, ctrInd, ctrCoord, **kargs):
        """Make a Wcs
        
        Inputs:
        - ctrInd: pixel index of center of WCS (LSST standard); used to compute CRPIX
        - ctrCoord: sky coordinate of center of WCS; used as CRVAL
        **kargs: FITS keyword arguments for WCS
        """
        ps = dafBase.PropertySet()
        crPix = [ind + 1.0 for ind in ctrInd]
        crVal = ctrCoord.getPosition(afwCoord.DEGREES)
        for i in range(2):
            ip1 = i + 1
            ps.add("CTYPE%1d" % (ip1,), self._ctypes[i])
            ps.add("CRPIX%1d" % (ip1,), crPix[i])
            ps.add("CRVAL%1d" % (ip1,), crVal[i])
        ps.add("RADECSYS", "ICRS")
        ps.add("EQUINOX", 2000)
        ps.add("CD1_1", -self._pixelScale)
        ps.add("CD2_1", 0.0)
        ps.add("CD1_2", 0.0)
        ps.add("CD2_2", self._pixelScale)
        return afwImage.makeWcs(ps)

RadPerDeg = math.pi / 180.0

class SkyMap(object):
    """Information about sky tiles
    """
    def __init__(self,
        overlap = 3.5 * RadPerDeg,
        pixelScale = _DefaultPlateScale,
        projection = "STG",
        withFacesOnPoles = False,
    ):
        """Inputs:
        - overlap: minimum overlap of each tile into the optimal region of the adjacent tiles (rad)
        - pixelScale: nominal pixel scale in degrees/pixel
        - projection: one of the FITS WCS projection codes, such as:
          - STG: stereographic projection
          - MOL: Molleweide's projection
          - TAG: tangent-plane projection
        - withFacesOnPoles: if True center a face on each pole, else put a vertex on each pole
        """
        self._overlap = float(overlap)
        self._pixelScale = float(pixelScale)
        self._projection = str(projection)
        self._dodecahedron = dodecahedron.Dodecahedron(withFacesOnPoles)
        self._skyTileInfoList = []
        self._wcsFactory = WcsFactory(self._pixelScale, self._projection)

        def coordFromVec(vec, defRA=None):
            if abs(vec[0]) < 1e-14 and abs(vec[1]) < 1e-14:
                if defRA == None:
                    raise RuntimeError("At pole and defRA==None")
                if vec[2] > 0:
                    dec = 90.0
                else:
                    dec = -90.0
                return afwCoord.makeCoord(afwCoord.ICRS, afwGeom.makePointD(defRA, dec), afwCoord.DEGREES)
            return afwCoord.makeCoord(afwCoord.ICRS, afwGeom.makePointD(*vec))
        
        for ind in range(12):
            faceVec = self._dodecahedron.getFace(ind)
            faceCoord = coordFromVec(faceVec)
            faceRA = faceCoord.getPosition(afwCoord.DEGREES)[0]
            vertexVecList = self._dodecahedron.getVertices(ind)
            
            self._skyTileInfoList.append(SkyTileInfo(
                ind = ind,
                ctrCoord = faceCoord,
                vertexCoordList = [coordFromVec(vec, defRA=faceRA) for vec in vertexVecList],
                overlap = self._overlap,
                wcsFactory = self._wcsFactory,
            ))

    def getSkyTileInd(self, coord):
        """Return the index of the sky tile whose optimal area includes the coord.
        
        If coord is on a boundary between two sky tiles then one of the two tiles will be returned
        without warning; it is explicitly not defined how the choice is made.
        """
        return self._dodecahedron.getFaceInd(coord.getVector())

    def getSkyTileInfo(self, ind):
        """Get information about a sky tile
        """
        return self._skyTileInfoList[ind]


if __name__ == "__main__":
    skyMap = SkyMap()
    skyInfo = skyMap._skyTileInfoList[0]
    wcs = skyInfo.getWcs()
    print "tile 0 has center ind of %s" % (skyInfo.getCtrIndex(),)
    print "tile 0 has center RA=%0.3f, Dec=%0.3f" % tuple(skyInfo.getCtrCoord().getPosition(afwCoord.DEGREES))
    def pixToSky(xPos, yPos):
        skyPos = wcs.pixelToSky(xPos, yPos).getPosition(afwCoord.DEGREES)
        print "xPos=%10.1f, yPos=%10.1f --> RA=%8.3f, Dec=%8.3f" % (xPos, yPos, skyPos[0], skyPos[1])
    def skyToPix(ra, dec):
        skyCoord = afwCoord.Coord(afwGeom.makePointD(ra, dec), afwCoord.DEGREES)
        pixPos = wcs.skyToPixel(skyCoord)
        print "xPos=%10.1f, yPos=%10.1f <-- RA=%8.3f, Dec=%8.3f" % (pixPos[0], pixPos[1], ra, dec)
        
    ctrInd = skyInfo.getCtrIndex()
    pixToSky(ctrInd[0], ctrInd[1])
    pixToSky(0, 0)
    skyToPix(0, 89)
    skyToPix(0, 90)
    skyToPix(180, 89)
    pixToSky(0,  900000)
    pixToSky(0,  980000)
    pixToSky(0, 1200000)
    for i in range(12):
        skyTileInfo = skyMap.getSkyTileInfo(i)
        dimensions = skyTileInfo.getDimensions()
        numPix = dimensions[0] * dimensions[1]
        wcs = skyTileInfo.getWcs()
        ctrCoord = skyTileInfo.getCtrCoord()
        print "sky tile %2d dimensions = %s; #pixels=%s" % (i, dimensions, numPix)
        for xInd in (0, dimensions[0]/2, dimensions[0]):
            xPos = afwImage.indexToPosition(xInd)
            for yInd in (0, dimensions[1]/2, dimensions[1]):
                yPos = afwImage.indexToPosition(yInd)
                edgeCoord = wcs.pixelToSky(xPos, yPos)
                angSep = ctrCoord.angularSeparation(edgeCoord, afwCoord.DEGREES)
                print "xInd=%8s, yInd=%8s, distFromCtr=%3.1f" % (xInd, yInd, angSep)
                
