"""
To do:
- Tweak pixel scale so the average scale is as specified, rather than the scale at the center
"""
import math
import lsst.afw.coord as afwCoord
import lsst.afw.geom as afwGeom
import detail
import skyTileInfo

# LSST plate scale is 50 um/arcsec
# LSST pixel size is 10 um
# Default sky pixel scale is 1/sqrt(2) of image pixel scale
_DefaultPlateScale = 10.0 / (50.0 * 3600.0 * math.sqrt(2.0))

_RadPerDeg = math.pi / 180.0

class SkyMap(object):
    """Information about sky tiles
    """
    def __init__(self,
        overlap = 3.5 * _RadPerDeg,
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
        self._dodecahedron = detail.Dodecahedron(withFacesOnPoles)
        self._skyTileInfoList = []
        self._wcsFactory = detail.WcsFactory(self._pixelScale, self._projection)

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
            
            self._skyTileInfoList.append(skyTileInfo.SkyTileInfo(
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
