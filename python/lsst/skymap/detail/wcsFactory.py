import lsst.daf.base as dafBase
import lsst.afw.coord as afwCoord
import lsst.afw.image as afwImage

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
