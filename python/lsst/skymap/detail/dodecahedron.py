#!/usr/bin/env python
import math
import numpy

class Dodecahedron(object):
    """A dodecahedron
    
    Contains positions of faces and associated vertices
    """
    def __init__(self, withFacesOnPoles=False):
        """Construct a Dodecahedron
        
        @param[in] withFacesOnPoles: if True center a face on each pole, else put a vertex on each pole
        """
        # Basis cartesian vectors describing the faces of a dodecahedron; the full set of vectors is obtained
        # by choosing both signs of each nonzero component of each vector.
        # The orientation of the resulting dodecahedron, while very convenient
        # for specifying face vectors, is not an orientation we want so it must be rotated.
        g = (1.0 + math.sqrt(5.0)) / 2.0
        faceBases = (
            (0, 1, g),
            (1, g, 0),
            (g, 0, 1),
        )
        unrotFaceVecList = _computeFullVecList(faceBases)
        unrotVertexVecList = _computeDodecahedronVertices(unrotFaceVecList)
        
        if withFacesOnPoles:
            # one face is centered on each pole
            vec0, vec1 = _findClosePair(unrotFaceVecList, 0)
            rotMat = _computeCoordTransform(vec0, vec1)
        else:
            # one vertex is on each pole
            vec0, vec1 = _findClosePair(unrotVertexVecList, 0)
            rotMat = _computeCoordTransform(vec0, vec1, vec1NegativeX=True)
        self.vertexVecList = [numpy.dot(rotMat, unrotVertexVec) for unrotVertexVec in unrotVertexVecList]
        unsortedFaceList = [numpy.dot(rotMat, unrotFaceVec) for unrotFaceVec in unrotFaceVecList]
        self.faceVecList = _sortedVectorList(unsortedFaceList)

    def getFaceVectorList(self):
        """Return a list of unit vectors of the center of each face
        """
        return self.faceVecList[:]
    
    def getFace(self, ind):
        """Return the position vector of the center of the specified face
        """
        return self.faceVecList[ind][:]

    def getVertices(self, ind):
        """Return the vertices for a given face, as a list of unit vectors
        """
        faceVec = self.getFace(ind)
        return _findCloseList(self.vertexVecList, faceVec)

    def getFaceInd(self, vec):
        """Return the index of the face containing the cartesian vector
        
        If the vector is on a border, picks one of the two faces in an undocumented way.
        """
        return numpy.argmax(numpy.sum(self.faceVecList * vec, 1))

def computeRotationMatrix(angle, axis):
    """Return a 3D rotation matrix for rotation by a specified amount around a specified axis
    
    Inputs:
    - angle: amount of rotation (rad)
    - axis: axis of rotation; one of 0, 1 or 2 for x, y or z
    """
    cosAng = math.cos(angle)
    sinAng = math.sin(angle)
    rotMat = numpy.zeros((3,3), dtype=float)
    rotMat[axis, axis] = 1
    rotMat[(axis + 1) % 3, (axis + 1) % 3] = cosAng
    rotMat[(axis + 2) % 3, (axis + 1) % 3] = sinAng
    rotMat[(axis + 1) % 3, (axis + 2) % 3] = -sinAng
    rotMat[(axis + 2) % 3, (axis + 2) % 3] = cosAng
    return rotMat

def _computeCoordTransform(vec0, vec1, vec1NegativeX=False):
    """Compute a rotation matrix that puts vec0 along z and vec1 along +x in the xz plane
    
    Inputs:
    - vec0: vector 0
    - vec1: vector 1
    - vec1NegativeX: if True then vec1 is rotated to face negative x
    """
    # rotate around x by angle of vec0 from z to y
    xAng = math.atan2(vec0[1], vec0[2])
    xRotMat = computeRotationMatrix(xAng, 0)

    # rotate around y by -angle of rotated vec0 from z to x
    vec0RotX = numpy.dot(xRotMat, vec0)
    yAng = -math.atan2(vec0RotX[0], vec0RotX[2])
    yRotMat = computeRotationMatrix(yAng, 1)
    xyRotMat = numpy.dot(yRotMat, xRotMat)
    
    # rotate around z by -angle of rotated vec1 from +/-x to y
    vec1RotXY = numpy.dot(xyRotMat, vec1)
    xVal = vec1RotXY[0]
    if vec1NegativeX:
        xVal = -xVal
    zAng = -math.atan2(vec1RotXY[1], xVal)
    zRotMat = computeRotationMatrix(zAng, 2)
    xyzRotMat = numpy.dot(zRotMat, xyRotMat)
    return xyzRotMat

def _computeDodecahedronVertices(faceVecList):
    """Given a vector of face positions of a Dodecahedron compute the vertices
    """
    closeIndSetList = []
    vertexDict = {}
    for i in range(len(faceVecList)):
        closeIndSet = _findCloseIndexSet(faceVecList, i)
        if len(closeIndSet) != 5:
            raise RuntimeError("Found %s vertices instead of 5 near %s: %s" % \
                (len(closeIndSet), faceVecList[i], closeIndSet))
        closeIndSetList.append(closeIndSet)
    for i, iCloseIndSet in enumerate(closeIndSetList):
        for j in iCloseIndSet:
            jCloseIndSet = closeIndSetList[j]
            sharedCloseIndSet = iCloseIndSet.intersection(jCloseIndSet)
            if len(sharedCloseIndSet) != 2:
                raise RuntimeError("Found %s vertices instead of 2 near %s and %s: %s" % \
                    (len(sharedCloseIndSet), faceVecList[i], faceVecList[j], sharedCloseIndSet))
            for k in sharedCloseIndSet:
                key = frozenset((i, j, k))
                if key in vertexDict:
                    continue
                vertexVec = faceVecList[i] + faceVecList[j] + faceVecList[k]
                vertexVec /= numpy.sqrt(numpy.sum(vertexVec**2))
                vertexDict[key] = vertexVec
    return vertexDict.values()

def _computeFullVecList(basisSet):
    """Given a collection of basis vectors, compute all permutations with both signs of all nonzero values
    
    For example: [(0, 1, 2)] -> [(0, 1, 2), (0, -1, 2), (0, 1, -2), (0, -1, -2)]
    """
    fullSet = []
    for basisVec in basisSet:
        vecLen = math.sqrt(numpy.sum(numpy.array(basisVec)**2))
        valueList = []
        for basisValue in basisVec:
            if basisValue == 0:
                valueList.append((0,))
            else:
                valueList.append((basisValue, -basisValue))
        fullSet += list(numpy.array((x, y, z))/vecLen
                        for z in valueList[2]
                        for y in valueList[1]
                        for x in valueList[0]
                        )
    return fullSet

def _findCloseIndexSet(vecList, ind):
    """Given a list of cartesian vectors, return a set of indices of those closes to one of them
    
    This is intended for regular grids where distances are quantized
    
    Inputs:
    - vecList: list of cartesian vectors
    - ind: index of vector to be nears
    """
    dotProductList = numpy.round(numpy.sum(vecList * vecList[ind], 1), 2)
    dotProductList[ind] = -9e99
    minDist = numpy.max(dotProductList)
    indList = numpy.arange(len(dotProductList))[dotProductList==minDist]
    return set(indList)

def _findCloseList(vecList, vec):
    """Given a list of cartesian vectors, return all those closest to a specified position
    
    This is intended for regular grids where distances are quantized
    
    Inputs:
    - vecList: list of cartesian vectors
    - vec: vector to be near
    """
    dotProductList = numpy.round(numpy.sum(vecList * vec, 1), 2)
    minDist = numpy.max(dotProductList)
    indList = numpy.arange(len(dotProductList))[dotProductList==minDist]
    retList = numpy.take(vecList, indList, 0)
    return retList

def _findClosePair(vecList, ind=0):
    """Given a list of cartesian vectors and an index, return the vector and one of its closest neighbors
    
    Inputs:
    - vecList: list of cartesian vectors
    - ind: index of first vector
    """
    vec = vecList[ind]
    otherVecList = vecList[0:ind] + vecList[ind+1:]
    ind1 = numpy.argmax(numpy.sum(otherVecList * vec, 1))
    return vec, otherVecList[ind1]
    
def _sortedVectorList(vecList):
    """Return a list of cartesian vectors sorted by decreasing latitude and increasing longitude
    """
    def vecToSort(vec):
        ang = round(math.atan2(vec[1], vec[0]), 2)
        if ang < 0:
            ang += 2.0 * math.pi
        return (-round(vec[2], 1), ang, vec)

    decoratedList = [vecToSort(v) for v in vecList]
    decoratedList.sort()
    return [d[2] for d in decoratedList]

if __name__ == "__main__":
    numpy.set_printoptions(precision=2, suppress=True, linewidth=120)
    import lsst.afw.coord as afwCoord
    import lsst.afw.geom as afwGeom

    print "Dodecahedron with vertices on poles"
    vertexDodec = Dodecahedron(withFacesOnPoles=False)
    for i in range(12):
        faceVec = vertexDodec.getFace(i)
        print "Face %2d: %s" % (i, faceVec)
