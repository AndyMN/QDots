from operator import itemgetter
from pylab import *
massElectron = 5.6778*10**(-13)  # in meV.s^2/cm^2
hbar = 6.58211928*10**(-13)  # in meV.s
wellWidth = 100*10**(-8)  # in cm
pi = math.pi
unitE = (hbar**2*pi**2)/(2*massElectron*wellWidth**2)



def fillMe(M,fillM,i,j,n):
    for x in xrange(i*n,(i+1)*n):
        for y in xrange(j*n, (j+1)*n):
            M[x][y] = fillM[x-i*n][y-j*n]

def fillTotal(M,fillM,n):
    i = 0
    for row in fillM:
        j = 0
        for item in row:
            fillMe(M, item, i, j, n)
            j += 1
        i += 1

class PotWellSolver:

    def __init__(self, compound, potWell, matrixDim=4):
        self.compound = compound
        self.potWell = potWell
        self.matrixDim = matrixDim
        self.unitV = self.potWell.getDepth()/unitE
        self.unitDelta = self.compound.getDelta()/unitE
        self.nGridPoints = 350
        self.xMax = 3
        self.xMin = -3
        self.xAxisVector = np.linspace(self.xMin, self.xMax, self.nGridPoints)
        self.stepSize = (self.xMax-self.xMin)/float(len(self.xAxisVector))
        self.potWellBoundary1 = len(self.xAxisVector)/2. - np.floor(self.potWell.getWidth()/(2.*self.stepSize))
        self.potWellBoundary2 = len(self.xAxisVector)/2. + np.ceil(self.potWell.getWidth()/(2.*self.stepSize))


    def setParameters(self, nGridPoints, xMin=-3, xMax=3):
        self.nGridPoints = nGridPoints
        self.xMin = xMin
        self.xMax = xMax
        self.xAxisVector = np.linspace(self.xMin, self.xMax, self.nGridPoints)
        self.stepSize = (self.xMax-self.xMin)/float(len(self.xAxisVector))
        self.potWellBoundary1 = len(self.xAxisVector)/2 - np.floor(self.potWell.getWidth()/(2*self.stepSize))
        self.potWellBoundary2 = len(self.xAxisVector)/2 + np.ceil(self.potWell.getWidth()/(2*self.stepSize))


    def getXAxisVector(self):
        return self.xAxisVector

    def setGridPoints(self, nGridPoints):
        self.setParameters(nGridPoints)

    def solveMatrix(self, k1):

        ky = 0
        diagP = None
        subdiagP = None
        superdiagP = None
        diagQ = None
        subdiagQ = None
        superdiagQ = None
        diagS = None
        subdiagS = None
        superdiagS = None
        diagR = None
        subdiagR = None
        superdiagR = None

        HKL = np.zeros((self.matrixDim*self.nGridPoints, self.matrixDim*self.nGridPoints), dtype=complex)

        if self.potWell.getNDirection() == 1:
            kz = k1
            diagP = (self.compound.getY1()/pi**2)*((kz**2+ky**2) + 2./self.stepSize**2)
            subdiagP = -self.compound.getY1()/(self.stepSize**2*pi**2)
            superdiagP = subdiagP

            diagQ = (self.compound.getY2()/pi**2) * ((-2.*kz**2+ky**2)+(2/self.stepSize**2))
            subdiagQ = -self.compound.getY2()/(self.stepSize**2*pi**2)
            superdiagQ = subdiagQ

            diagR = (np.sqrt(3)*self.compound.getY2()/pi**2) * (-2/self.stepSize**2 + ky**2)
            subdiagR = (np.sqrt(3)/pi**2)*(self.compound.getY2()/self.stepSize**2 - self.compound.getY3()*ky/self.stepSize)
            superdiagR = (np.sqrt(3)/pi**2)*(self.compound.getY2()/self.stepSize**2 + self.compound.getY3()*ky/self.stepSize)

            diagS = -2.*self.compound.getY3()*np.sqrt(3)*1j*ky*kz/pi**2
            subdiagS = np.sqrt(3)*self.compound.getY3()*kz*1j/(pi**2*self.stepSize)
            superdiagS = -np.sqrt(3)*self.compound.getY3()*kz*1j/(pi**2*self.stepSize)

        
        w, v = eigh(HKL)
        return w*unitE, v

    def getEigenValues(self, kVec, nSmallest):
        if type(kVec) == int:
            eigenValues, eigenVectors = self.solveMatrix(kVec)
            eigenValues = sorted(eigenValues)
            EArray = np.zeros((1,nSmallest))
            for i in xrange(0, nSmallest):
                EArray[0][i] = eigenValues[i*2].real
            return EArray
        else:
            EMatrix = np.zeros((nSmallest, len(kVec)))
            column = 0
            for k in kVec:
                eigenValues, eigenVectors = self.solveMatrix(k)
                eigenValues = sorted(eigenValues)
                for i in xrange(0, nSmallest):
                    EMatrix[i][column] = eigenValues[i*2].real
                column += 1
            return EMatrix

    def getEigenvectors(self, k, state):
        w, v = self.solveMatrix(k)
        data = [(w[i], v[:,i]) for i in xrange(0, self.matrixDim*self.nGridPoints)]
        data = sorted(data, key=itemgetter(0))
        return data[state][1]

    def getMixing(self, k, state):
        w, v = self.solveMatrix(k)
        data = [(w[i], v[:,i]) for i in xrange(0, self.matrixDim*self.nGridPoints)]
        data = sorted(data, key=itemgetter(0))
        eigenVector = data[state][1]
        splitVectors = np.zeros((self.nGridPoints, self.matrixDim),dtype=complex)
        for i in xrange(0,self.matrixDim):
            splitVectors[:,i] = eigenVector[i*self.nGridPoints:(i+1)*self.nGridPoints]

        normSQ = np.zeros((1,self.matrixDim))
        for i in xrange(0,self.matrixDim):
            normSQ[0][i] = norm(splitVectors[:,i])**2
        totalDensity = sum(normSQ[0])
        fractions = []
        for i in xrange(0,self.matrixDim):
            fractions.append(normSQ[0][i]/totalDensity)
        return fractions


    def getStepSize(self):
        return self.stepSize

    def getGridPoints(self):
        return self.nGridPoints

    def getMatrixDim(self):
        return self.matrixDim

    def setXMax(self, xMax):
        self.setParameters(self.nGridPoints, self.xMin, xMax)

    def setXMin(self, xMin):
        self.setParameters(self.nGridPoints, xMin, self.xMax)

    def setXRange(self, xMin, xMax):
        self.setParameters(self.nGridPoints, xMin, xMax)

    def getXMax(self):
        return self.xMax

    def getXMin(self):
        return self.xMin





