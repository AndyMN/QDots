from operator import itemgetter
from pylab import *
massElectron = 5.6778*10**(-13)  # in meV.s^2/cm^2
hbar = 6.58211928*10**(-13)  # in meV.s
wellWidth = 100*10**(-8)  # in cm
pi = math.pi
unitE = (hbar**2*pi**2)/(2*massElectron*wellWidth**2)


class PotWellSolver:

    def __init__(self, compound, potWell, matrixDim=4):
        self.compound = compound
        self.potWell = potWell
        self.matrixDim = matrixDim
        self.unitV = self.potWell.getDepth()/unitE
        self.unitDelta = self.compound.getDelta()/unitE
        self.nGridPoints = 100
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
        self.potWellBoundary1 = len(self.xAxisVector)/2. - np.floor(self.potWell.getWidth()/(2.*self.stepSize))
        self.potWellBoundary2 = len(self.xAxisVector)/2. + np.ceil(self.potWell.getWidth()/(2.*self.stepSize))


    def getXAxisVector(self):
        return self.xAxisVector

    def setGridPoints(self, nGridPoints):
        self.setParameters(nGridPoints)

    def makeMatrix(self, k):
        ky = 0

        potVec = zeros(self.nGridPoints)
        potVec[0:self.potWellBoundary1] = self.unitV
        potVec[self.potWellBoundary2:self.nGridPoints] = self.unitV
        V = diag(potVec)
        P = zeros((self.nGridPoints, self.nGridPoints), dtype=complex)
        Q = zeros((self.nGridPoints, self.nGridPoints), dtype=complex)
        R = zeros((self.nGridPoints, self.nGridPoints), dtype=complex)
        S = zeros((self.nGridPoints, self.nGridPoints), dtype=complex)
        diagP = 0
        subdiagP = 0
        superdiagP = 0
        diagQ = 0
        subdiagQ = 0
        superdiagQ = 0
        diagS = 0
        subdiagS = 0
        superdiagS = 0
        diagR = 0
        subdiagR = 0
        superdiagR = 0

        if self.potWell.getNDirection() == 1:
            diagP = (self.compound.getY1()/pi**2)*(k**2 + ky**2 + 2./self.stepSize**2)
            subdiagP = -self.compound.getY1()/(self.stepSize**2*pi**2)
            superdiagP = subdiagP

            diagQ = (self.compound.getY2()/pi**2) * ((ky**2-2*k**2)+(2/self.stepSize**2))
            subdiagQ = -self.compound.getY2()/(self.stepSize**2*pi**2)
            superdiagQ = subdiagQ

            diagR = (np.sqrt(3)*self.compound.getY2()/pi**2) * (-2/self.stepSize**2 + ky**2)
            subdiagR = (np.sqrt(3)/pi**2)*(self.compound.getY2()/self.stepSize**2 - self.compound.getY3()*ky/self.stepSize)
            superdiagR = (np.sqrt(3)/pi**2)*(self.compound.getY2()/self.stepSize**2 + self.compound.getY3()*ky/self.stepSize)


            diagS = -2.*self.compound.getY3()*np.sqrt(3)*1j*ky*k/pi**2
            subdiagS = np.sqrt(3)*self.compound.getY3()*k*1j/(pi**2*self.stepSize)
            superdiagS = -subdiagS

        elif self.potWell.getNDirection() == 3:
            diagP = (self.compound.getY1()/pi**2)*(k**2 + ky**2 + 2./self.stepSize**2)
            subdiagP = -self.compound.getY1()/(self.stepSize**2*pi**2)
            superdiagP = subdiagP

            diagQ = (self.compound.getY2()/pi**2) * (k**2 + ky**2 - (4/self.stepSize**2))
            subdiagQ = 2*self.compound.getY2()/(self.stepSize**2*pi**2)
            superdiagQ = subdiagQ

            diagR = (1/pi**2) * (-np.sqrt(3)*self.compound.getY2()*(k**2 - ky**2) + 1j*2*np.sqrt(3)*self.compound.getY3()*k*ky)

            subdiagS = (1j*self.compound.getY3()*np.sqrt(3))/(pi**2*self.stepSize) * (k-1j*ky)
            superdiagS = -subdiagS




        i, j = indices(P.shape)
        P[i == j] = diagP
        P[i == j-1] = superdiagP
        P[i == j+1] = subdiagP

        i, j = indices(Q.shape)
        Q[i == j] = diagQ
        Q[i == j-1] = superdiagQ
        Q[i == j+1] = subdiagQ

        i, j = indices(R.shape)
        R[i == j] = diagR
        R[i == j-1] = superdiagR
        R[i == j+1] = subdiagR

        i, j = indices(S.shape)
        S[i == j] = diagS
        S[i == j-1] = superdiagS
        S[i == j+1] = subdiagS

        HKL = None
        if self.matrixDim == 6:
            Delta = diag((ones(self.nGridPoints)*self.unitDelta))
            HKL = np.bmat([[P+Q+V, -S, R, zeros((self.nGridPoints, self.nGridPoints)), -S/np.sqrt(2), np.sqrt(2)*R],
                            [-S.conj().T, P-Q+V, zeros((self.nGridPoints, self.nGridPoints)), R, -np.sqrt(2)*Q, np.sqrt(3./2.)*S],
                            [R.conj().T, zeros((self.nGridPoints, self.nGridPoints)), P-Q+V, S, np.sqrt(3./2.)*S.conj().T, np.sqrt(2)*Q],
                            [zeros((self.nGridPoints, self.nGridPoints)), R.conj().T, S.conj().T, P+Q+V, -np.sqrt(2)*R.conj().T, -S.conj().T/np.sqrt(2)],
                            [-S.conj().T/np.sqrt(2), -np.sqrt(2)*Q.conj().T, np.sqrt(3./2.)*S, -np.sqrt(2)*R, P+Delta+V, zeros((self.nGridPoints, self.nGridPoints))],
                            [np.sqrt(2)*R.conj().T, np.sqrt(3./2.)*S.conj().T, np.sqrt(2)*Q.conj().T, -S/np.sqrt(2), zeros((self.nGridPoints, self.nGridPoints)), P+Delta+V]])
        elif self.matrixDim == 4:
            HKL = np.bmat([[P+Q+V, -S, R, zeros((self.nGridPoints, self.nGridPoints))],
                            [-S.conj().T, P-Q+V, zeros((self.nGridPoints, self.nGridPoints)), R],
                            [R.conj().T, zeros((self.nGridPoints, self.nGridPoints)), P-Q+V, S],
                            [zeros((self.nGridPoints, self.nGridPoints)), R.conj().T, S.conj().T, P+Q+V]])
        return HKL


    def calcEigs(self, k):
        HKL = self.makeMatrix(k)
        w, v = eigh(HKL)
        return w*unitE, v

    def calcEigVals(self, k):
        HKL = self.makeMatrix(k)
        w = eigvalsh(HKL)
        return w*unitE

    def getEigenValues(self, kVec, nSmallest):
        if type(kVec) == int:
            eigenValues = self.calcEigVals(kVec)
            eigenValues = sorted(eigenValues)
            EArray = np.zeros((1,nSmallest))
            for i in xrange(0, nSmallest):
                EArray[0][i] = eigenValues[i*2].real
            return EArray
        else:
            EMatrix = np.zeros((nSmallest, len(kVec)))
            column = 0
            for k in kVec:
                eigenValues = self.calcEigVals(k)
                eigenValues = sorted(eigenValues)
                for i in xrange(0, nSmallest):
                    EMatrix[i][column] = eigenValues[i*2].real
                column += 1
            return EMatrix

    def getEigenvectors(self, k, state):
        w, v = self.calcEigs(k)
        data = [(w[i], v[:,i]) for i in xrange(0, self.matrixDim*self.nGridPoints)]
        data = sorted(data, key=itemgetter(0))
        return data[state][1]

    def getMixing(self, k, state):
        w, v = self.calcEigs(k)
        data = [(w[i], v[:,i]) for i in xrange(self.matrixDim*self.nGridPoints)]
        data = sorted(data, key=itemgetter(0))
        eigenVector = data[state][1]
        splitVectors = np.zeros((self.nGridPoints, self.matrixDim), dtype=complex)
        for i in xrange(self.matrixDim):
            splitVectors[:,i] = np.squeeze(np.array(eigenVector[i*self.nGridPoints:(i+1)*self.nGridPoints]))

        normSQ = np.zeros(self.matrixDim)
        for i in xrange(self.matrixDim):
            normSQ[i] = norm(splitVectors[:,i])**2

        totalDensity = sum(normSQ)

        fractions = []
        for i in xrange(self.matrixDim):
            fractions.append(normSQ[i]/totalDensity)
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





