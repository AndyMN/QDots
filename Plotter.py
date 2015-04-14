from pylab import *


class Plotter:

    def __init__(self):
        self.title = None
        self.xLabel = None
        self.yLabel = None
        self.legend = []


    def plotEigenvalues(self,potWellSolver, kVec, nSmallest):
        EMatrix = potWellSolver.getEigenValues(kVec, nSmallest)
        EMatrix = np.array(EMatrix)
        if type(kVec) == int:
            for i in xrange(0, nSmallest):
                plot(kVec, EMatrix[0][i],'+')
        else:
            plot(kVec, EMatrix.T)
        self.setEigenvaluePlotLabels(potWellSolver, nSmallest)

    def plotEigenvectors(self, potWellSolver, k=0,  state=0):
        xVec = potWellSolver.getXAxisVector()
        eigenVector = potWellSolver.getEigenvectors(k, state)
        vectorMatrix = zeros((potWellSolver.getGridPoints(), potWellSolver.getMatrixDim()))
        for i in xrange(0, potWellSolver.getMatrixDim()):
            vectorMatrix[:,i] = np.squeeze(np.array(eigenVector[i*potWellSolver.getGridPoints():(i+1)*potWellSolver.getGridPoints()].real))
        plot(xVec, vectorMatrix.real)
        self.setEigenvectorPlotLabels(potWellSolver)

    def plotMixing(self, potWellSolver, kVec, state=0):
        fractionMatrix = zeros((potWellSolver.getMatrixDim()/2, len(kVec)))
        column = 0
        for k in kVec:
            fractions = potWellSolver.getMixing(k, state)
            for i in xrange(potWellSolver.getMatrixDim()/2):
                if i == 0:
                    fractionMatrix[i][column] = fractions[0] + fractions[3]
                elif i == 1:
                    fractionMatrix[i][column] = fractions[1] + fractions[2]
                elif i == 2:
                    fractionMatrix[i][column] = fractions[4] + fractions[5]
            column += 1
        plot(kVec, fractionMatrix.T)
        self.setMixingPlotLabels(potWellSolver)

    def setTitle(self, title):
        self.title = title

    def setYLabel(self, yLabel):
        self.yLabel = yLabel

    def setXLabel(self, xLabel):
        self.xLabel = xLabel


    def displayPlots(self):
        title(self.title)
        xlabel(self.xLabel)
        ylabel(self.yLabel)
        legend(self.legend)
        show()

    def setEigenvectorPlotLabels(self, potWellSolver):
        potWellDirection = potWellSolver.potWell.getDirection()
        matrixType = None
        if potWellSolver.getMatrixDim() == 4:
            matrixType = "4x4"
        elif potWellSolver.getMatrixDim() == 6:
            matrixType = "6x6"
        if not self.title:
            self.title = potWellSolver.compound.getName()+": Golffuncties: put in "+potWellDirection+" richting ("+matrixType+")"
        else:
            self.title += " + Golffuncties: put in "+potWellDirection+" richting ("+matrixType+")"
        self.xLabel = potWellDirection
        self.yLabel = "Psi("+potWellDirection+")"

        if potWellSolver.getMatrixDim() == 4:
            self.legend.append("HH1")
            self.legend.append("LH1")
            self.legend.append("LH2")
            self.legend.append("HH2")
        elif potWellSolver.getMatrixDim() == 6:
            self.legend.append("HH1")
            self.legend.append("LH1")
            self.legend.append("LH2")
            self.legend.append("HH2")
            self.legend.append("SO1")
            self.legend.append("SO2")

    def setEigenvaluePlotLabels(self, potWellSolver, nSmallest):
            directionOfK = None
            matrixType = None
            if potWellSolver.getMatrixDim() == 4:
                matrixType = "4x4"
            elif potWellSolver.getMatrixDim() == 6:
                matrixType = "6x6"
            if potWellSolver.potWell.getNDirection() == 1:
                directionOfK = "z"
            elif potWellSolver.potWell.getNDirection() == 3:
                directionOfK = "x"
            if not self.title:
                self.title = potWellSolver.compound.getName()+": E(k"+directionOfK+") "+matrixType
            else:
                self.title += " + E(k" + directionOfK + ") "+matrixType
            for i in xrange(0, nSmallest):
                self.legend.append("E"+str(i+1)+" "+matrixType)

            self.yLabel = "E (meV)"
            self.xLabel = "k"+directionOfK+" (1/cm)"

    def setMixingPlotLabels(self, potWellSolver):
        directionOfK = None
        if potWellSolver.potWell.getNDirection() == 1:
            directionOfK = "z"
        elif potWellSolver.potWell.getNDirection() == 3:
            directionOfK = "x"
        matrixType = None
        if potWellSolver.getMatrixDim() == 4:
            matrixType = "4x4"
        elif potWellSolver.getMatrixDim() == 6:
            matrixType = "6x6"
        self.yLabel = "Fraction"
        self.xLabel = "k"+directionOfK+" (1/cm)"
        potWellDirection = potWellSolver.potWell.getDirection()
        if not self.title:
            self.title = potWellSolver.compound.getName()+": Golffuncties Mixing: put in "+potWellDirection+" richting ("+matrixType+")"
        else:
            self.title += " + Golffuncties Mixing: put in "+potWellDirection+" richting ("+matrixType+")"
        self.legend.append("HH"+potWellDirection)
        self.legend.append("LH"+potWellDirection)
        if potWellSolver.getMatrixDim() == 6:
            self.legend.append("SO"+potWellDirection)
