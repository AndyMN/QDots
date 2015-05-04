from pylab import *


class Plotter:

    def __init__(self):
        self.title = None
        self.xLabel = None
        self.yLabel = None
        self.legend = []
        self.lineStyle = '-'
        self.nPlots = 0


    def plotEigenvalues(self,potWellSolver, kVec, nSmallest):
        EMatrix = potWellSolver.getEigenValues(kVec, nSmallest)
        EMatrix = np.array(EMatrix)
        if type(kVec) == int:
            for i in xrange(0, nSmallest):
                plot(kVec, EMatrix[0][i],'+')
        else:
            if self.nPlots > 0:
                self.lineStyle = '--'
            plot(kVec, EMatrix.T, self.lineStyle)
        self.setEigenvaluePlotLabels(potWellSolver, nSmallest)
        self.nPlots += 1

    def plotEigenvectors(self, potWellSolver, k=0,  state=0):
        xVec = potWellSolver.xAxisVector
        eigenVector = potWellSolver.getEigenvectors(k, state)
        vectorMatrix = zeros((potWellSolver.nGridPoints, potWellSolver.matrixDim))
        for i in xrange(0, potWellSolver.matrixDim):
            vectorMatrix[:,i] = np.squeeze(np.array(eigenVector[i*potWellSolver.nGridPoints:(i+1)*potWellSolver.nGridPoints].real))
        if self.nPlots > 0:
            self.lineStyle = '--'
        plot(xVec, vectorMatrix.real, self.lineStyle)
        self.setEigenvectorPlotLabels(potWellSolver)
        self.nPlots += 1

    def plotMixing(self, potWellSolver, kVec, state=0):
        fractionMatrix = zeros((potWellSolver.matrixDim/2, len(kVec)))
        column = 0
        for k in kVec:
            fractions = potWellSolver.getMixing(k, state)
            for i in xrange(potWellSolver.matrixDim/2):
                if i == 0:
                    fractionMatrix[i][column] = fractions[0] + fractions[3]
                elif i == 1:
                    fractionMatrix[i][column] = fractions[1] + fractions[2]
                elif i == 2:
                    fractionMatrix[i][column] = fractions[4] + fractions[5]
            column += 1
        if self.nPlots > 0:
            self.lineStyle = '--'
        plot(kVec, fractionMatrix.T, self.lineStyle)
        self.setMixingPlotLabels(potWellSolver)
        self.nPlots += 1

    def plotGridPointAnalysis(self, potWellSolver, gridPointVec, nSmallest=1, k=0,  state=0):
        eigenValues = zeros((nSmallest, len(gridPointVec)))
        startNGridPoints = potWellSolver.nGridPoints
        column = 0
        for gridPoint in gridPointVec:
            print gridPoint
            if gridPoint < 2000:
                potWellSolver.setDense(1)
            else:
                potWellSolver.setDense(0)
            potWellSolver.setGridPoints(gridPoint)
            w = potWellSolver.getEigenValues(k, nSmallest)
            for i in xrange(nSmallest):
                eigenValues[i][column] = w[i]
            column += 1
        potWellSolver.setGridPoints(startNGridPoints)
        if self.nPlots > 0:
            self.lineStyle = '--'
        plot(gridPointVec, eigenValues.T, self.lineStyle)
        self.setGridPointAnalysisPlotLabels(potWellSolver, nSmallest, k)
        self.nPlots += 1

    def plotRotatedMixing(self, potWellSolver, kVec):
        rotValues = zeros((2, len(kVec)))
        rotateTo = None
        if potWellSolver.potWell.nDirection == 1:
            rotateTo = "z"
        elif potWellSolver.potWell.nDirection == 3:
            rotateTo = "x"
        column = 0
        for k in kVec:
            rotVals = potWellSolver.rotateMixing(k, rotateTo)
            rotValues[0][column] = rotVals[0].real
            rotValues[1][column] = rotVals[1].real
            column += 1
        plot(kVec, rotValues.T, '*')
        self.setRotatedMixingPlotLabels(potWellSolver)

    def plotBulkEigenValues(self, potWellSolver, kVec):
        EigenvaluesMatrix = potWellSolver.getEigenValues(kVec, potWellSolver.matrixDim, True)
        plot(kVec, EigenvaluesMatrix.T)
        self.setBulkPlotLabels(potWellSolver)

    def setBulkPlotLabels(self, potWellSolver):
        self.xLabel = "k (1/cm)"
        self.yLabel = "E (meV)"
        if potWellSolver.matrixDim == 6:
            self.legend.append("SO")
            self.legend.append("LH")
            self.legend.append("HH")
        elif potWellSolver.matrixDim == 4:
            self.legend.append("LH")
            self.legend.append("HH")
        if self.title:
            self.title += " + Bulk E(k)"
        else:
            self.title = potWellSolver.compound.name + ": Bulk E(k)"


    def setRotatedMixingPlotLabels(self, potWellSolver):
        rotateTo = None
        if potWellSolver.potWell.nDirection == 1:
            rotateTo = "z"
        elif potWellSolver.potWell.nDirection == 3:
            rotateTo = "x"
        if not self.title:
            self.title = potWellSolver.compound.name+": Rotated mixing for well in "+potWellSolver.potWell.direction+" direction"
        else:
            self.title += " + Rotated mixing for well in "+potWellSolver.potWell.direction+" direction"
        self.xLabel = "k"+potWellSolver.potWell.direction+" (1/cm)"
        self.yLabel = "Fractions"
        self.legend.append("HHrot")
        self.legend.append("LHrot")




    def displayPlots(self):
        title(self.title)
        xlabel(self.xLabel)
        ylabel(self.yLabel)
        legend(self.legend)
        show()
        self.clearPlot()

    def savePlots(self, fileName):
        title(self.title)
        xlabel(self.xLabel)
        ylabel(self.yLabel)
        legend(self.legend)
        savefig(fileName)
        self.clearPlot()

    def clearPlot(self):
        self.title = None
        self.xLabel = None
        self.yLabel = None
        self.legend = []
        self.lineStyle = '-'
        self.nPlots = 0
        clf()

    def setGridPointAnalysisPlotLabels(self, potWellSolver, nSmallest, k):
        self.xLabel = "# Gridpoints"
        self.yLabel = "E(meV)"
        self.title = potWellSolver.compound.name+": E(N) k = "+str(k)+" depth = "+str(potWellSolver.potWell.depth)+" meV"
        for i in xrange(nSmallest):
            self.legend.append("E"+str(i+1))

    def setEigenvectorPlotLabels(self, potWellSolver):
        potWellDirection = potWellSolver.potWell.direction
        matrixType = None
        if potWellSolver.matrixDim == 4:
            matrixType = "4x4"
        elif potWellSolver.matrixDim == 6:
            matrixType = "6x6"
        if not self.title:
            self.title = potWellSolver.compound.name+": Wavefunctions: Well in "+potWellDirection+" direction ("+matrixType+")"
        else:
            self.title += " + Wavefunctions: Well in "+potWellDirection+" direction ("+matrixType+")"
        self.xLabel = potWellDirection
        self.yLabel = "Psi("+potWellDirection+")"

        if potWellSolver.matrixDim == 4:
            self.legend.append("HH1")
            self.legend.append("LH1")
            self.legend.append("LH2")
            self.legend.append("HH2")
        elif potWellSolver.matrixDim == 6:
            self.legend.append("HH1")
            self.legend.append("LH1")
            self.legend.append("LH2")
            self.legend.append("HH2")
            self.legend.append("SO1")
            self.legend.append("SO2")

    def setEigenvaluePlotLabels(self, potWellSolver, nSmallest):
            directionOfK = None
            matrixType = None
            if potWellSolver.matrixDim == 4:
                matrixType = "4x4"
            elif potWellSolver.matrixDim == 6:
                matrixType = "6x6"
            if potWellSolver.potWell.nDirection == 1:
                directionOfK = "z"
            elif potWellSolver.potWell.nDirection == 3:
                directionOfK = "x"
            if not self.title:
                self.title = potWellSolver.compound.name+": E(k"+directionOfK+") "+matrixType
            else:
                self.title += " + E(k" + directionOfK + ") "+matrixType
            for i in xrange(0, nSmallest):
                self.legend.append("E"+str(i+1)+" "+matrixType)

            self.yLabel = "E (meV)"
            self.xLabel = "k"+directionOfK+" (1/cm)"

    def setMixingPlotLabels(self, potWellSolver):
        directionOfK = None
        if potWellSolver.potWell.nDirection == 1:
            directionOfK = "z"
        elif potWellSolver.potWell.nDirection == 3:
            directionOfK = "x"
        matrixType = None
        if potWellSolver.matrixDim == 4:
            matrixType = "4x4"
        elif potWellSolver.matrixDim == 6:
            matrixType = "6x6"
        self.yLabel = "Fraction"
        self.xLabel = "k"+directionOfK+" (1/cm)"
        potWellDirection = potWellSolver.potWell.direction
        if not self.title:
            self.title = potWellSolver.compound.name+": Mixing: Well in "+potWellDirection+" direction ("+matrixType+")"
        else:
            self.title += " + Mixing: well in "+potWellDirection+" direction ("+matrixType+")"
        self.legend.append("HH"+potWellDirection)
        self.legend.append("LH"+potWellDirection)
        if potWellSolver.matrixDim == 6:
            self.legend.append("SO"+potWellDirection)
