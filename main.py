from baseCompounds import *
from PotentialWell import PotentialWell
from PotWellSolver import PotWellSolver
from Plotter import Plotter
from pylab import *

nGridPoints = 10

potWellZ = PotentialWell("z")
potWellX = PotentialWell("x")

GaAs = Compound(GaAsValues)
SmallPotWellSolverGaAsZ = PotWellSolver(GaAs, potWellZ, 4)
SmallPotWellSolverGaAsX = PotWellSolver(GaAs, potWellX, 4)

SmallPotWellSolverGaAsX.setGridPoints(nGridPoints)
SmallPotWellSolverGaAsZ.setGridPoints(nGridPoints)

nGridPointsZ = SmallPotWellSolverGaAsZ.nGridPoints
nGridPointsX = SmallPotWellSolverGaAsX.nGridPoints
matrixDim = 4



myPlotter = Plotter()


kVec = arange(-3,3,0.1)

fractionsHHXinZ = []
fractionsHH1 = []
fractionsHH2 = []
fractionsLHinZ = []
fractionsLH1 = []
fractionsLH2 = []
for k in kVec:

    EigenvectorsZ = SmallPotWellSolverGaAsZ.getEigenvectors(k)
    EigenvectorsX = SmallPotWellSolverGaAsX.getEigenvectors(k)

    splitVectorsZ = np.zeros((nGridPointsZ, matrixDim), dtype=complex)
    for i in xrange(matrixDim):
        splitVectorsZ[:,i] = np.squeeze(np.array(EigenvectorsZ[i*nGridPointsZ:(i+1)*nGridPointsZ]))

    splitVectorsX = np.zeros((nGridPointsX, matrixDim), dtype=complex)
    for i in xrange(matrixDim):
        splitVectorsX[:,i] = np.squeeze(np.array(EigenvectorsX[i*nGridPointsX:(i+1)*nGridPointsX]))


    H1rot = np.sqrt(2)/4 * (splitVectorsX[:,0] + splitVectorsX[:,3]) + np.sqrt(6)/4 * (splitVectorsX[:,1] + splitVectorsX[:,2])
    L1rot = np.sqrt(6)/4 * (splitVectorsX[:,3] - splitVectorsX[:,0]) + np.sqrt(2)/4 * (-splitVectorsX[:,1] - splitVectorsX[:,2])
    L2rot = np.sqrt(6)/4 * (splitVectorsX[:,0] + splitVectorsX[:,3]) + np.sqrt(2)/4 * (-splitVectorsX[:,1] - splitVectorsX[:,2])
    H2rot = np.sqrt(2)/4 * (splitVectorsX[:,3] - splitVectorsX[:,0]) + np.sqrt(6)/4 * (splitVectorsX[:,1] - splitVectorsX[:,2])

    H1normsq = norm(H1rot)**2
    H2normsq = norm(H2rot)**2
    L1normsq = norm(L1rot)**2
    L2normsq = norm(L2rot)**2

    total = H1normsq+ H2normsq + L1normsq + L2normsq

    HHnormsq = H1normsq + H2normsq
    LHnormsq = L1normsq + L2normsq

    fractionsHH1.append(H1normsq/total)
    fractionsHH2.append(H2normsq/total)
    fractionsLH1.append(LHnormsq/total)
    fractionsHHXinZ.append(HHnormsq/total)


plot(kVec, fractionsLH1, '*')
plot(kVec, fractionsHHXinZ, '*')
myPlotter.plotMixing(SmallPotWellSolverGaAsX, kVec)
show()

