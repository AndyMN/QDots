from baseCompounds import *
from PotentialWell import PotentialWell
from PotWellSolver import PotWellSolver
from Plotter import Plotter
from pylab import *


potWellZ = PotentialWell("z")
potWellX = PotentialWell("x")

GaAs = Compound(GaAsValues)
BigPotWellSolverGaAsZ = PotWellSolver(GaAs, potWellZ, 6)
SmallPotWellSolverGaAsZ = PotWellSolver(GaAs, potWellZ, 4)
BigPotWellSolverGaAsX = PotWellSolver(GaAs, potWellX, 6)
SmallPotWellSolverGaAsX = PotWellSolver(GaAs, potWellX, 4)


Si = Compound(SiValues)
BigPotWellSolverSiZ = PotWellSolver(Si, potWellZ, 6)
SmallPotWellSolverSiZ = PotWellSolver(Si, potWellZ, 4)
BigPotWellSolverSiX = PotWellSolver(Si, potWellX, 6)
SmallPotWellSolverSiX = PotWellSolver(Si, potWellX, 4)

myPlotter = Plotter()

kVec = arange(-3, 3, 0.1)
gridPoints = arange(100, 5600, 500)

myPlotter.plotEigenvalues(BigPotWellSolverGaAsX, kVec, 6)
myPlotter.savePlots("GaAsBigXE.png")
myPlotter.plotEigenvalues(BigPotWellSolverGaAsZ, kVec, 6)
myPlotter.savePlots("GaAsBigZE.png")
myPlotter.plotEigenvalues(SmallPotWellSolverGaAsX, kVec, 4)
myPlotter.savePlots("GaAsSmallXE.png")
myPlotter.plotEigenvalues(SmallPotWellSolverGaAsZ, kVec, 4)
myPlotter.savePlots("GaAsSmallZE.png")
myPlotter.plotEigenvalues(BigPotWellSolverGaAsZ, kVec, 6)
myPlotter.plotEigenvalues(SmallPotWellSolverGaAsZ, kVec, 4)
myPlotter.savePlots("GaAsBigANDSmallE.png")
myPlotter.plotMixing(BigPotWellSolverGaAsZ, kVec)
myPlotter.savePlots("GaAsMixingZ.png")
myPlotter.plotMixing(BigPotWellSolverGaAsX, kVec)
myPlotter.savePlots("GaAsMixingX.png")
myPlotter.plotEigenvectors(BigPotWellSolverGaAsX)
myPlotter.savePlots("GaAsWaveX.png")
myPlotter.plotEigenvectors(BigPotWellSolverGaAsZ)
myPlotter.savePlots("GaAsWaveZ.png")
myPlotter.plotRotatedMixing(SmallPotWellSolverGaAsZ, kVec)
myPlotter.plotMixing(SmallPotWellSolverGaAsX, kVec)
myPlotter.savePlots("GaAsMixingRotationOfZ.png")
myPlotter.plotRotatedMixing(SmallPotWellSolverGaAsX, kVec)
myPlotter.plotMixing(SmallPotWellSolverGaAsZ, kVec)
myPlotter.savePlots("GaAsMixingRotationOfX.png")



myPlotter.plotEigenvalues(BigPotWellSolverSiX, kVec, 6)
myPlotter.savePlots("SiBigXE.png")
myPlotter.plotEigenvalues(BigPotWellSolverSiZ, kVec, 6)
myPlotter.savePlots("SiBigZE.png")
myPlotter.plotEigenvalues(SmallPotWellSolverSiX, kVec, 4)
myPlotter.savePlots("SiSmallXE.png")
myPlotter.plotEigenvalues(SmallPotWellSolverSiZ, kVec, 4)
myPlotter.savePlots("SiSmallZE.png")
myPlotter.plotEigenvalues(BigPotWellSolverSiZ, kVec, 6)
myPlotter.plotEigenvalues(SmallPotWellSolverSiZ,  kVec, 4)
myPlotter.savePlots("SiBigANDSmallE.png")
myPlotter.plotMixing(BigPotWellSolverSiZ, kVec)
myPlotter.savePlots("SiMixingZ.png")
myPlotter.plotMixing(BigPotWellSolverSiX, kVec)
myPlotter.savePlots("SiMixingX.png")
myPlotter.plotEigenvectors(BigPotWellSolverSiX)
myPlotter.savePlots("SiWaveX.png")
myPlotter.plotEigenvectors(BigPotWellSolverSiZ)
myPlotter.savePlots("SiWaveZ.png")
myPlotter.plotRotatedMixing(SmallPotWellSolverSiZ, kVec)
myPlotter.plotMixing(SmallPotWellSolverSiX, kVec)
myPlotter.savePlots("SiMixingRotationOfZ.png")
myPlotter.plotRotatedMixing(SmallPotWellSolverSiX, kVec)
myPlotter.plotMixing(SmallPotWellSolverSiZ, kVec)
myPlotter.savePlots("SiMixingRotationOfX.png")


kVec = arange(-10**7, 10**7, 10**5)

myPlotter.plotBulkEigenValues(BigPotWellSolverSiZ, kVec)
myPlotter.savePlots("SiBulkBig.png")
myPlotter.plotBulkEigenValues(SmallPotWellSolverSiZ, kVec)
myPlotter.savePlots("SiBulkSmall.png")

myPlotter.plotBulkEigenValues(BigPotWellSolverGaAsZ, kVec)
myPlotter.savePlots("GaAsBulkBig.png")
myPlotter.plotBulkEigenValues(SmallPotWellSolverGaAsX, kVec)
myPlotter.savePlots("GaAsBulkSmall.png")