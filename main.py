from baseCompounds import *
from PotentialWell import PotentialWell
from PotWellSolver import PotWellSolver
from Plotter import Plotter
from pylab import *



GaAs = Compound(GaAsValues)
myPotWell = PotentialWell("x", 130, 1)
myPotWellSolver = PotWellSolver(GaAs, myPotWell, 6)
myPotWellSolver.setGridPoints(100)
myPlotter = Plotter()



myPlotter.plotMixing(myPotWellSolver, arange(0,1,0.1))
myPlotter.displayPlots()










