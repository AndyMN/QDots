from baseCompounds import *
from PotentialWell import PotentialWell
from PotWellSolver import PotWellSolver
from Plotter import Plotter
from pylab import *

GaAs = Compound(GaAsValues)
myPotWell = PotentialWell("z", 130, 1)
myPotWell2 = PotentialWell("x", 130, 1)
myPotWellSolver = PotWellSolver(GaAs, myPotWell, 4)
myPotWellSolver2 = PotWellSolver(GaAs, myPotWell2, 4)
myPotWellSolver.setGridPoints(100)
myPotWellSolver2.setGridPoints(100)
myPlotter = Plotter()


kVec = arange(-3, 3, 0.1)
def rotate(fractions, rotateTo="x"):
    a = np.sqrt(fractions[0])
    b = np.sqrt(fractions[1])
    c = np.sqrt(fractions[2])
    d = np.sqrt(fractions[3])

    aRot = None
    bRot = None
    cRot = None
    dRot = None
    if rotateTo == "x":
        aRot = np.sqrt(2)/4*a - np.sqrt(6)/4*b + np.sqrt(6)/4*c - np.sqrt(2)/4*d
        bRot = np.sqrt(6)/4*a - np.sqrt(2)/4*b - np.sqrt(2)/4*c + np.sqrt(6)/4*d
        cRot = np.sqrt(6)/4*a - np.sqrt(2)/4*b - np.sqrt(2)/4*c - np.sqrt(6)/4*d
        dRot = np.sqrt(2)/4*a + np.sqrt(6)/4*b + np.sqrt(6)/4*c + np.sqrt(2)/4*d
    elif rotateTo == "z":
        aRot = np.sqrt(2)/4*a + np.sqrt(6)/4*b + np.sqrt(6)/4*c + np.sqrt(2)/4*d
        bRot = -np.sqrt(6)/4*a - np.sqrt(2)/4*b - np.sqrt(2)/4*c + np.sqrt(6)/4*d
        cRot = np.sqrt(6)/4*a - np.sqrt(2)/4*b - np.sqrt(2)/4*c + np.sqrt(6)/4*d
        dRot = -np.sqrt(2)/4*a + np.sqrt(6)/4*b - np.sqrt(6)/4*c + np.sqrt(2)/4*d

    TotalDensity = aRot**2 + bRot**2 + cRot**2 + dRot**2
    HeavyHole = (aRot**2+dRot**2)/TotalDensity
    LightHole = (bRot**2 + cRot**2)/TotalDensity

    both = []
    both.append(HeavyHole)
    both.append(LightHole)
    return both


xValues = zeros((2, len(kVec)))
rotXValues = zeros((2, len(kVec)))
zValues = zeros((2, len(kVec)))
rotZValues = zeros((2, len(kVec)))

column = 0
for k in kVec:
    valsX = myPotWellSolver2.getMixing(k, 0)
    rotValsX = rotate(valsX, "z")

    valsZ = myPotWellSolver.getMixing(k, 0)
    rotValsZ = rotate(valsZ, "x")

    xValues[0][column] = valsX[0]+valsX[3]
    xValues[1][column] = valsX[2]+valsX[1]
    rotXValues[0][column] = rotValsX[0]
    rotXValues[1][column] = rotValsX[1]

    zValues[0][column] = valsZ[0]+valsZ[3]
    zValues[1][column] = valsZ[2]+valsZ[1]
    rotZValues[0][column] = rotValsZ[0]
    rotZValues[1][column] = rotValsZ[1]
    column += 1

plot(kVec, rotZValues.T)
plot(kVec, xValues.T)
show()