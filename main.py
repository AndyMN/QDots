from baseCompounds import *
from PotentialWell import PotentialWell
from PotWellSolver import PotWellSolver
from Plotter import Plotter
from pylab import *
from scipy.sparse import diags



#GaAs = Compound(GaAsValues)
#myPotWell = PotentialWell("x", 130, 1)
#myPotWellSolver = PotWellSolver(GaAs, myPotWell, 6)
#myPotWellSolver.setGridPoints(100)
#myPlotter = Plotter()

#kVec = arange(-3, 3, 0.1)

#myPlotter.plotMixing(myPotWellSolver, kVec)
#myPlotter.displayPlots()

minisize = 350
pi = np.pi
y1 = 3.12
y2 = 5.13
y3 = 9.5
delta = 44
kz = 1
stepsize = 0.01
ky = 0

diagR = (np.sqrt(3)*y2/pi**2) * (-2/stepsize**2)
subdiagR = (np.sqrt(3)/pi**2)*(y2/stepsize**2 - y3*ky/stepsize)
superdiagR = (np.sqrt(3)/pi**2)*(y2/stepsize**2 + y3*ky/stepsize)
R = zeros((minisize, minisize), dtype=complex)
i, j = indices(R.shape)
R[i == j] = diagR
R[i == j-1] = superdiagR
R[i == j+1] = subdiagR

diagP = (y1/pi**2)*(kz**2) + 2./(stepsize**2)
subdiagP = -y1/(stepsize**2*pi**2)
superdiagP = subdiagP
P = zeros((minisize, minisize), dtype=complex)
i, j = indices(P.shape)
P[i == j] = diagP
P[i == j-1] = superdiagP
P[i == j+1] = subdiagP

diagQ = (y2/pi**2) * ((-2.*kz**2)+(2/stepsize**2))
subdiagQ = -y2/(stepsize**2*pi**2)
superdiagQ = subdiagQ
Q = zeros((minisize, minisize), dtype=complex)
i, j = indices(Q.shape)
Q[i == j] = diagQ
Q[i == j-1] = superdiagQ
Q[i == j+1] = subdiagQ

diagS = -2.*y3*np.sqrt(3)*1j*ky*kz/pi**2
subdiagS = np.sqrt(3)*y3*kz*1j/(pi**2*stepsize)
superdiagS = -np.sqrt(3)*y3*kz*1j/(pi**2*stepsize)
S = zeros((minisize, minisize), dtype=complex)
i, j = indices(S.shape)
S[i == j] = diagS
S[i == j-1] = superdiagS
S[i == j+1] = subdiagS

A = np.bmat([[P+Q, -S, R, zeros((minisize,minisize)), -S/np.sqrt(2), np.sqrt(2)*R],
             [-S.conj().T, P-Q, zeros((minisize,minisize)), R, -np.sqrt(2)*Q, np.sqrt(3./2.)*S],
             [R.conj().T, zeros((minisize,minisize)), P-Q, S, np.sqrt(3./2.)*S.conj().T, np.sqrt(2)*Q],
             [zeros(minisize,minisize), R.conj().T, S.conj().T, P+Q, -np.sqrt(2)*R.conj().T, -S.conj().T/np.sqrt(2)],
             [-S.conj().T/np.sqrt(2), -np.sqrt(2)*Q.conj().T, np.sqrt(3./2.)*S, -np.sqrt(2)*R, P, zeros((minisize,minisize))],
             [np.sqrt(2)*R.conj().T, np.sqrt(3./2.)*S.conj().T, np.sqrt(2)*Q.conj().T, -S/np.sqrt(2), zeors((minisize,minisize)), P]])

w, v = eigh(A)










