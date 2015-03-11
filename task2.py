import numpy as np
from ase import *
from ase.io import read
#from ase.calculators.eam import EAM
#from gpaw.eigensolvers import Davidson
#from gpaw import GPAW, PW
from eam_calculator import get_calc
import scipy.optimize as opt

# Gets the foces of an atoms object (at) using the (p) params
def getEMForces(p, at):
	calc = get_calc(p)
	at.set_calculator(calc)
	forces = at.get_forces()
		
	return(forces)
# Calculates the residual of the foreces of (apos) with (p) params compared to DFTF 
def residual(p, aPos,DFTF):
	EAMF = getEMForces(p, aPos)
	#res = [a - b for a, b in zip(DFTF, EAMF)]
	res = np.subtract(DFTF, EAMF)
	return(res)

	
def main():
	# Pseudokod som antagligen inte fungerar. Har inte vagat testa =(
	atDFT = read('res_POSCAR_1.0.traj')
	at = read('POSCAR_1.0')
	forcesDFT = atDFT.get_forces()

	A = 10000
	lmbda = 3
	D = 5
	mu = 1
	pInitial = (A, lmbda, D, mu)
	params = opt.leastsq(residual, pInitial, args=(at, forcesDFT))
		

main()
