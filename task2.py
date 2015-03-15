import numpy as np
from ase import *
from ase.io import read
from ase.calculators.eam import EAM
from gpaw.eigensolvers import Davidson
from gpaw import GPAW, PW
from eam_calculator import get_calc
import scipy.optimize as opt
import sys

# Calculates the residual of the foreces of (apos) with (p) params compared to DFTF 
def residual(p, aPos,DFTF, w1, w2):
	
	pI = (p[0], p[1], p[2], p[3])
	calc = get_calc(pI)
	aPos.set_calculator(calc)
	EAMF = aPos.get_forces()
	#Get the latParam and cohesive energy 
	latParam = aPos.get_cell()[0][0]/3
	cohE = aPos.get_potential_energy()/aPos.get_number_of_atoms()

	EAMF = EAMF.ravel()
	EAMF = EAMF.reshape((len(EAMF), 1))
	EAMF = EAMF.flatten()
	np.append(EAMF,[latParam,cohE])	
	
	res = list(np.subtract(DFTF, EAMF))	
	res[-1] = res[-1]*w2
	res[-2] = res[-2]*w1
	return(res)

	
def main(argv):

	atDFT = read('res_POSCAR_1.0.traj')
	at = read('POSCAR_1.0')
	forcesDFT = atDFT.get_forces()
	forcesDFT = forcesDFT.ravel()
	forcesDFT = forcesDFT.reshape((len(forcesDFT), 1))
	forcesDFT = forcesDFT.flatten()	

	latParam =  4.032 
	cohE = 3.36	
	np.append(forcesDFT, [latParam,cohE ])	

	#Initial guesses
	A = 1000
	lmbda = 2
	D = 5
	mu = 1

	
	if len(argv) == 2:
		w1 = float(argv[0])**0.5
		w2 = float(argv[1])**0.5
	else:
		w1 = 1**0.5 #weight for lattice param
		w2 = 1**0.5 #weight for cohE


	pInitial = [A, lmbda, D, mu]
	
	# Residual minimasation
	params = opt.leastsq(residual, pInitial, args=(at, forcesDFT, w1,w2))

	result = str(params[0][0])+'\t'+str(params[0][1])+'\t'+str(params[0][2])+'\t'+str(params[0][3])+'\t'+str(w1**2)+'\t'+str(w2**2)+'\n'
	with open('params', 'a') as myfile:
		myfile.write(result)

if __name__ == "__main__":
   main(sys.argv[1:])
