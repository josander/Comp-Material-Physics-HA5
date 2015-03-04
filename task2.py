import numpy as np
from ase import *
from ase.io import read
from ase.calculators.eam import EAM
from gpaw.eigensolvers import Davidson
from gpaw import GPAW, PW
from eam_calculator import get_calc

def getEMForces(p)
	at = read('POSCAR_1.0')
	calc = get_calc(p)
	at.set_calculator(calc)
	forces = at.get_forces()
		
	return(forces)
	
def main():
	
	atDFT = read('res_POSCAR_1.0')
	forcesDFT = atDFT.get_forces;

	pInitial = [1 2 3 4]
	
		









main()
