import numpy as np
from ase import *
from ase.io import read
from ase.calculators.eam import EAM
from gpaw.eigensolvers import Davidson
from gpaw import GPAW, PW

# Find convergence for k, Ecut and h
al = read('POSCAR_1.1')

outputFile = 'convH1.1out.txt'

# Find convergence for h
k = 4
Ecut = 100

gridSpaces = [0.2]##, 0.3, 0.4, 0.5, 0.6]

for gridSpace in gridSpaces:
	calc = GPAW(kpts = (k,k,k),
		mode = PW(Ecut),
		h = gridSpace,
		nbands = 320,
		xc = 'PBE',
		txt = outputFile)

	al.set_calculator(calc)

	forces = al.get_forces()

	meanF = 0

	for atom in range(1, 108):
		Forces[atom] = sqrt(forces[atom,1]*forces[atom,1]+forces[atom,2]*forces[atom,2]+forces[atom,3]*forces[atom,3])
		meanF = meanF + Forces[atom]

	print("%f \t %d" %(meanF/108, gridSpace))

