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
Ecut = 90

gridSpaces = [0.2 0.3 0.4 0.5]

for gridSpace in gridSpaces:
	calc = GPAW(kpts = (k,k,k),
		mode = PW(Ecut),
		h = gridSpace,
		nbands = 320,
		xc = 'PBE',
		txt = outputfile)

	al.set_calculator(calc)

	forces = al.get_forces()
	print('%f \ %f',(forces, gridSpace));

