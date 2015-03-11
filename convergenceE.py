import numpy as np
from ase import *
from ase.io import read
from ase.calculators.eam import EAM
from gpaw.eigensolvers import Davidson
from gpaw import GPAW, PW

# Find convergence for Ecu
al = read('POSCAR_1.1')

outputFile = 'convE1.1out.txt'

# Change these befonre running the script
gridSpace = 0.3
k = 4

EnergiesCutOff = [50]## 70 90 110 130 150 170 190 210 230 250 270]

for Ecut in EnergiesCutOff:
	calc = GPAW(kpts = (k,k,k),
		mode = PW(Ecut),
		h = gridSpace,
		nbands = 320,
		xc = 'PBE',
		txt = outputFile)

	al.set_calculator(calc)
	forces = al.get_forces()

	print('%f \ %f',(forces, Ecut));


