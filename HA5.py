import numpy as np
from ase import *
from ase.io import read
from ase.calculators.eam import EAM
from gpaw.eigensolvers import Davidson
from gpaw import GPAW, PW

#at09 = read('POSCAR_0.9')
#at10 = io.read('POSCAR_1.0')
ata = read('POSCAR_1.1')

f = open('forcefileK', 'w')

for k in range (2, 10):
	calc = GPAW(kpts = (k,k,k),
		mode = PW(500),
		h = 0.2,
		nbands = 162,
		xc = 'PBE',
		txt = 'HA5out.txt')

	ata.set_calculator(calc)

	forces = ata.get_forces()

	f.write('%f \ %f',(forces, k))

f.closed
