import numpy as np
from ase import *
from ase.io import read
from ase.calculators.eam import EAM
from gpaw.eigensolvers import Davidson
from gpaw import GPAW, PW

#at09 = read('POSCAR_0.9')
#at10 = io.read('POSCAR_1.0')
ata = read('POSCAR_1.1')


for k in range (2, 3):
	calc = GPAW(kpts = (k,k,k),
		mode = PW(50),
		h = 0.2,
		nbands = 162,
		xc = 'PBE',
		txt = 'HA5out1.txt')

	ata.set_calculator(calc)

	forces = ata.get_forces()
	with open('forcefileK','ab') as f: f.write('%f \ %f',(forces, k))

