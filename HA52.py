import numpy as np
from ase import *
from ase.io import read
from ase.calculators.eam import EAM
from gpaw.eigensolvers import Davidson
from gpaw import GPAW, PW
from scipy.optimize import leastsq


DFT = read('res_POSCAR_1.1.traj')
EAM = read('POSCAR_1.1')
p = [1000, 2, 5, 1] # eV, Å^-1, Å, Å^-1

calc = eam_calc(p) 
EAM.set_calculator(calc)

DFTforce = DFT.get_forces()
EAMforce = EAM.get_forces()

DFTce = DFT.get_
EAMce =	EAM.get_

