import numpy as np
from ase import *
from ase.io import read
from eam_calculator import get_calc

#Orthorhombic shear
def ortShear(atoms, e = None):
	if e is None:
		cell = np.array([[ 12.15,   0.  ,   0.  ],
			[  0.  ,  12.15,   0.  ],
			[  0.  ,   0.  ,  12.15]])
	else:
		
		cell = atoms.get_cell()
		shear = np.array([[ e,   0.  ,   0.  ],
			[  0.  ,  -e,   0.  ],
			[  0.  ,   0.  , e**2/(1-e**2)]])
		cell = np.dot(cell,np.add(shear,np.identity(3))) 
		
	atoms.set_cell(cell, True)

	return(atoms)
#Monoclinic
def monoCShear(atoms, e = None):
	if e is None:
		cell = np.array([[ 12.15,   0.  ,   0.  ],
		 	[  0.  ,  12.15,   0.  ],
			[  0.  ,   0.  ,  12.15]])
	else:
		
		cell = atoms.get_cell()
		shear = np.array([[ 0.,   e/2.0  ,   0.  ],
			[  e/2.0  ,  0.,   0.  ],
			[  0.  ,   0.  , e**2/(4-e**2)]])
		cell = np.dot(cell,np.add(shear,np.identity(3)))
	atoms.set_cell(cell, True)
	
	return(atoms)
#Tetragonal 
def tetShear(atoms, e = None):
	if e is None:
		cell = np.array([[ 12.15,   0.  ,   0.  ],
			[  0.  ,  12.15,   0.  ],
			[  0.  ,   0.  ,  12.15]])
	else:
		
		cell = atoms.get_cell()
		shear = np.array([[ e,   0.  ,   0.  ],
			[  0.  ,  e,   0.  ],
			[  0.  ,   0.  , 0.]])
		cell = np.dot(cell,np.add(shear,np.identity(3))) 
		
	atoms.set_cell(cell, True)

	return(atoms)


def main():

	A = 1000
	lmbda = 2
	D = 5
	mu = 1
	p = [A, lmbda, D, mu]
	calc = get_calc(p)
	alAtoms = read('POSCAR_1.0')
	alAtoms.set_calculator(calc)
	conversion = 160.2177 #1eV = 160 GPa  (gigapascals) 

	
	strain = list(np.arange(0.0, 0.45, 0.01))
	E_0 = alAtoms.get_potential_energy()
	E_orth = np.array([])	
	E_mono = np.array([])
	E_tet = np.array([])	

	for i in range(0, len(strain)):
		
		alAtoms = ortShear(alAtoms, strain[i])
		E = alAtoms.get_potential_energy()*conversion
		vol = alAtoms.get_volume()    
		E_orth = np.append(E_orth,(E_0 - E)/vol)
		alAtoms = ortShear(alAtoms)	
		
		alAtoms = monoCShear(alAtoms, strain[i])
		E = alAtoms.get_potential_energy()*conversion
		vol = alAtoms.get_volume()     
		E_mono = np.append(E_mono,(E_0 - E)*2/(vol))
		alAtoms = ortShear(alAtoms)	
		
		alAtoms = tetShear(alAtoms, strain[i])
		E = alAtoms.get_potential_energy()*conversion
		vol = alAtoms.get_volume()     
		E_tet = np.append(E_tet,(E_0 - E)*2/(vol))
		alAtoms = ortShear(alAtoms)	
		
	
		se = str(E_orth[-1])+'\t'+str(E_mono[-1])+'\t'+str(E_tet[-1])+'\n'
		with open('shearEnergies', 'a') as myfile:
			myfile.write(se)
	# Fit the energy diff to the strain
	ort = np.polyfit(strain , E_orth, 2)	
	monoC = np.polyfit(strain , E_mono, 2)
	tet = np.polyfit(strain , E_tet, 2)	

	result = str(ort[0])+'\t'+str(ort[1])+'\t'+str(ort[2])+'\n'+str(monoC[0])+'\t'+str(monoC[1])+'\t'+str(monoC[2])+'\n'+str(tet[0])+'\t'+str(tet[1])+'\t'+str(tet[2])+'\n \n'
	with open('cRes', 'a') as myfile:
		myfile.write(result)
	
# Calculate C11, C12, C44 and
	C11 = (ort[0] + tet[0])/2.
	C12 = tet[0] - C11 
	C44 =  monoC[0]
	
	B = (C11 + 2*C12)/3.0
	C_prime = (C11 - C12)/2.0


	CC = str(C11)+'\t'+str(C12)+'\t'+str(C44)+'\t'+str(B)+'\t'+str(C_prime)+'\n'
	with open('cConst', 'a') as myfile:
		myfile.write(CC)
	

main()
