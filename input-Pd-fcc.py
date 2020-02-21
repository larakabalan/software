import numpy as np
from ase import Atoms
from ase.io.trajectory import Trajectory
from ase.calculators.aims import Aims
from ase.visualize import view
from ase.build import fcc111, bulk, surface
from ase.constraints import FixAtoms, FixBondLengths
from ase.io import write
from ase.io import read 

a = 3.909 # approximate lattice constant
b = a / 2   

calc = Aims(xc='pbe',
           spin='none',
           k_grid=(9,9,9),
           vdw_correction_hirshfeld="True",
           relativistic=('atomic_zora','scalar'),
           #use_dipole_correction='True',
           compute_forces="true",
           output=['mulliken'],
          # elsi_restart=("write",1)
           )

bulk = Atoms('Pd',
           cell=[(0, b, b), (b, 0, b), (b, b, 0)],
           pbc=1,
           calculator=calc)  # use EMT potential

bulk.get_potential_energy()
#traj.write(bulk)
energy_bulk=bulk.get_potential_energy() 
print ("potential-energy-bulk", energy_bulk)
#print (energy_bulk)

############SURFACE################################
from ase.io import write

## Edit these
atomic_species='Pd'
a2 = a                          
unit_cell_depth=3
unit_cell_width=3
slab_depth=4
vacuum_region_size=10.0

## Create surface
slab = fcc111(atomic_species, a=a2, size=(unit_cell_width,unit_cell_depth,slab_depth))
#slab = fcc100(atomic_species, a=lattice_parameter*sqrt(2), size=(unit_cell_width,unit_cell_depth,slab_depth))

## Add vacuum
slab.center(vacuum=vacuum_region_size, axis=2)
mask0 = [atom.tag > 2 for atom in slab]
constraint0 = FixAtoms(mask=mask0)
slab.set_constraint([constraint0])

calc2 = Aims(xc='pbe',
           spin='none',
           k_grid=(9,9,9),
           vdw_correction_hirshfeld="True",
           relativistic=('atomic_zora','scalar'),
          #use_dipole_correction='True',
           compute_forces="true",
           output=['mulliken'],
          # elsi_restart=("write",1)
           )
slab.set_calculator(calc2)
slab.get_potential_energy()

energy_slab=slab.get_potential_energy()
print(slab.get_potential_energy())
#print (energy_slab)

######claculating the energy surface through  the 2 models: the first one is through energy per atom and the second one is energy by surface area which is a*a*sinO with O is 90 for 100 and 110 faces while it is 60 for 111 face

print ("energy of surface Esur1: ", (0.5*(energy_slab-36*energy_bulk))/9,  "ev/atom")
print ("energy of surface Esur2: ", (0.5*((energy_slab-36*energy_bulk)/(a*a*np.sin(np.deg2rad(60)))),  "ev/((A^2)"))
#traj.write(slab.traj)

