from openmm.app.modeller import Modeller
from openmm.app.pdbfile import PDBFile
from openmm.app.forcefield import ForceField

pdb = PDBFile("test.pdb")
mod = Modeller(pdb.topology, pdb.positions)
ff = ForceField("force-field.offxml")
mod.addExtraParticles(ff)
