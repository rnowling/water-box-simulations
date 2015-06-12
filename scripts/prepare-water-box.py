import argparse
from sys import stdout, exit, stderr

from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *

def parse_args():
    parser = argparse.ArgumentParser(description="Generate and equilibrate a water box")
    parser.add_argument("-T", "--temperature", action="store", default=298.15, type=float,
                        required=False, help="Simulation temperature")
    parser.add_argument("-S", "--steps", action="store", default=100000, type=int,
                        required=False, help="Equilibration simulation length")
    parser.add_argument("-d", "--dim", action="store", nargs=3, type=float, required=True,
                        help="Specify widths of x, y, and z dimensions of box in nm")
    parser.add_argument("-o", "--output", action="store", type=str, required=True,
                        help="Resulting PDB file")
    parser.add_argument("-s", "--structure", action="store", type=str, required=True,
                        help="PDB structure file")

    args = parser.parse_args()
    
    return vars(args)

def create_system(structure_file, dim):
    pdb = PDBFile(structure_file)
    modeller = Modeller(pdb.topology, pdb.positions)
    water_forcefield = ForceField("amber96.xml", "tip3p.xml")
    modeller.addSolvent(water_forcefield, boxSize=dim)
    
    n_atoms = len(list(modeller.topology.atoms()))
    print n_atoms, "atoms"
    density = dim[0] * dim[1] * dim[2] * (10.0 * 10.0 * 10.0) / n_atoms
    print density, "A^3 / atom"

    return modeller.topology, modeller.positions

def equilibrate(topology, positions, temp, steps):
    forcefield = ForceField("amber96.xml", "tip3p.xml")
    system = forcefield.createSystem(topology, nonbondedMethod=CutoffPeriodic,
                                     rigidWater=False)
    integrator = LangevinIntegrator(temp * kelvin, 1.0 / picoseconds, 0.001 * picoseconds)
    simulation = Simulation(topology, system, integrator)
    simulation.context.setPositions(positions)
    simulation.context.setVelocitiesToTemperature(temp * kelvin)
    simulation.reporters.append(StateDataReporter(stdout, 1000, step=True,
        potentialEnergy=True, temperature=True))
    simulation.step(steps)

    positions = simulation.context.getState(getPositions=True).getPositions()

    return positions
    

if __name__ == "__main__":
    args = parse_args()

    print "Args", args

    topology, positions = create_system(args["structure"], args["dim"])

    eq_positions = equilibrate(topology, positions, args["temperature"], args["steps"])

    PDBFile.writeFile(topology, eq_positions, open(args["output"], "w"))
