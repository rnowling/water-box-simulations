from sys import stdout, exit, stderr

from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *

import simtk.unit

def read_sim_fl(flname):
    globals = vars(simtk.unit)
    globals["NoCutoff"] = NoCutoff
    globals["CutoffNonPeriodic"] = CutoffNonPeriodic
    globals["CutoffPeriodic"] = CutoffPeriodic
    
    fl = open(flname)
    string = "\n".join(fl)
    fl.close()

    conf = eval(string, globals)

    return conf

def lowercase_keys(map):
    new_map = dict()
    for key, value in map.iteritems():
        if isinstance(value, dict):
            value = lowercase_keys(value)
        new_map[key.lower()] = value

    return new_map

def convert_steps(conf):
    sim_steps = int(conf["simulation"]["length"] / conf["simulation"]["timestep"])
    output_steps = int(conf["output"]["period"] / conf["simulation"]["timestep"])
    conf["simulation"]["steps"] = sim_steps
    conf["output"]["steps"] = output_steps

    return conf
    
# TODO
# add conf validation
# add factories

def simulate(conf):
    sim_conf = conf["simulation"]
    int_conf = conf["simulation"]["integrator"]
    bc_conf = sim_conf["boundaryconditions"]
    output_conf = conf["output"]

    pdb = PDBFile(conf["input"])

    forcefield = ForceField(*sim_conf["forcefields"])
    system = forcefield.createSystem(pdb.topology, nonbondedMethod=bc_conf["type"],
                                     nonbondedCutoff=bc_conf["cutoff"],
                                     rigidWater=sim_conf["rigidwater"])

    if int_conf["type"] == "langevin":
        integrator = LangevinIntegrator(int_conf["temperature"], int_conf["gamma"], sim_conf["timestep"])
    else:
        print "Unknown integrator", int_conf["type"]
        sys.exit(1)
        
    simulation = Simulation(pdb.topology, system, integrator)
    simulation.context.setPositions(pdb.positions)
    simulation.context.setVelocitiesToTemperature(sim_conf["initialtemperature"])

    output_period = output_conf["steps"]
    simulation.reporters.append(StateDataReporter(stdout, output_period, step=True, time=True,
                                                kineticEnergy=True, potentialEnergy=True,
                                                temperature=True))
    for output_type, flname in output_conf.iteritems():
        print output_type, flname
        if output_type == "period":
            continue
        elif output_type == "steps":
            continue
        elif output_type == "structure":
            continue
        elif output_type == "trajectory":
            reporter = DCDReporter(flname, output_period)
        elif output_type == "statistics":
            reporter = StateDataReporter(open(flname, "w"), output_period, step=True, time=True, kineticEnergy=True,
                                         speed=True, potentialEnergy=True, temperature=True)
        else:
            print "Print unknown output", output_type

        simulation.reporters.append(reporter)

    print len(simulation.reporters), "reporters"


    simulation.step(sim_conf["steps"])

    positions = simulation.context.getState(getPositions=True).getPositions()
    PDBFile.writeFile(pdb.topology, pdb.positions, open(output_conf["structure"], "w"))
    
if __name__ == "__main__":
    sim_fl = sys.argv[1]

    conf = read_sim_fl(sim_fl)
    conf = lowercase_keys(conf)
    conf = convert_steps(conf)

    print conf

    simulate(conf)

    
