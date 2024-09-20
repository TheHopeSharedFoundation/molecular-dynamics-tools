############################################################################
# Copyright 2024 The Hope Shared Foundation                                #
#                                                                          #
# Licensed under the Apache License, Version 2.0 (the "License");          #
# you may not use this file except in compliance with the License.         #
# You may obtain a copy of the License at                                  #
#                                                                          #
# http://www.apache.org/licenses/LICENSE-2.0                               #
#                                                                          #
# Unless required by applicable law or agreed to in writing, software      #
# distributed under the License is distributed on an "AS IS" BASIS,        #
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. #
# See the License for the specific language governing permissions and      #
# limitations under the License.                                           #
#                                                                          #
# The Hope Shared Foundation, September 2024                               #
# https://hopesharedfoundation.org/                                        #
# Towards Alleviating Suffering                                            #
############################################################################

########### Import OpenMM libraries
from openmm.app import *
from openmm import *
from openmm.unit import *
########### Import molecule parameterization libraries
from openff.toolkit.topology import Molecule, Topology
from openmmforcefields.generators import SMIRNOFFTemplateGenerator # For examples generating small molecule parameters using SMIRNOFF: https://github.com/openmm/openmmforcefields
########### Import other useful libraries
from sys import stdout
from datetime import datetime
import plotly.express as px
import pandas as pd
import argparse
import os

parser = argparse.ArgumentParser(prog='gpu-accelerated molecular dynamics through OpenMM',description='This program runs the following: minimization, NVT ensemble equilibration (1ns),  NPT ensemble equilibration (1ns), and NPT ensemble production',epilog='Copyright 2024 The Hope Shared Foundation')

now= datetime.now()
current_time = now.strftime("%Y%b%dUTC%H%M%S").capitalize()


parser.add_argument('--temperature',default=300,help="The simulation temperature (thermostat and barostat) in units of kelvin.  The default is 300K.",required=False)
parser.add_argument('--pressure',default=1,help="The simulation pressure in units of bar.  The default is 1 bar.",required=False)
parser.add_argument('--time',default=120,help="The simulation time in units of hour.  The default is 120 hours.",required=False)
parser.add_argument('--basename',default=now,help="The file basename for all simulation output files.  The default is a date-timestamp",required=True)
parser.add_argument('--receptor',help="A protein data bank (PDB) formatted file of the receptor.  This is a required argument.",required=True)
parser.add_argument('--ligand',default="NULL",help="A structure-data formatted file (SDF) of the ligand.",required=False)
parser.add_argument('--scratch',default="/home/user/scratch/",help="This is the path to the very fast scatch directory on the server.  These files are not automatically backed up.  Make sure to move these files BEFORE deprovisioning the server.",required=False)
args = parser.parse_args()

############################### PROJECT SPECIFIC PARAMETERS.  MODIFY THE FOLLOWING AS NECESSARY.

base_file_name = str(args.basename)
receptor_file = str(args.receptor)
pdb = PDBFile(receptor_file) # mixture.pdb is produced by PackMol and surrounds a protein with 0.2 M benzene.

if args.ligand != 'NULL':
    ligand_file = str(args.ligand)
    molecule = Molecule.from_file(ligand_file) # https://docs.openforcefield.org/projects/toolkit/en/stable/users/molecule_cookbook.html#from-sdf-file-object
    smirnoff = SMIRNOFFTemplateGenerator(molecules=molecule)

#####################################################

now= datetime.now()
current_time = now.strftime("%Y%b%dUTC%H%M%S").capitalize()
conda_explicit_command = 'conda list --explicit > ' + args.scratch + base_file_name + '_CondaSpecFile.txt'
os.system(conda_explicit_command)
log_file_all = args.scratch + base_file_name + "_log_" + current_time + ".txt"
log_file_NVT_equilibration = args.scratch + base_file_name + "_NVT_equilibration_" + current_time + ".csv"  
log_file_NPT_equilibration = args.scratch + base_file_name + "_NPT_equilibration_" + current_time + ".csv"
log_file_NPT_production = args.scratch + base_file_name + "_NPT_production_" + current_time + ".csv"
#log_file_entire_simulation
chk_file_NPT_production = args.scratch + base_file_name + "_NPT_production_" + current_time + ".chk" 
xml_file_NPT_production = args.scratch + base_file_name + "_NPT_production_" + current_time + ".xml"
dcd_file_NPT_production = args.scratch + base_file_name + "_NPT_production_" + current_time + ".dcd"
xtc_file_NPT_production = args.scratch + base_file_name + "_NPT_production_" + current_time + ".xtc"
pdb_file_initial = args.scratch + base_file_name + "_initialStructure_" + current_time + ".pdb"
pdb_file_minimized = args.scratch + base_file_name + "_minimizedStructure_" + current_time + ".pdb"
pdb_file_final = args.scratch + base_file_name + "_finalStructure_" + current_time + ".pdb"
write_interval_dcd = 12500 # Units of steps.  Write to the dcd file every 0.004 picoseconds * 12500 steps = 50 picoseconds.
write_interval_xtc = 12500 # Units of steps.  Write to the dcd file every 0.004 picoseconds * 12500 steps = 50 picoseconds.
write_interval_chk = 1 # Units of hours.  Write to the checkpoint file every 1 hour.
write_interval_NVT_equilibration = 250 # Units of steps.  Write to the NVT output file every 0.004 picoseconds * 250 steps = 1 picosecond.
write_interval_NPT_equilibration = 250 # Units of steps.  Write to the NPT output file every 0.004 picoseconds * 250 steps = 1 picosecond.
write_interval_NPT_production = 250 # Units of steps.  Write to the NPT output file every 0.004 picoseconds * 250 steps = 1 picosecond.
stdout_interval_NVT_equilibration = 250 # Units of steps.  Write to the NVT output file every 0.004 picoseconds * 250 steps = 1 picosecond.
stdout_interval_NPT_equilibration = 250 # Units of steps.  Write to the NPT output file every 0.004 picoseconds * 250 steps = 1 picosecond.
stdout_interval_equilibration = 250# Units of steps.  Write to the equilibration output file every 0.004 picoseconds * 250 steps = 1 picosecond.
stdout_interval_NPT_production = 250 # Units of steps.  Write to the NPT output file every 0.004 picoseconds * 250 steps = 1 picosecond.
total_steps_NVT_equilibration = 250000 # Units of steps.  Run NVT equilibration stage for 0.004 picoseconds * 250000 steps = 1 nanosecond.
total_steps_NPT_equilibration = 250000 # Units of steps.  Run NPT equilibration stage for 0.004 picoseconds * 250000 steps = 1 nanosecond.
total_steps_equilibration = total_steps_NVT_equilibration + total_steps_NPT_equilibration 
total_run_time_NPT_production = args.time # Units of hours.  
simulation_temperature = args.temperature # Units of kelvin.
simulation_pressure = args.pressure # Units of bar.
simulation_stepsize = 0.004 # Units of picosecond

#####################################################

now= datetime.now()
current_time = now.strftime("%Y%b%dUTC%H%M%S").capitalize()
print(f"Setting up the simulation...{current_time}")
with open(log_file_all, 'a') as log_output:
    print(f"Setting up the simulation...{current_time}",file=log_output)


################## Specify forcefield

forcefield = ForceField('amber14-all.xml', 'amber14/tip3pfb.xml')
if args.ligand != 'NULL':
    forcefield.registerTemplateGenerator(smirnoff.generator)

################## Create model (topology and positions; protein, ligand, and solvent)

modeller = Modeller(pdb.topology, pdb.positions) # Place topology and positions of protein and benzenes into "modeller"
modeller.addSolvent(forcefield,'tip3p',padding=1.0*nanometer,boxShape='dodecahedron') # http://docs.openmm.org/development/api-python/generated/openmm.app.modeller.Modeller.html

################### Create system

system = forcefield.createSystem(modeller.topology, nonbondedMethod=PME,removeCMMotion=True,rigidWater=True,nonbondedCutoff=0.9*nanometer, constraints=HBonds,hydrogenMass=1.5*amu,ewaldErrorTolerance=0.0005)

################### Define integrator

integrator = LangevinMiddleIntegrator(simulation_temperature*kelvin, 1/picosecond, simulation_stepsize*picoseconds)

################### Create simulation (model, system, and itegrator)

simulation = Simulation(modeller.topology, system, integrator)
simulation.context.setPositions(modeller.positions)

################### Write initial model file

now= datetime.now()
current_time = now.strftime("%Y%b%dUTC%H%M%S").capitalize()
with open(log_file_all, 'w') as log_output_1:
    print(f"file basename: {args.basename}",file=log_output_1)
    print(f"receptor file: {args.receptor}",file=log_output_1)
    print(f"ligand file: {args.ligand}",file=log_output_1)
    print(f"Start of script..{current_time}",file=log_output_1)
    print(f"thermostat temperature: {simulation_temperature} K",file=log_output_1)
    print(f"barostat temperature: {simulation_temperature} K",file=log_output_1)
    print(f"barostat pressure: {simulation_pressure} bar",file=log_output_1)
    print(f"NPT ensemble equilibration steps: {total_steps_NVT_equilibration} steps",file=log_output_1)
    print(f"NPT ensemble equilibration steps: {total_steps_NPT_equilibration} steps",file=log_output_1)
    print(f"NPT ensemble production runtime: {total_run_time_NPT_production} hours",file=log_output_1)

initial_positions = simulation.context.getState(getPositions=True,enforcePeriodicBox=True).getPositions()
PDBFile.writeFile(modeller.topology,initial_positions,file=pdb_file_initial,keepIds=True)

now= datetime.now()
current_time = now.strftime("%Y%b%dUTC%H%M%S").capitalize()
print(f"Setting up of the simulation is complete.  Starting energy minimization...{current_time}")
with open(log_file_all, 'a') as log_output_2:
    print(f"Setting up of the simulation is complete.  Starting energy minimization...{current_time}",file=log_output_2)

simulation.minimizeEnergy()

################## Write minimized model file

minimized_positions = simulation.context.getState(getPositions=True,enforcePeriodicBox=True).getPositions()
PDBFile.writeFile(modeller.topology,minimized_positions,file=pdb_file_minimized,keepIds=True)

now= datetime.now()
current_time = now.strftime("%Y%b%dUTC%H%M%S").capitalize()
print(f"Energy minimization is complete.  Starting NVT equilibration...{current_time}")
with open(log_file_all, 'a') as log_output_3:
    print(f"Energy minimization is complete.  Starting NVT equilibration...{current_time}",file=log_output_3)

simulation.reporters.append(StateDataReporter(stdout, stdout_interval_equilibration, step=True, time=True, potentialEnergy=True, kineticEnergy=True, totalEnergy=True, temperature=True, progress=True, volume=True,density=True, remainingTime=True, speed=True, totalSteps=total_steps_NPT_equilibration, separator="\t")) # Progress = 100% means NVT equilibration is complete, progress = 200% means NPT equilibration is complete, and anything above 200% represents the % above total_setps_NPT_equilibration

simulation.reporters.append(StateDataReporter(log_file_NVT_equilibration, write_interval_NVT_equilibration, step=True, time=True, potentialEnergy=True, kineticEnergy=True, totalEnergy=True, temperature=True, progress=True, volume=True,density=True, remainingTime=True, speed=True, totalSteps=total_steps_NVT_equilibration, separator=","))


simulation.step(total_steps_NVT_equilibration)

now= datetime.now()
current_time = now.strftime("%Y%b%dUTC%H%M%S").capitalize()
print(f"NVT equilibration is complete.  Starting NPT...{current_time}")
with open(log_file_all, 'a') as log_output_4:
    print(f"NVT equilibration is complete.  Starting NPT...{current_time}",file=log_output_4)


system.addForce(MonteCarloBarostat(simulation_pressure*bar,simulation_temperature*kelvin))
simulation.context.reinitialize(preserveState=True) # http://docs.openmm.org/7.2.0/api-python/generated/simtk.openmm.openmm.Context.html#simtk.openmm.openmm.Context.reinitialize

simulation.step(total_steps_NPT_equilibration)

simulation.reporters.append(StateDataReporter(log_file_NPT_equilibration, write_interval_NPT_equilibration, step=True, time=True, potentialEnergy=True, kineticEnergy=True, totalEnergy=True, temperature=True, progress=True, volume=True,density=True, remainingTime=True, speed=True, totalSteps=total_steps_NPT_equilibration, separator=","))

now= datetime.now()
current_time = now.strftime("%Y%b%dUTC%H%M%S").capitalize()
print(f"NPT equilibration is complete. Starting NPT production...{current_time}")
with open(log_file_all, 'a') as log_output_5:
    print(f"NPT equilibration is complete. Starting NPT production...{current_time}",file=log_output_5)

simulation.reporters.append(StateDataReporter(log_file_NPT_production, write_interval_NPT_production, step=True, time=True, potentialEnergy=True, kineticEnergy=True, totalEnergy=True, temperature=True, progress=True, volume=True,density=True, remainingTime=True, speed=True, totalSteps=total_steps_equilibration, separator=","))

simulation.reporters.append(DCDReporter(dcd_file_NPT_production,write_interval_dcd,enforcePeriodicBox=True))
simulation.reporters.append(XTCReporter(xtc_file_NPT_production,write_interval_xtc,enforcePeriodicBox=True))

final_positions = simulation.context.getState(getPositions=True,enforcePeriodicBox=True).getPositions()
PDBFile.writeFile(modeller.topology,final_positions,file=pdb_file_final,keepIds=True)

simulation.runForClockTime(total_run_time_NPT_production*hours,checkpointFile=chk_file_NPT_production,stateFile=xml_file_NPT_production,checkpointInterval=write_interval_chk)

#simulation.reporters.append(StateDataReporter(stdout, stdout_interval_NPT_production, step=True, time=True, potentialEnergy=True, kineticEnergy=True, totalEnergy=True, temperature=True, volume=True,density=True, speed=True, separator="\t"))

now= datetime.now()
current_time = now.strftime("%Y%b%dUTC%H%M%S").capitalize()
print(f"NPT production and entire run are complete...{current_time}")
with open(log_file_all, 'a') as log_output_6:
    print(f"NPT production and entire run are complete...{current_time}",file=log_output_6)
