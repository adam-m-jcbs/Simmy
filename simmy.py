# This is intended to be a generalized python infrastructure for managing a
# suite or grid of simulations.  The design to to have high-level abstract
# classes that take care of basics and that template the management, while
# implementations of the abstract classes will handle specific machines, codes,
# etc.  This is based in part on the infrastructure I developed for managing a
# suite of 3D Maestro models that ran on multiple supercomputers as well as
# debug runs on local workstations or laptops.  Those models are described in
# http://adsabs.harvard.edu/abs/2016ApJ...827...84J

# Author:        Adam Jacobs
# Creation date: September 20, 2017

#Notes on conventions, programming style:
#  1) All classes are in CamelCase with the first letter capitalized.  Class methods
#     are also in CamelCase with the first letter lower-case, e.g. myMethod().
#  2) Non-class functions and all variables are lowercase with underscores (_) 
#     acting as delimiters when needed/wanted.
#  3) An underscore prefix, e.g. _variable, means the named item is intended to be
#     private (Python doesn't provide for easily enforcing data hiding, it's all by 
#     convention).
#  4) Names that are in ALL_CAPS are intended as constants.

# Usage: load as a module
# Requirements
#    + Python 3
#    + Common scientific python tools: NumPy, matplotlib
#    TODO: fill in other requirements

#Global TODO
#   + As currently designed, simmy imposes a directory structure convention.
#   This can be a bit dangerous as I rely on this assumption quite a bit (e.g. I
#   will infer simulation labels from directory labels).  A user manually
#   changing a directory name can break pretty much everything.  Is it worth it
#   to redesign?  I personally like systematically enforcing a directory
#   structure/organization, but I predict this to be a point of fragility for
#   users that aren't me.

###########################################
### Global Imports, Data, and Constants ###
###########################################
import sys

# Global checks, assertions
if not sys.version_info >= (3,):
    #TODO if you can, perhaps make this work for 2.7+ and 3 with from __future__
    raise RuntimeError("simmy requires Python 3")

#ASCII escape sequences for color text
START_RED   = '\033[31m'
START_GREEN = '\033[32m'
START_BLUE  = '\033[34m'
RESET       = '\033[0m'

##################################
### Simulation Management Code ###
##################################
# The code in this section is all about managing simulations.  I've sectioned it
# off as I may one day break it off into its own module.  TODO If this gets
# hefty enough, separate into a module.

class SimulationGrid(object):
    """A suite or grid of related simulations and the different operations one
    might want to perform on such such a suite.
    
    This class serves as a base class for a grid of a specific code's
    simulations.  Specific codes should subclass this and add functionality
    particular to them while making use of the interface and generically useful
    methods/properties provided by SimulationGrid."""
    ##Global Data##
 
    ##Constructor##
    def __init__(self, label, stage_base, scratch_base=None):
        """Initialize the SimulationGrid object.
 
           Why a staging directory and a scratch directory?
           On many machines there is a large filesystem (scratch) for running
           simulations that is purged and a smaller filesystem (stage) that is
           backed up.  This is the motivation for having a staging and scratch
           directory.
           
           The setup for simulations and their reduced (small size) output will
           be kept in staging, while the actual run and large output will be
           executed in scratch.  Large output that can't fit in staging should be
           stored in high performance storage that is available on most machines
           that have the stage/scratch configuration.  Even for systems that
           don't have a purged scratch directory, this is a nice way to keep a
           clean space with simulation essentials and a place to tweak and
           experiment on a run.
 
           Arguments:
           self         --> implicitly passed reference to this instance of SimulationGrid
           label        --> this grid's label
           stage_base   --> path to base directory containing all simulations in the grid
           scratch_base --> path to base directory where the work is done but data's purged (optional)
                            If not provided, it is assumed simulations are run in
                            the staging directories.
           """
        self._stage_base = stage_base
        if scratch_base is None:
            self._scratch_base = self._stage_base
            self._use_scratch = False
        else:
            self._scratch_base = scratch_base
            self._use_scratch = True
        self._label = label
        self._my_sims = self._getActiveSims()
 
    ##Public methods##
    def listSimulations(self):
        """Print a list of simulations in this grid.
        
        This will only list 'active' simulations that the user is still
        exploring.  Simulations that are no longer being explored can be archived
        so that they will not pollute the list.
        
        This basic implementation just lists the labels.  It is advised to
        override this in a subclass to add other relevant information.
        """
        #from supercomputer import TermColors, Supercomputer #TODO delete
        from os.path import isfile, isdir, join
        from glob import glob
 
        active_sims = self._getActiveSimDirs()
        
        if self._use_scratch:
            heading = '{0:29s}|{1:14s}'.format('Label', 'In scratch?')
            list_format = '{0:29s}|{1:14s}'
        else:
            heading = '{0:29s}'.format('Label')
            list_format = '{0:29s}'
        yep    = START_GREEN + '{0:14s}'.format("Yes!")    + RESET
        nope   = START_RED   + '{0:14s}'.format("No!")     + RESET
        purged = START_BLUE  + '{0:14s}'.format("Purged!") + RESET
        
        print(heading)
        for r in active_runs:
            if self._use_scratch:
                #Check scratch TODO mention purge check for subclasses?
                rundir = join(self._scratch_base, r)
                if isdir(rundir):
                    if self._inScratch(rundir):
                        sc_str = yep
                    else:
                        sc_str = purged
                    #TODO add to Maestro or SubChandra subclass method
                    #found_all_expected_files = (
                    #      len(glob( join(rundir, 'main.*') ))         > 0 and
                    #      len(glob( join(rundir, 'inputs*') ))        > 0 and
                    #      len(glob( join(rundir, 'helm_table.dat') )) > 0
                    #      )
                    #if found_all_expected_files:
                    #    sc_str = yep
                    #else:
                    #    sc_str = purged
                else:
                    sc_str = nope
                outstr = list_format.format(r, sc_str) 
            else:
                outstr = list_format.format(r) 

            #Check queue TODO Add this to an e.g. Maestro simulation subclass
            #  ASSUMPTION: run directory is same as queue label
            #if curcomp.isQueued(r):
            #    q_str = yep
            #else:
            #    q_str = nope
 
            print(outstr)
 
    def printStatus(self):
        """Print the current status of all simulations in the grid"""
        #Loop over all runs in staging directory assuming <run label>/[output, run, plots]
        #structure
        sims = self._getActiveSims()
        for s in sims:
            print(s._getStatusSummary()
 
        return
 
    ##Private methods##
    #TODO The below is for SubChandar sim, move to a subclass then delete
    #def _getRStat(self, run_label):
    #   """Determine run status of run_label"""
    #   from supercomputer import Supercomputer
    #   
    #   curcomp = Supercomputer.getCurrentSC()
    #   (job_status, qcount) = curcomp.getQStatus(run_label)
 
    #   return job_status + '*{0:02d}'.format(qcount)
 
    #def _getTPeak(self, run_label, temp_tol):
    #   """Return the largest temperature found in diagnostic .out files with a preference for
    #      files found in the work directory"""
    #   import os
    #   import numpy as np
    #   from supercomputer import Supercomputer, TermColors
    #  
    #   stg_dir, wrk_dir = self._stage_dir, self._scratch_dir
 
    #   #Check work directory
    #   if(os.path.isfile(wrk_dir + '/' + run_label + '/subchandra_temp_diag.out')):
    #      diag_file = wrk_dir + '/' + run_label + '/subchandra_temp_diag.out'
    #      temps = np.loadtxt(diag_file, usecols=(1,), unpack=True)
    #      maxt = max(temps)
    #      if(maxt > temp_tol): 
    #         return TermColors.START_RED + '{:6.2f}'.format(maxt/1.e6) + TermColors.RESET
    #      else:
    #         return '{:6.2f}'.format(maxt/1.e6) 
    #   #If no luck, check staging directory
    #   elif(os.path.isfile(stg_dir + '/' + run_label + '/output/subchandra_temp_diag.out')):
    #      diag_file = stg_dir + '/' + run_label + '/output/subchandra_temp_diag.out'
    #      temps = np.loadtxt(diag_file, usecols=(1,), unpack=True)
    #      maxt = max(temps)
    #      if(maxt > temp_tol): 
    #         return TermColors.START_RED + '{:6.2f}'.format(maxt/1.e6) + TermColors.RESET
    #      else:
    #         return '{:6.2f}'.format(maxt/1.e6) 
    #   else:
    #      return 'no files'
    #
    #def _getMPeak(self, run_label, mach_tol):
    #   """Return the largest Mach number found in diagnostic .out files with a preference for
    #      files found in the work directory"""
    #   import os
    #   import numpy as np
    #   from supercomputer import TermColors
    #   
    #   stg_dir, wrk_dir = self._stage_dir, self._scratch_dir
    #   
    #   #Check work directory
    #   if(os.path.isfile(wrk_dir + '/' + run_label + '/subchandra_vel_diag.out')):
    #      diag_file = wrk_dir + '/' + run_label + '/subchandra_vel_diag.out'
    #      machnums = np.loadtxt(diag_file, usecols=(3,), unpack=True)
    #      maxm = max(machnums)
    #      if(maxm > mach_tol): 
    #         return TermColors.START_RED + '{:5.3f}'.format(maxm) + TermColors.RESET
    #      else:
    #         return '{:5.3f}'.format(maxm)
    #   #If no luck, check staging directory
    #   elif(os.path.isfile(stg_dir + '/' + run_label + '/output/subchandra_vel_diag.out')):
    #      diag_file = stg_dir + '/' + run_label + '/output/subchandra_vel_diag.out'
    #      machnums = np.loadtxt(diag_file, usecols=(3,), unpack=True)
    #      maxm = max(machnums)
    #      if(maxm > mach_tol): 
    #         return TermColors.START_RED + '{:5.3f}'.format(maxm) + TermColors.RESET
    #      else:
    #         return '{:5.3f}'.format(maxm)
    #   else:
    #      return 'no files'
    #
    #def _getTime(self, run_label):
    #   """Return the latest time found in diagnostic .out files with a preference for
    #      files found in the work directory"""
    #   import os
    #   
    #   stg_dir, wrk_dir = self._stage_dir, self._scratch_dir
 
    #   #Check work directory
    #   diag_file = None
    #   if(os.path.isfile(wrk_dir + '/' + run_label + '/subchandra_temp_diag.out')):
    #      diag_file = wrk_dir + '/' + run_label + '/subchandra_temp_diag.out'
    #   #If no luck, check staging directory
    #   elif(os.path.isfile(stg_dir + '/' + run_label + '/output/subchandra_temp_diag.out')):
    #      diag_file = stg_dir + '/' + run_label + '/output/subchandra_temp_diag.out'
    #   else:
    #      return 'no files'
    #   
    #   with open(diag_file, 'r') as f:
    #      f.seek(-2,2)              #Jump to near the end of the file (second-to-last byte)
    #      while f.read(1) != '\n':  #Read to EOL
    #         f.seek(-2, 1)          #Jump back a bit
    #      last_line = f.readline()  #Read current line
    #   tokens = last_line.split()
    #   return '{:06.2f}'.format(float(tokens[0]))
    
    def _getNote(self, run_label):
        """Return the any short note found in the staging directory's out directory
           as note.txt."""
        import os
        
        stg_dir, wrk_dir = self._stage_dir, self._scratch_dir
 
        #Check for a note, return it if it exists
        if(os.path.isfile(stg_dir + '/' + run_label + '/output/note.txt')):
            note_file = stg_dir + '/' + run_label + '/output/note.txt'
            with open(note_file, 'r') as nf:
                for i, l in enumerate(nf):
                    assert (i == 0), 'Error! note.txt should have only one line!'
                    txt = l
            return txt.strip('\n')
        
        #Otherwise, return a blank
        return ' '
 
    def _getActiveSimDirs(self):
        """Return a list of directories for each of the active simulations in the grid."""
        from os import listdir
        ret = []
 
        for d in listdir(self._stage_base):
            #archived runs are put into an inactive directory and ignored
            #TODO Implement archive method
            #TODO ASSUMPTION that all directories other than `inactive` are
            #simulation directories.  Make sure this is consistent, or refactor
            if d == 'inactive':
                continue
            ret.append(d)
        
        return ret

    def _getActiveSims(self):
        """Return a list of Simulation objects representing the simulations in this grid.
        
        This should be implemented by subclasses, making use of their subclass
        of Simulation."""
       
        #TODO Any way to make a reasonable version of this in SimulationGrid?
        #     Maybe pass in the desired Simulation subclass?
        raise NotImplementedError("A subclass of SimulationGrid did not implement
        this method or you're directly instantiating SimulationGrid.  Either way,
        NO!")

# Simulation class to represent a specific simulation.
class Simulation(object):
    """Represents a particular simulation's configuration and model parameters
    as well as the different management and analysis actions you'd like to carry
    out on a simulation."""

    def __init__(self, label, base_dir):
        """Initializes a simulation object using an existing setup in a
        directory. To generate a new simulation, use this class's static factory
        methods.
        
        Arguments:
            label      --> label for this simulation, will also be name of dir
                           where it's stored
            base_dir   --> path to the base directory this simulation is stored in
        """
        self.label = label
        self._full_path = join(base_dir, label)
        self._initFromDir(join(base_dir, label))

    def _initFromDir(self, simdir):
        """Initialize the simulation from an existing directory containing
        configuration and output."""
        self._config = self._genConfig(simdir)
        self._output = self._genOutput(simdir)

    def _genConfig(self, simdir):
        """Generate a SimConfig object for this simulation based on existing
        configuration."""

        raise NotImplementedError("A subclass of Simulation did not implement
        this method or you're directly instantiating Simulation.  Either way,
        NO!")

    def _genOutput(self, simdir):
        """Generate a SimOutput object for this simulation based on existing
        configuration."""

        raise NotImplementedError("A subclass of Simulation did not implement
        this method or you're directly instantiating Simulation.  Either way,
        NO!")



class SimConfig(object):
    """Represents all of the configuration needed to specify a particular
    simulation.  This includes inputs files, any initial models, and the files
    needed to execute the simulation (binaries, batch scripts, etc)."""

    def __init__(self, simdir):
        """Constructs a SimConfig object using an existing configuration in the
        given directory."""
        self._initFromDir(simdir)

    def _initFromDir(self, simdir):
        """Initialize this object using an existing configuration.  Subclasses
        must define how to do this for their particular code."""

        raise NotImplementedError("A subclass of SimConfig did not implement
        this method or you're directly instantiating SimConfig.  Either way,
        NO!")

class SimOutput(object):
    """Represents the products of a simulation, such as checkpoint files, data
    files, diagnostic data, etc."""

    def __init__(self, simdir):
        """Constructs a SimOutput object using an existing configuration in the
        given directory."""
        self._initFromDir(simdir)

    def _initFromDir(self, simdir):
        """Initialize this object using an existing configuration.  Subclasses
        must define how to do this for their particular code."""

        raise NotImplementedError("A subclass of SimOutput did not implement
        this method or you're directly instantiating SimOutput.  Either way,
        NO!")



###############################
### Machine Management Code ###
###############################
# The code in this section is all about abstracting a machine, broadly defined
# as a computational system (e.g. supercomputer, laptop, maybe even mobile?).
# I've sectioned it off as I may one day break it off into its own module.  
# TODO If this gets hefty enough, separate into a module.

# Machine class to represent the machine and filesystem we're currently on.


# TODO Write a small driver here that tests functionality.
