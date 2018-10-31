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
#   + For initial development I'm dumping most things here.  Once things are
#   reasonable prototyped, I need to break into modules.  Current module
#   categories I'm imagining are Simulation, Machine, and Util (e.g. for
#   TemplateFile, DefinedDict).
#   + Check PEP8 cromulency.
#   + For now, I'm using a strategy of copying needed files from the main
#   codebase into sim directories.  This means the files need to be built.  I
#   want to add functionality that builds executables and such instead of
#   copying manually built ones out.
#   + Add module docstring

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
    methods/properties provided by SimulationGrid.  SimulationGrid can also be
    seen as a simulation manager, so it's useful even if you're only running one
    model."""
    ##Global Data##
 
    ##Constructor##
    def __init__(self, label, code_base, stage_base, exe_base=None):
        """Initialize the SimulationGrid object.
 
           Arguments:
           self         --> implicitly passed reference to this instance of SimulationGrid
           label        --> this grid's label
           code_base    --> path to base directory containing source code the
                            simulations are built from
           stage_base   --> path to base directory containing all simulations in the grid
           exe_base     --> path to base directory where the simulations are executed (optional)
                            If not provided, it is assumed simulations are run in
                            the staging directories.

           Why a staging directory and an exe directory?
           On many machines there is a large filesystem (scratch/exe) for running
           simulations that is purged and a smaller filesystem (stage) that is
           backed up.  This is the motivation for having a staging and exe
           directory.
           
           The setup for simulations and their reduced (small size) output will
           be kept in staging, while the actual run and large output will be
           executed in exe.  Large output that can't fit in staging should be
           stored in high performance storage that is available on most machines
           that have the stage/scratch configuration.  Even for systems that
           don't have a purged scratch directory, this is a nice way to keep a
           clean space with simulation essentials and a place to tweak and
           experiment on a run.
           """
        self._code_base = code_base
        self._stage_base = stage_base
        if exe_base is None:
            self._exe_base = self._stage_base
            self._use_exe = False
        else:
            self._exe_base = exe_base
            self._use_exe = True
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
        
        if self._use_exe:
            heading = '{0:29s}|{1:14s}'.format('Label', 'In exe?')
            list_format = '{0:29s}|{1:14s}'
        else:
            heading = '{0:29s}'.format('Label')
            list_format = '{0:29s}'
        yep    = START_GREEN + '{0:14s}'.format("Yes!")    + RESET
        nope   = START_RED   + '{0:14s}'.format("No!")     + RESET
        purged = START_BLUE  + '{0:14s}'.format("Purged!") + RESET
        
        print(heading)
        for r in active_runs:
            if self._use_exe:
                #Check exe TODO mention purge check for subclasses?
                rundir = join(self._exe_base, r)
                if isdir(rundir):
                    if self._inEXE(rundir):
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
            print(s._getStatusSummary())
 
 
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
        
        stg_dir, wrk_dir = self._stage_dir, self._exe_dir
 
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
        raise NotImplementedError("""A subclass of SimulationGrid did not implement
        this method or you're directly instantiating SimulationGrid.  Either way,
        NO!""")

# Simulation class to represent a specific simulation.
class Simulation(object):
    """Represents a particular simulation's configuration and model parameters
    as well as the different management actions you'd like to carry out on a
    simulation.
    """

    def __init__(self, label, base_dir):
        """Initializes a simulation object using an existing setup in a
        directory. To generate a new simulation, use this class's static factory
        methods.
        
        Arguments:
            label      --> label for this simulation, will be used for naming
                           (e.g. directories and files)
            base_dir   --> path to the base directory this simulation is stored in
        """
        self.label = label
        self._full_path = join(base_dir, label)
        self._initFromDir(join(base_dir, label))

    def _initFromDir(self, simdir):
        """Initialize the simulation from an existing directory containing
        configuration and output."""
        self._simconfig = self._genSimConfigFromDir(simdir)
        self._runconfigs = self._genRunConfigsFromDir(simdir)
        self._output = self._genOutputFromDir(simdir)

    def _genSimConfigFromDir(self, simdir):
        """Generate a SimConfig object for this simulation based on existing
        configuration."""

        raise NotImplementedError("""A subclass of Simulation did not implement
        this method or you're directly instantiating Simulation.  Either way,
        NO!""")

    def _genRunConfigsFromDir(self, simdir):
        """Generate a list of RunConfig objects for this simulation based on existing
        configuration."""

        raise NotImplementedError("""A subclass of Simulation did not implement
        this method or you're directly instantiating Simulation.  Either way,
        NO!""")

    def _genOutputFromDir(self, simdir):
        """Generate a SimOutput object for this simulation based on existing
        configuration."""

        raise NotImplementedError("""A subclass of Simulation did not implement
        this method or you're directly instantiating Simulation.  Either way,
        NO!""")

    @classmethod
    def buildSim(cls, label, base_dir, **kwargs):
        """Build a new Simulation object configured with kwargs.  The kwargs
        will be determined by subclasses."""

        return cls._buildMe(label, base_dir, **kwargs)

    @classmethod
    def _buildMe(cls, label, base_dir, **kwargs):
        """Build a new Simulation object configured with kwargs.  The kwargs
        will be determined by subclasses."""

        raise NotImplementedError("""A subclass of Simulation did not implement
        this method or you're directly instantiating Simulation.  Either way,
        NO!""")



class SimConfig(object):
    """Represents all of the configuration needed to specify a particular
    simulation.  This could include inputs files, initial models, and/or the
    location of any relevant codebases."""
    #TODO 
    #   + Add utility method to inspect and explain config dictionaries
    #   + Add ability to modify an existing configuration

    def __init__(self, simdir, config_recs=None):
        """Constructs a SimConfig object representing a simulation found in
        simdir.
        
        If a simulation already exists in simdir, it will be used to create this
        class.  If not, you need to provide config_recs.  The simulation's
        label will be the name of the base directory of simdir.
        
        config_recs is list of ConfigRecords.  Each ConfigRecord specifies
        an aspect of the simulation's configuration.
        """
        from os.path import basename
        self._label = basename(simdir)
        self._simdir = simdir
        if config_recs is None:
            self._config_recs = self._initFromDir()
        else:
            self._config_recs = config_recs
            self._initFromRecs()

    def printConfig(self):
        """Print out the current configuration."""
        #This is a basic implementation.  Subclasses may want to override.
        for rec in self._config_recs:
            print(rec)

    def _initFromDir(self):
        """Initialize this object using an existing configuration found in
        self._simdir.  Subclasses must define how to do this for their
        particular code. A list of ConfigRecords defining the configuration should
        be returned."""

        raise NotImplementedError("""A subclass of SimConfig did not implement
        this method or you're directly instantiating SimConfig.  Either way,
        NO!""")

    def _initFromRecs(self):
        """Initialize this object using the ConfigRecords in self._config_recs.
        Subclasses must define how to do this for their particular code.
        """

        raise NotImplementedError("""A subclass of SimConfig did not implement
        this method or you're directly instantiating SimConfig.  Either way,
        NO!""")

class SimOutput(object):
    """Represents the products of a simulation, such as checkpoint files, data
    files, diagnostic data, etc, as well as tools for managing the data."""

    def __init__(self, simdir):
        """Constructs a SimOutput object using an existing configuration in the
        given directory."""
        self._initFromDir(simdir)

    def _initFromDir(self, simdir):
        """Initialize this object using an existing configuration.  Subclasses
        must define how to do this for their particular code."""

        raise NotImplementedError("""A subclass of SimOutput did not implement
        this method or you're directly instantiating SimOutput.  Either way,
        NO!""")

class ConfigRecord(object):
    """A single "configuration record" for a simulation.
    
    For example, a simulation could have a file that configures its initial
    model, a file that configures its parameters, and a set of commands to carry
    out the simulation.  Each of these could be represented with a ConfigRecord.

    The core of a ConfigRecord is a dictionary with a set of fixed keys which
    must be defined.
    
    A ConfigRecord may also have an associated file that they know how to
    convert into a new ConfigRecord as well as how to write based on an existing
    ConfigRecord.
    """
    #TODO Initially I thought of having this implement MutableMapping, but in
    #   trying I wasn't sure it was buying me what I wanted and wasn't sure how to
    #   safely customized the parsing of *args and **kwargs sent to __init__.  Might
    #   want to revisit this.

    #Design: 
    #   + limited dictionary consisting of a set of keys defined at init.  keys
    #     can't be changed after this.
    #   + parallel dictionary with description of keys
    #   + optional associated file
    #TODO Should ConfigRecord be able to initialize self from file, or leave
    #   that to callers?

    def __init__(self, fields_dict, label, description, fieldmap=None):
        """Initialize ConfigRecord with a dict of the valid fields for this
        record as well as a description of each field.
       
        Arguments:
            fields_dict: The keys of fields_dict will be the valid fields, and the values will be
                         the descriptions.
            label:       String label for this configuration.
            description: String description of this configuration.

        Keyword Arguments:
            fieldmap: An optional alias map that maps
                      aliases to valid fields.  The fieldmap is useful for mapping variable
                      names used in files to field names.
                      fieldmap --> {'file_variable': 'field'}
        """
        self._label = label
        self._desc = description
        self._fields = tuple(fields_dict.keys())
        self._fields_desc = fields_dict.copy()
        self._config_dict = {key: None for key in fields_dict}
        self._myfile = None
        self._fieldmap = fieldmap

    @classmethod
    def genFromFile(cls):
        """Construct a ConfigRecord from file."""
        #TODO Does this make sense?
        pass

    def associateFile(self, tempfile):
        """Associate a TemplateFile with this ConfigRecord."""
        #TODO Verify all fields are accounted for
        #TODO Verify tempfile is TemplateFile?
        #TODO In practice you always need this ConfigRecord's config_dict to
        #create the tempfile that gets passed in.  So we should probably move the
        #logic of creating the tempfile into here instead of having tempfile
        #passed.
        #TODO Also, in practice the tempfile's replacement dict has to always
        #match (modulo any fieldmap) ConfigRecord's config_dict.  As of now, we
        #rely on the user making sure this is true, but should instead build it
        #into the design of ConfigRecord + TemplateFile.
        self._myfile = tempfile

    def setFields(self, **kwargs):
        """Set fields found in **kwargs."""
        for key, val in kwargs.iteritems():
            self.setField(key, val)

    def setField(self, field, data):
        """Set a single field."""
        #Map the field, if needed
        cur_key = field
        if self._fieldmap is not None and cur_key in self._fieldmap:
            cur_key = self._fieldmap[field]
        if cur_key not in self._fields:
            raise KeyError("{} is not a valid field for this ConfigRecord".format(cur_key))
        self._config_dict[cur_key] = data

    def getField(self, field):
        """Get the data stored in field."""
        if field not in self._fields:
            raise KeyError("{} is not a valid field for this ConfigRecord".format(field))
        return self._config_dict[field]

    def getFieldDict(self):
        """Get the underlying field dictionary."""
        #TODO Changes to this by user will be reflected here, should do deep copy?
        return self._config_dict

    def saveFile(self, savepath):
        """Save a configured file from this ConfigRecord's data."""
        self._myfile.saveFile(savepath)

    def __str__(self):
        ret = self._label + '\n'
        ret += 'Description: {}\n'.format(self._desc)
        ret += 'Fields:\n'
        for field, desc in self._fields_desc.items():
            ret += '    {}\n        Description: {}\n        Current value: {}\n'.format(field, desc, self._config_dict[field])
        return ret


###############################
### Machine Management Code ###
###############################
# The code in this section is all about abstracting a machine, broadly defined
# as a computational system (e.g. supercomputer, laptop, maybe even mobile?).
# I've sectioned it off as I may one day break it off into its own module.  
# TODO If this gets hefty enough, separate into a module.
# TODO Originally I designed this as an ABC meta class. Now I just mark abstract
# methods with NotImplementedError.  Is there any payoff for implementing as ABC
# worth the extra complication?

# Machine class to represent the machine and filesystem we're currently on.
class Machine(object):
    """
    Abstracts a computational machine, from a laptop to a supercomputer.
    
    The abstraction is specifically designed for simulation management.  The
    intent is to enable the design of scripts for a problem/simulation that can
    be run with little-to-no modification from one machine to another.  Ideally,
    it would be "no" modification, but it's hard to completely hide the
    hardware/software/filesystem differences between, say, a laptop and
    leadership-class computing facility.
    """

    def __init__(self):
        """
        Initialize the Machine object.
        
        The basic design as of now is to define a bunch of base DefinedDicts
        that subclasses then populate with data.  The bases serve as Machine's
        declaration of the required data all subclasses must implement.
        Subclasses are free to add their own supplemental DefinedDicts with data
        particular to that subclass.
        """
        #+Meta data
        #TODO Some systems have lots of possible hostnames (e.g. iCER). I may
        #   want to change hosts to contain regexes or globs, or perhaps change how
        #   I use the hostname to determine where Machine instance to return with
        #   the factor method
        meta_desc = {'name':  "str, A name for this Machine",
                     'hosts': "list, List of strings of valid hostnames for this Machine"}
        meta_dict = {'name': None,
                     'hosts': None}
        #This base metadata is required to be defined by all Machine
        # implementations.  Implementations are free to add additional fields to
        # their own meta ddicts.  
        #If the base ddict defaults to None, it is required to be defined by
        # implementations. If it defaults to a "null" value ("" for strings,
        # NaN for numbers, etc), then it is optional.
        self._meta_base = DefinedDict(meta_dict, meta_desc)

        #iCER example
        #meta_dict = {'name': 'iCER',
        #             'hosts': ['dev-intel14',
        #                       'dev-intel14-k20',
        #                       'dev-intel14-phi',
        #                       'dev-intel16',
        #                       'dev-intel16-k80',
        #                       'dev-intel18']}

        #+Computational resources available
        self._partitions = {} #dictionary of {'partition_label': partition DefinedDict}
                              #for all partitions on the system.
        cpart_desc = {'arch': 'str, informal description of the computing ' +
                              'architecture or processor and any attached ' +
                              'accelerators, e.g. KNL, haswell, POWER9 + '  +
                              'Volta V100, i7, Xeon E5 + K80, etc',
                      'node_count': 'int, number of nodes in the partition',
                      'hw_cores_pn': 'int, hardware/physical cores per node',
                      'logical_core_fac': 'int, number of logical cores per ' +
                                          'hardware core (e.g. hardware ' +
                                          'threads per physical core)',
                      'gpus_pn': 'int, GPUs per node',
                      'mem_domains_pn': 'int, memory domains per node (e.g. ' +
                                        'NUMA nodes per node)',
                      'mem_per_domain': 'str, human-readable size of memory ' +
                                         'in each domain, e.g. 256 GB'}
        cpart_dict = {'arch': '',
                      'node_count': None,
                      'hw_cores_pn': None,
                      'logical_core_fac': None,
                      'gpus_pn': None,
                      'mem_domains_pn': None,
                      'mem_per_domain': None}
        #TODO: Should these base definitions be put elsewhere? Not sure they
        #make sense in __init__
        self._partition_base = DefinedDict(cpart_dict, cpart_desc)

        #+Storage resources available
        #Local filesystem, includes scratch but NOT HPSS
        filesys_desc = {'user_root': 'str, Full path to the user root ' +
                                     'directory, e.g. $HOME on most *nix systems.',
                        'scratch_root': 'str, Full path to the base directory ' +
                                        'you would like to use in any scratch filesystem.  If ' +
                                        'None, it will be assumed no scratch is available.'}
        filesys_dict = {'user_root': None, 'scratch_root': None}
        self._filesys_base = DefinedDict(filesys_dict, filesys_desc)

        #HPSS
        #TODO Implement HPSS functionality

        #+Define what's required for generating a batch submission script.  The
        #   idea is for this to be generic, and subclasses can then implement SLURM,
        #   PBS, or whatever the machine uses.
        #TODO: need to figure out what to do re: a bunch of batch scripts that
        #   run one executable and scripts that run multiple executables in one job.
        #   For prototyping, stick with single executable per batch script
        #   (though I'm allowing use of array jobs, since these are very close
        #   to single executable per batch and can make things easier)
        #TODO: SLURM takes ntasks, easiest PBS equivalent is nodes.  Need to
        #   figure out how to abstract this to work for both, if possible.  For now,
        #   just let ntasks=nodes
        batch_desc = {'job_name': 'str, Name to be used for the job, no spaces!',
                      'walltime': 'str, The requested walltime, in form [DD:]HH:MM:SS',
                      'array_str': 'str, optional, If you want an array job, give a valid array str here, e.g. "1-24", "1,4,8-9". "" if you do not want an array job (default)',
                      'nodes': 'int, Number of nodes requested',
                      'tasks_pn': 'int, Number of tasks (roughly, processes or MPI tasks) per node, defaults to 1',
                      'cores_pt': 'int, Number of cores (cpus) per task',
                      'mem_pn': 'str, human-readable memory requested per node (e.g. 20 GB), defaults to most memory guaranteed to be available (system-dependent)',
                      'user_email': 'str, optional, email to send any alerts to',
                      'exe_script': 'str, filename for the script to be executed in the batch submission, put all simulation-specific execution instructions here.'}
        batch_dict = {'job_name': None,
                      'walltime': None,
                      'array_str': "",
                      'nodes': None,
                      'tasks_pn': 1,
                      'mem_pn': None,
                      'user_email': "",
                      'exe_script': None}
        
        self._batch_base = DefinedDict(batch_dict, batch_desc)

    def genBatch(self, batch_path, batch_ddict, exe_text, batch_template=None):
        """
        Generate a batch script.

        Arguments:
            batch_path     --> str, Full path to the file the generated batch script
                               will be saved to
            batch_ddict    --> DefinedDict of the data needed for a batch
                               submission.  Use static method getBaseBatchDDict() to get a base
                               ddict.
            exe_text       --> str, executable text to append to the batch
                               script.  These are the simulation/problem-specific shell
                               statements that will be executed.
            batch_template --> TemplateFile, optional, the TemplateFile to be
                               used to generate the batch script.  By default will be generated
                               with static method getBatchTemplateFile().
        """

        #TODO Add batch_base consistency check/error check
        #TODO User needs to build and send batch_ddict, not subclass.  Should
        #   make machinery to facilitate this.
        raise NotImplementedError("""A subclass of Machine did not implement
        this method or you're directly instantiating Machine.  Either way,
        NO!""")


    def _getHome(self):
        """Get the root directory for user's home on this machine."""
        #TODO Error checking, verify exists, etc 
        return self._filesys_base['user_root']
        #raise NotImplementedError("""A subclass of Machine did not implement
        #this method or you're directly instantiating Machine.  Either way,
        #NO!""")

    def _getScratch(self):
        """
        Get the root directory for scratch space on this machine.
        
        Scratch space is a part of the filesystem or an independent filesystem
        reserved for large amounts of data.  Such spaces are often purged of old
        files at regular intervals.  If this machine has no scratch space,
        returns None.
        """
       
        return self._filesys_base['scratch_root']
        #raise NotImplementedError("""A subclass of Machine did not implement
        #this method or you're directly instantiating Machine.  Either way,
        #NO!""")

    @staticmethod
    def getBatchTemplateFile():
        """
        Return the base TemplateFile for this Machine.
        
        Subclasses must implement this.
        """
        
        raise NotImplementedError("""A subclass of Machine did not implement
        this method or you're directly instantiating Machine.  Either way,
        NO!""")

    @staticmethod
    def getBaseBatchDDict():
        """
        Return the base DefinedDict for a batch job submission.
        
        Subclasses are free to provide their own expanded batch ddict, but must
        include all fields given here.
        """
        batch_desc = {'job_name': 'str, Name to be used for the job, no spaces!',
                      'walltime': 'str, The requested walltime, in form [DD:]HH:MM:SS',
                      'array_str': 'str, optional, If you want an array job, give a valid array str here, e.g. "1-24", "1,4,8-9". "" if you do not want an array job (default)',
                      'nodes': 'int, Number of nodes requested',
                      'tasks_pn': 'int, Number of tasks (roughly, processes or MPI tasks) per node, defaults to 1',
                      'cores_pt': 'int, Number of cores (cpus) per task',
                      'mem_pn': 'str, human-readable memory requested per node (e.g. 20 GB), defaults to most memory guaranteed to be available (system-dependent)',
                      'user_email': 'str, optional, email to send any alerts to',
                      'exe_script': 'str, filename for the script to be executed in the batch submission, put all simulation-specific execution instructions here.'}
        batch_dict = {'job_name': None,
                      'walltime': None,
                      'array_str': "",
                      'nodes': None,
                      'tasks_pn': 1,
                      'mem_pn': None,
                      'user_email': "",
                      'exe_script': None}
        
        return DefinedDict(batch_dict, batch_desc)

    @staticmethod
    def getCurrentMachine():
        """Return a Machine representing the current host machine.
        
        This requires that the detected system or host has an implementation
        in the python path.
        """
        import socket
        from platform import system
    
        #Based on the hostname, import an implementation of Machine
        #TODO: I like the idea here, but it'll clearly get unwieldy once enough
        #   machines are implemented.  Need to design it to do the reflection better.
        #   AFTER prototyping is done!
        host = socket.gethostname()
        if host is None:
            sys_type = system().lower()
            if sys_type.startswith('linux'):
                machmodule = __import__("linux_generic") # The __import__ builtin is
                                                         # a way to import modules not named until runtime
                return machmodule.LinuxMachine(scomp_config)
            elif not sys_type:
                raise OSError("Cannot determine platform.")
            else:
                raise OSError("Unimplemented system platform({}) with no hostname information.".format(sys_type))
        elif host.count('titan') == 1:
            try:
                machmodule = __import__("titan")        # The __import__ builtin is
                                                        # a way to import modules not named until runtime
                return machmodule.TitanSC(scomp_config)
            except ImportError:
                print('Cannot find module for host {}'.format(host))
                return None
        elif host.count('dev-intel14') == 1: #TODO Add all the other iCER hosts, dev-*, gateway, etc
            try:
                machmodule = __import__("icer")         # The __import__ builtin is
                                                        # a way to import modules not named until runtime
                return machmodule.ICER(scomp_config)
            except ImportError:
                print('Cannot find module for host {}'.format(host))
                return None
        #elif host.count('sn.astro') == 1:
        #    try:
        #        machmodule = __import__("sn")
        #        return machmodule.snSC(scomp_config)
        #    except ImportError:
        #        print('Cannot find module for host {}'.format(host))
        #        return None
        #elif host.count('siona') == 1:
        #    try:
        #        machmodule = __import__("siona")
        #        return machmodule.sionaSC(scomp_config)
        #    except ImportError:
        #        print('Cannot find module for host {}'.format(host))
        #        return None
        else:
            raise NotImplementedError("Unimplemented host: {}".format(host))
    

#TODO Get rid of this, just use ConfigRecord as generic config object.
class RunConfig(object):
    """Represents the configuration and files needed to execute a simulation on
    a particular machine.  
    
    An instance of this class is stored in the Machine class.  Subclasses should
    implement this for a particular machine (or class of machine, e.g. a generic
    linux cluster) to be run on.
    """

    def __init__(self, simdir, config_dict=None):
        """Constructs a RunConfig object for a simulation found in
        simdir.
        
        If a simulation already exists in simdir, it will be used to create this
        class.  If not, you need to provide config_dict.  The run's
        label will be the name of the base directory of simdir.
        
        config_dict is a dictionary of all parameters needed to configure a run
        script for the simulation.  These may be scripts that are run directly
        or submitted to a job queue.  Properties such as resources requires
        (cores, memory), that name of the executable and any arguments, etc, may
        be stored in the dictionary.
        """
        from os.path import basename
        self._label = basename(simdir)
        self._simdir = simdir
        if config_dict is None:
            self._initFromDir()
        else:
            self._config_dict = config_dict
            self._initFromDict()

    def _initFromDir(self):
        """Initialize this object using an existing configuration found in
        self._simdir.  Subclasses must define how to do this for their
        particular code. A dictionary defining the configuration should
        be initialized as self._config_dict."""

        #TODO Fragile to rely on implementers to define self._config_dict. Could
        #redesign or use ABS?
        raise NotImplementedError("""A subclass of RunConfig did not implement
        this method or you're directly instantiating RunConfig.  Either way,
        NO!""")

    def _initFromDict(self):
        """Initialize this object using the configuration dictionary found in
        self._config_dict.  Subclasses must define how to do this for their
        particular machine.
        """

        raise NotImplementedError("""A subclass of RunConfig did not implement
        this method or you're directly instantiating RunConfig.  Either way,
        NO!""")

####################
### Utility Code ###
####################
# The code in this section supports other code sections/modules, but is not
# specific to them.  Like with other sections, it may be broken off into its own
# module.
from collections.abc import Mapping, MutableMapping

## DefinedDict
# This is my stab at a class I've been wanted for a long while: DefinedDict.
#
# DefinedDict is a restricted variation on Python's dict.  It has a defined set
# of keys/fields that is established on creation and cannot be changed
# afterward.  However, the data the keys point to CAN be updated.  I subclass
# Mapping instead of MutableMapping because MM suggests the ability to add,
# change, remove keys, which I don't want to do.
#
# In developing this and discussing with others, I've come across possible
# alternatives to writing my own class.  None quite gave me what I wanted, but
# for reference check out: `attrs`, `recordclass`, and the new builtin
# dataclasses.  The issues with these include: not quite achieving what I want,
# introducing dependencies without enough of a payoff, and/or requiring as much
# effort to achieve what I want as writing DefinedDict.  In particular, I value
# how a DefinedDict 
#   + restricts the fields(keys), allowing one to write a template DefinedDict
#     that guarantees a certain set of data to be available when passed to other
#     code aware of the template,
#   + gives fields little docstring-like descriptions, and
#   + is lightweight and dict-like.
#
# In some ways, I see these as facilitating the rapid prototyping of and
# experimenting with potential classes, classes that may have interesting
# inheritance relationships.
#
# TODO: use that __hasattr__ to make mydefdict['key'] and mydefdict.key both
# work?
#
# TODO: __slots__ is a magic that signals a fixed set of attributes will be
# used, instead of the usual __dict__ that allows for adding an arbitrary
# number of attributes and implies a decent amount of overhead.  If this class
# proves a good case, might use this instead to be more memory efficient.
# However, Mapping may already implement this for us.
#
# Python's documentation and source provides useful info on implementing
# containers, though sadly not a single StackOverflow answer pointed to these
# docs.  I stumbled upon them after much frustration.
# See S 3.3.7 of:
# https://docs.python.org/3/reference/datamodel.html#objects
# In addition, see the source of UserDict, Mapping, MutableMapping. _Environ in
# os module is an example of a MutableMapping implementation.
#
# Recommended methods to implement (* provided by MutableMapping (** abstract)):
#   **keys(), *values(), *items(), *get(), clear(), setdefault(), pop(), popitem(),
#   copy(), and update()
#   __*contains__(), **__iter__(), **__len__(), **__getitem__(),__missing__(),
#   **__setitem__(), **__delitem__(), __repr__(), __str__(), __format__()
#
# Some details on some methods:
#object.__getitem__(self, key)
#Called to implement evaluation of self[key]. [negative keys handled by __getitem__()]. 
#   If key is of an inappropriate type, TypeError may be raised; if of a value outside the set of indexes for the sequence (after any special interpretation of negative values), IndexError should be raised. 
#   For mapping types, if key is missing (not in the container), KeyError should be raised.
#
#   Note for loops expect that an IndexError will be raised for illegal indexes to allow proper detection of the end of the sequence.
#
#object.__missing__(self, key)
#Called by dict.__getitem__() to implement self[key] for dict subclasses when key is not in the dictionary.
#
#object.__setitem__(self, key, value)
#Called to implement assignment to self[key]. Same note as for __getitem__(). 
#This should only be implemented for mappings if the objects support changes to the values for keys, or if new keys can be added, or for sequences if elements can be replaced. 
#The same exceptions should be raised for improper key values as for the __getitem__() method.
#
#object.__delitem__(self, key)
#Called to implement deletion of self[key]. Same note as for __getitem__(). 
#This should only be implemented for mappings if the objects support removal of keys, or for sequences if elements can be removed from the sequence. 
#The same exceptions should be raised for improper key values as for the __getitem__() method.
#
#object.__iter__(self)
#This method is called when an iterator is required for a container. This method should return a new iterator object that can iterate over all the objects in the container. For mappings, it should iterate over the keys of the container.
#
#The membership test operators (in and not in) are normally implemented as an iteration through a sequence. However, container objects can supply the following special method with a more efficient implementation, which also does not require the object be a sequence.
#
#object.__contains__(self, item)
#Called to implement membership test operators. Should return true if item is in self, false otherwise. 
#For mapping objects, this should consider the keys of the mapping rather than the values or the key-item pairs.
#
#For objects that don’t define __contains__(), the membership test first tries iteration via __iter__(), then the old sequence iteration protocol via __getitem__(), see this section in the language reference.
#
# 
# NOTE: Implementation strategy, start with what Mapping gives, take what you
# need from MutableMapping, UserDict
#
# TODO: Implement type hinting for dict items?
# TODO: Enforce key to be string?
# TODO: Implement __add__ so that you can do addition similar to tuples.  Like
#       with tuples, you shouldn't be able to change an initialized
#       DefinedDict's fields, but it would be nice to facilitate building new
#       ones out of existing.
class DefinedDict(Mapping):
    """
    A dict-like object which has a restricted set of keys/fields defined on
    creation, each of which includes a human-readable description of the key
    similar to a docstring.
    """

    def __init__(self, init_dict, desc_map):
        """
        DefinedDict(D, desc)
            D: dictionary to initialize DefinedDict with.
            desc: mapping of D's keys to human-readable strings describing each key

        Example:
            D = {'field1': 1, 'field2': 2}
            desc = {'field1': 'The first field, an integer.',
                    'field2': 'The second field, an integer.'}
            mydd = DefinedDict(D,desc)

            f1 = mydd['field1']  # Works
            mydd['field3'] = f1  # Raises KeyError
            mydd['field2'] = f1  # Works
            mydd.desc('field2')  # Prints 'The second field, an integer.'
        """
        #NOTE: key private attributes are:
        #self._data = {} The underlying data dictionary
        #self._desc = {} A mapping of keys from _data to descriptions
        #self._fields = frozenset() Immutable set of the keys for both desc and data.
        #   This is included because the keys() sets for either _data or _desc
        #   can be changed without raising errors. This set provides an
        #   immutable reference for ~O(1) lookups to verify fields.  __init__
        #   guarantees _fields == _desc.keys() == _data.keys()
        self._data = {}

        #Define the fields, must be done before __setitem__ can work.
        if isinstance(desc_map, Mapping):
            for key in desc_map:
                if key not in init_dict:
                    error_message = "keys of init_dict and desc_map don't match"
                    raise KeyError(error_message)
            #TODO: Check that desc_map is all strings
            self._fields = frozenset(desc_map.keys())
            if not (self._fields == set(desc_map.keys()) == set(init_dict.keys())):
                raise ValueError('init_dict and desc_map do not have the same keys!')
            #Should be no need for a deep copy
            self._desc = desc_map.copy()
        else:
            error_message = "desc_map is not a Mapping: {}".format(desc_map)
            raise TypeError(error_message)

        #Initialize the local dict, borrowing practices from MutableMapping update()
        if isinstance(init_dict, Mapping):
            for key in init_dict:
                self[key] = init_dict[key]
        elif hasattr(init_dict, "keys"):
            for key in init_dict.keys():
                self[key] = init_dict[key]
        else:
            for key, value in init_dict:
                self[key] = value

    def desc(self, key):
        """Get a string describing `key`."""
        if key not in self._fields:
            raise KeyError(key)
        return self._desc[key]

    def fields_str(self):
        """Get a string listing and describing all fields."""
        return self.__str__(verbose=False)

    #TODO: setdefault(), copy(), deepcopy()

    ## magic/dunder methods
    ## the key methods required by Mapping:
    # __getitem__ __iter__, __len__
    def __getitem__(self, key):
        if key in self._fields:
            #TODO: this will raise KeyError if somehow key not in _data even
            #though it was in self._fields, should do a try except?
            return self._data[key]
        else:
            raise KeyError(
                  '`{}` is not a valid field in this DefinedDict'.format(key))

    def __iter__(self):
        return iter(self._data)

    def __len__(self):
        return len(self._data)

    ## abstract MutableMapping methods and dict methods to consider
    #  __setitem__ __delitem__ __add__ __eq__ __ne__ keys() 
    #  __contains__ (esp if I ever need to deal with __missing__ protocols, see UserDict)
    def __setitem__(self, key, value):
        if key in self._fields:
            self._data[key] = value
        else:
            raise KeyError(
                  '"{}" is not a valid field in this DefinedDict'.format(key))

    def __delitem__(self, key):
        raise TypeError('DefinedDict fields (keys) are immutable, cannot be deleted')

    def __contains__(self, key):
        return self._data.__contains__(key)

    ## object magics
    # Potentially useful class dunders
    #__init,str,repr,getattr,setattr,delattr()__
    def __str__(self, verbose=True):
        field_str = '{} - {}'
        join_str = '\n'
        if verbose:
            #If verbose, include data and separate fields with ---
            field_str += '\n    {}'
            join_str += '---\n'
        ret_str = "DefinedDict class instance\n{}".format(
            join_str.join(
                (field_str.format(field, self._desc[field], self._data[field])
                for field in self._desc)))
        return ret_str

    #TODO: __repr__, which as much as possible should be valid Python code that
    #   could be used to recreate the instance.

## TemplateFile
# TODO: determine dependencies/version reqs
#   Uses OrderedDict, which means 2.7+ without more work.  Note for 3.7+ all
#   dicts are ordered.

# NOTE: While mature Python Templating engines exist, when I did a quick review
# the general convention seems to have a set of data (e.g. in a dict) that is
# substituted into a text template (rendered).  That is half of what I want to
# do.  I also want to be able to read in a file that has already been rendered
# and use it to populate the data.  This may exist, but it wasn't clear how to
# do this in less time with Python builtins or mature templating engines.  May
# want to revisit this after prototyping.

# TODO: I've just made a basic class, but I could subclass the builtin File
# type.  After prototyping, check if this offers benefits.

# TODO: How should line regexes work on lines with more than one sub descriptor?
#   As of now, it is assumed all descriptors on the same line have the same line
#   regex

# TODO: I kinda abuse Formatter in this class.  Maybe worth it to write a
# Formatter subclass?

# TODO: Allow users to give a format-spec for each field for rendering

class TemplateFile(object):
    """
    A TemplateFile implements a very common practice of having a template of a
    file's text that is then populated with some user data and rendered, or a
    rendered file can be read in to extract the data.
    """

    def __init__(self, template_text):
        """Construct a TemplateFile

        Arguments:
            template_text    --> A template of the file's content, with
                                 substitution descriptors of the form 
                                 {field[line_regex]:field_regex} (see below).

        Substitution descriptors:
            The template_text can contain substitution descriptors used to both
            render a file and to interpret a rendered file. They are of the form
                {field[line_regex]:field_regex}
            `field` is the name of the variable/field that stores the
            substitution data, `field_regex` is a regular expression that is
            guaranteed to match the form of the substitution data, `line_regex`
            is an optional regex that will match the line the field is in.  If
            not given, all non-empty lines are matched (i.e. line_regex = `.+`)

        Rendering a file:
            A new rendering of the TemplateFile can be written after
            initializing data for all fields.  This is done with one of
            initData(), setField(), or passing the data dict directly to
            render() (see each method for details).

        Interpreting a file:
            To interpret a file, use the class factory method fromfile() to make
            a new initialized TemplateFile.
            For an existing, initialized TemplateFile `tfile`,
            use `tfile.readFile()`

        Example of creating and rendering new file:
            template_text = ""\"
            # Comment describing {filename[.*#.*]:[a-zA-Z0-9]+.dat}
            {resvar:resolution} = {resdata:[0-9]{2,4}}
            ""\"

            mytfile = TemplateFile(template_text)

            init_dict = {'filename': 'myfile.dat',
                'resvar' = 'resolution',
                'resdata' = '128'}

            mytfile.initData(init_dict)

        Example of loading in a file:
            template_text = ""\"
            # Comment describing {filename[.*#.*]:[a-zA-Z0-9]+.dat}
            {resvar:resolution} = {resdata:[0-9]{2,4}}
            ""\"

            # Assuming test.txt lines are 
            # "# Comment describing myfile.dat"
            # "resolution = 128",
            # then
            mytfile = fromfile('test.txt', template_text)
            # will yield a data dict with
            print(mytfile.getField('filename')) # myfile.dat
            print(mytfile.getField('resvar'))   # resolution
            print(mytfile.getField('resdata'))  # 128
        """
        #Key instance variables:
        #   self._template_text: The file's template text
        #   self._fields: Frozen set of the file fields.
        #   self._regexes: ordered mapping of _fields to (field_regex, line_regex)
        #Initialized here, but initially populated with None data:
        #   self._data: mapping of _fields to substitution data
        self._template_text = template_text
        self._fields, self._regexes = self._getFieldsAndRegs()
        self._data = {}
        for field in self._fields:
            self._data[field] = None

    def writeFile(self, savepath):
        """Write the file to the given full path, making directories if needed."""
        from os.path import dirname
        from os import makedirs
        if dirname(savepath):
            makedirs(dirname(savepath), exist_ok=True)
        with open(savepath, 'w') as f:
            f.write(self.render())

    def setField(self, field, value):
        """Set `field` to store `value`."""
        if field not in self._fields:
            raise KeyError('{} is not a valid field for this TemplateFile!'.format(field))
        self._data[field] = value

    def getField(self, field):
        """Get the data stored in `field`."""
        if field not in self._fields:
            raise KeyError('{} is not a valid field for this TemplateFile!'.format(field))
        return self._data[field]

    def initData(self, data_dict):
        """Initialize the dictionary mapping fields to substitution data."""
        if data_dict.keys != self._fields:
            raise KeyError("data_dict keys do not match this TemplateFile's fields!")
        self._data = data_dict

    def render(self, data_dict=None):
        """Render the template into usable file text, return the text."""
        if self._data is None and data_dict is None:
            raise RuntimeError('TemplateFile data has not been initialized!')
        rendered_text = ""
        for line in self._template_text.splitlines(keepends=True):
            rendered_line = self._renderLine(line)
            rendered_text += rendered_line
        return rendered_text

    def _renderLine(self, line):
        """Render a single line."""
        if line.count('{') == 0:
            #No substitution descriptors, just render as is
            return line
        if line.count('}') < 1:
            error_message = "TemplateFile: Template text contains invalid line!\n"
            error_message += "    Offending line: {}".format(line)
            raise ValueError(error_message)

        #Convert any substitution descriptors into valid format strings
        remove = False
        open_count = 0
        format_line = ""
        for i,c in enumerate(line):
            if c == '{':
                open_count += 1
            if c == '}':
                if open_count-1 == 0:
                    remove = False
                open_count -= 1
            if open_count > 0 and not remove:
                #Look for [ or :
                if c == '[' or c == ':':
                    remove = True
            if not remove:
                format_line += c

        return format_line.format(**self._data)

    def readFile(self, filepath):
        """Initialize data from file at `filepath`."""
        import re
        #NOTE: This algorithm makes use of the fact that fields are found in
        #   order, from the top to the bottom of the file.  It assumes data for
        #   all fields is available in the file.

        #If no fields, just return
        if len(self._fields) == 0:
            return

        #Build a stack to pop fields off of. TODO: is this inefficient?
        field_stack = []
        for field in self._regexes:
            field_regex, line_regex = self._regexes[field]
            field_stack.append((field, field_regex, line_regex))
        field_stack.reverse()

        #Iterate over lines in the file, searching for fields in order
        with open(filepath, 'r') as f:
            cur_field, cur_fregex, cur_lregex = field_stack.pop() #(field, field_regex, line_regex)
            fre = re.compile(cur_fregex)
            lre = re.compile(cur_lregex)
            found_all = False
            for line in f:
                lmatch = lre.match(line)
                if lmatch:
                    fmatch = fre.search(line)
                    #We use a while here to search for multiple fields on a
                    #single line. ASSUMPTION: all fields on a line have the
                    #same line_regex, if given
                    while fmatch:
                        field_data = line[fmatch.start():fmatch.end()]
                        self.setField(cur_field, field_data)
                        if len(field_stack) == 0:
                            found_all = True
                            break
                        cur_field, cur_fregex, cur_lregex = field_stack.pop()
                        fre = re.compile(cur_fregex)
                        lre = re.compile(cur_lregex)
                        fmatch = fre.search(line)
                if found_all:
                    break
       
        if not found_all:
            #We didn't find the current item, so add it back to the stack
            field_stack.insert(0,(cur_field,cur_fregex, cur_lregex))
            err_str = """All fields not found! Remaining fields in search order:
            {}""".format(field_stack)
            raise RuntimeError(err_str)

    @classmethod
    def fromfile(cls, filepath, template_text):
        """Factory method for creating a TemplateFile from an existing file."""
        ret = cls(template_text)
        ret.readFile(filepath)
        return ret

    def _getFieldsAndRegs(self):
        """Return fields and regexes in self._template_text.
        
        Parses self._template_text and returns a frozenset of fields as well as
        a mapping of those fields to their field regex and any line regex
            {<field_str>: (<field_regex>, <line_regex | None>)}"""
        #TODO Above docstring is out of date
        from string import Formatter
        from collections import OrderedDict
        #TODO: Error checking
        #TODO: What happens if field_regex contains :?
        fmt = Formatter()
        field_ret = []
        field_index = 1
        fregex_index = 2
        regex_ret = OrderedDict({})
        for parse_tuple in fmt.parse(self._template_text):
            field_text = parse_tuple[field_index]
            if field_text is None:
                continue
            #TODO: Currently, line_regex defaults to any non-blank line (.+).  I think it
            #   would be better to default to any literal text preceding the first sub
            #   descriptor
            line_regex = ".+"
            #Get the optional line_regex if provided
            #TODO maybe make sure user knows a literal '[' in field will break things
            if field_text.count('[') > 0:
                tokens = field_text.partition('[')
                field_text = tokens[0].strip()
                line_regex = tokens[2].rpartition(']')[0]
            field_ret.append(field_text)
            field_regex = parse_tuple[fregex_index]
            regex_ret[field_text] = (field_regex, line_regex)

        return frozenset(field_ret), regex_ret

###################
### Test Driver ###
###################
#
# TODO Write a small driver here that tests functionality.
