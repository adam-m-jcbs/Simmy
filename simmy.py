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
#   + Would it be better to use Namespaces instead of dicts for config?
#   + For initial development I'm dumping most things here.  Once things are
#   reasonable prototyped, I need to break into modules.  Current module
#   categories I'm imagining are Simulation, Machine, and Util (e.g. for
#   TemplateFile, ConfigRecord).
#   + Check PEP8 cromulency.
#   + For now, I'm using a strategy of copying needed files from the main
#   codebase into sim directories.  This means the files need to be built.  I
#   want to add functionality that builds executables and such instead of
#   copying manually built ones out.

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

class TemplateFile(object):
    """Represents a template file. 
    
    The essential elements are a string representing the content of a file with
    .format()-style format-spec replacement tokens through the file as well as a
    dictionary containing the values to be substituted in.
    """
    #TODO Add a regex dict of some sort to recognize data lines TemplateFiles
    #want to interact with.  This will facilitate reading data in instead of just
    #writing it as we currently do.

    def __init__(self, replacement_dict, template_string, lead_space):
        """Construct a TemplateFile

        Arguments:
            replacement_dict --> Dictionary of data to be substituted into the template file.
            template_string --> A template of the file's content, with
                               .format()-style format-spec tokens to be replaced with
                               replacement_dict values
            lead_space --> The number of leading spaces to be removed from all
                           lines in template_string except for the first line (similar to how
                           python doc strings are processed).
        """
        self._replacement_dict = replacement_dict
        self._template_string = template_string
        self._lead_space = lead_space
        self._initFileText()

    def saveFile(self, savepath):
        """Save the file to the given full path."""
        from os.path import dirname
        from os import makedirs
        makedirs(dirname(savepath), exist_ok=True)
        with open(savepath, 'w') as f:
            f.write(self._filetext)

    def getFileText(self):
        """Get the final, processed text of the file."""
        return self._filetext

    def _initFileText(self):
        """Initialize the file's text."""
        template_lines = self._template_string.splitlines(keepends=True)
        file_lines = [template_lines[0]]
        #Trim leading spaces, except on the first line
        for line in template_lines[1:]:
            file_lines.append(line[self._lead_space:])
        filetext = "".join(file_lines)

        #Populate template with data from dictionary
        #print('rep dict:   ', self._replacement_dict)
        #print('rep string: ', filetext)
        self._filetext = filetext.format(**self._replacement_dict)
 
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

# Machine class to represent the machine and filesystem we're currently on.
class Machine(object):
    """Represents a computational machine, from a laptop to a supercomputer."""

    def __init__(self):
        """Initialize the Machine object."""
        self.home_root = self._getHome()
        self.scratch_root = self._getScratch()

    def _getHome(self):
        """Get the root directory for user's home on this machine."""
        
        raise NotImplementedError("""A subclass of Machine did not implement
        this method or you're directly instantiating Machine.  Either way,
        NO!""")

    def _getScratch(self):
        """Get the root directory for scratch space on this machine.
        
        Scratch space is a part of the filesystem or an independent filesystem
        reserved for large amounts of data.  Such spaces are often purged of old
        files at regular intervals.  If this machine has no scratch space,
        returns None."""
        
        raise NotImplementedError("""A subclass of Machine did not implement
        this method or you're directly instantiating Machine.  Either way,
        NO!""")

    @staticmethod
    def getCurrentMachine():
        """Return a Machine representing the current host machine.
        
        This requires that the detected system or host has an implementation
        in the python path.
        """
        import socket
        from platform import system
    
        #Based on the hostname, import an implementation of Machine
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




# TODO Write a small driver here that tests functionality.
