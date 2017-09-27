# This code implements the various classes from the simmy framework to represent
# a sub-chandra grid of models.  The code is based in part of code I previously
# wrote to carry out the models described in
# http://adsabs.harvard.edu/abs/2016ApJ...827...84J

# Author:        Adam Jacobs
# Creation date: September 22, 2017

# Usage: load as a module
# Requirements
#    + Python 3
#    + Common scientific python tools: NumPy, matplotlib
#    TODO: fill in other requirements

#TODO Implement intermediate packages, e.g. for Maestro?



###########################################
### Global Imports, Data, and Constants ###
###########################################
import sys
from simmy import SimulationGrid, Machine, Simulation, SimConfig, SimOutput


# Global checks, assertions
if not sys.version_info >= (3,):
    #TODO if you can, perhaps make this work for 2.7+ and 3 with from __future__
    raise RuntimeError("subchandra requires Python 3")

class SubChandraGrid(SimulationGrid):
    """A grid of Maestro simulations of sub-Chandrasekhar mass CO WDs with
    helium shells."""

    def __init__(self, label, stage_base, scratch_base):
        """Construct the SubChandraGrid object.

        Arguments:
        label        --> string describing the grid
        stage_base   --> base directory where simulations will be staged
        scratch_base --> base directory for the scratch space
                         where runs will be executed
        """
        super(label, stage_base, scratch_base)

    def listSimulations(self):
        """Print a list of simulations in this grid.
        
        This will only list 'active' simulations that the user is still
        exploring.  Simulations that are no longer being explored can be archived
        so that they will not pollute the list.
        
        Details listed: label, if it's in scratch space, if it's in the queue.
        """
        from simmy import START_GREEN, START_RED, START_BLUE, RESET
        from simmy import Machine
        from os.path import isfile, isdir, join
        from glob import glob

        active_sims = self._getActiveSimDirs()
        curcomp = Machine.getCurrentMachine()

        #Define heading/formatting for simulation list
        heading = '{0:29s}|{1:14s}'.format('Label', 'In scratch?')
        list_format = '{0:29s}|{1:14s}'
        if curcomp.has_queue:
            heading = '{0:29s}|{1:14s}|{2:14s}'.format('Label', 'In scratch?', 'In queue?')
            list_format = '{0:29s}|{1:14s}|{2:14s}'
        yep    = START_GREEN + '{0:14s}'.format("Yes!")    + RESET
        nope   = START_RED   + '{0:14s}'.format("No!")     + RESET
        purged = START_BLUE  + '{0:14s}'.format("Purged!") + RESET
       
        print(heading)
        for s in active_sims:
            #Check scratch: is it there, not, or there but purged?
            simdir = join(self._scratch_base, s)
            sc_str = nope
            if isdir(simdir):
                sc_str = purged
                found_all_expected_files = (
                      len(glob( join(simdir, 'main.*') ))         > 0 and
                      len(glob( join(simdir, 'inputs*') ))        > 0 and
                      len(glob( join(simdir, 'helm_table.dat') )) > 0
                      )
                if found_all_expected_files:
                    sc_str = yep
         
            #Check if the simulation is queued and define summary string
            if curcomp.has_queue:
                #Check queue 
                #  ASSUMPTION: sim directory is same as queue label
                q_str = nope
                if curcomp.isQueued(s):
                    q_str = yep
                outstr = list_format.format(s, sc_str, q_str) 
            else:
                outstr = list_format.format(s, sc_str) 

            print(outstr)

    def _getActiveSims(self):
        """Return a list of Simulation objects representing the simulations in this grid."""
        ret = []

        for d in listdir(self._stage_base):
            if d == 'inactive':
                continue
            #TODO This uses Simulation, which should generally be subclassed.
            #Think about if we should do some fancy reflection of some sort to
            #use the proper subclass' static factory method
            newSimObject = SCSimulation.genFromDir(d)
            ret.append(newSimObject)

        return ret

class SCSimulation(Simulation):
    """Represents a particular sub-Chandra simulation."""

    def __init__(self, label, base_dir):
        """Construct an SCSimulation object using an existing simulation.
        
        Arguments:
            label      --> label for this simulation, will also be name of dir
                           where it's stored
            base_dir   --> Path to the base directory this simulation is stored in.
                           By convention this will also be the label for the
                           grid of models this simulation belongs to.
        """
        super(label, base_dir)

    def _genConfig(self, simdir):
        """Generate an  SCConfig object for this simulation based on existing
        configuration."""
        return SCConfig(simdir)

    def _genOutput(self, simdir):
        """Generate an SCOutput object for this simulation based on existing
        configuration."""
        return SCOutput(simdir)

class SCConfig(SimConfig):
    """Represents all of the configuration needed to specify a sub-Chandra
    simulation.  This includes inputs files, initial models, and the location of
    needed binary files."""

    def __init__(self, simdir):
        """Constructs an SCConfig object using an existing configuration in the
        given directory."""
        super(simdir)

    def _initFromDir(self, simdir):
        """Initialize this object using an existing configuration."""
        from os.path import isfile, isdir, join, basename

        #Should include
        #   + inputs file data
        #   + initial model data
        #   + location of template files?
        #   + job configuration
        #
        self._initFiles(simdir) 
        self._initInputsDict()
        self._initIMDict()
        self._initRunDict()

    def _initFiles(self, simdir):
        """Initialize variables with paths to inputs and config files."""
        from os.path import join, basename
        #NOTE Convention established here: 
        #     A simulation directory contains a `model` directory with inputs
        #     file in it.  inputs file looks like
        #         inputs<dim>d.<simulation label>
        #TODO Handle multiple inputs files in model dir?
        #TODO Use regex instead of glob, it's safer and better-defined

        #Input file that is passed to executable as arg
        inputs_list = glob(join(simdir, 'model', 'inputs*'))
        if len(inputs_list) > 1:
            raise NotImplementedError("Having multiple inputs files in model dir
                    not currently implemented")
        self._inputs_file = inputs_list[0] #NOTE full path
      
        #Parameters file describing initial model parameters
        params_list = glob(join(simdir, 'model', '_params*'))
        if len(params_list) > 1:
            raise NotImplementedError("Having multiple _params files in model dir
                    not currently implemented")
        self._params_file = params_list[0]

        #Initial model data, pointed to by inputs file
        hse_list = glob(join(simdir, 'model', 'sub_chandra.*.hse.*'))
        if len(hse_list) > 1:
            raise NotImplementedError("Having multiple hse initial model files in model dir
                    not currently implemented")
        extras_list = glob(join(simdir, 'model', 'sub_chandra.*.extras.*'))
        if len(hse_list) > 1:
            raise NotImplementedError("Having multiple extras initial model files in model dir
                    not currently implemented")
        self._imdata_files = (hse_list[0], extras_list[0])

        #Script executed to run the simulation,
        #will be submitted to supercomputer queue
        runscript_list = glob(join(simdir, 'model', '*.run'))
        if len(runscript_list) > 1:
            raise NotImplementedError("Having multiple run scripts in model dir
                    not currently implemented")
        self._runscript_file = runscript_list[0]

    def _initInputsDict(self):
        """Initialize a dictionary of inputs variables based on the inputs file.
        
        Creates self._inputs_dict
        """
        #An inputs file consists of definitions of the form "var = value".
        #Here we convert this into a dictionary that will allow easy programmatic
        #access to these variables.
        self._inputs_dict = {}
        with open(self._inputs_file, 'r') as f:
            for line in f:
                tokens = line.partition('=')
                if tokens[1]: #Only do anything if a '=' was found
                    key = tokens[0].strip()
                    strval = tokens[2].strip()
                    self._inputs_dict[key] = strval

    def _initIMDict(self):
        """Initialize a dictionary of initial model data, config.
        
        Creates self._im_dict
        """
        from numpy import loadtxt, array
        #TODO Make sure loadtxt is robust for things like blank lines, bad lines, etc
        self._im_dict = {}

        #Store initial model parameters
        #TODO I use this logic multiple times, move to helper function?
        with open(self._params_file, 'r') as f:
            for line in f:
                tokens = line.partition('=')
                if tokens[1]: #Only do anything if a '=' was found
                    key = tokens[0].strip()
                    strval = tokens[2].strip()
                    self._im_dict[key] = strval

        #Store initial model data
        hse_file = self._imdata_files[0]
        extras_file = self._imdata_files[1]
        rad, rho, temp, pressure, Xhe4, Xc12, Xo16, Xfe56 = loadtxt(
                hse_file, unpack=True)
        rad, cs, ent = loadtxt(extras_file, unpack=True)
        self._im_dict['radius'] = rad
        self._im_dict['density'] = rho
        self._im_dict['temperature'] = temp
        self._im_dict['pressure'] = pressure
        self._im_dict['soundspeed'] = cs
        self._im_dict['entropy'] = ent
        self._ihe4, self._ic12, self._io16 = 0, 1, 2
        self._im_dict['species'] = array([Xhe4, Xc12, Xo16])

    def _initRunDict(self):
        """Initialize a dictionary of key properties describing the runscript
        submitted to the supercomputer queue."""
        #TODO Generalize this to run on either supercomputer/cluster or local machine
        self._run_dict = {}
        self._run_dict['allocation_label'] = 
        self._run_dict['job_label'] = 
        self._run_dict['walltime'] = 
        self._run_dict['nodes'] = 
        self._run_dict['threads_per_node'] = 
        self._run_dict['process_script'] = 
        self._run_dict['exe_label'] = 

    def _inputsFile(self, sim_label, **kwargs):
        """Return a string representing the content of an inputs file populated
        with data from **kwargs.
        
        Data is inserted into a template based on the passed keyword args.
        If no kwargs are given or expected values are missing,
        default values defined in this method are used.

        Arguments:
            sim_label --> label for the simulation, no spaces (e.g. 10044-090-210-4lev-full-512)

        Keyword arguments:
            im_file   --> Name of the "hse" initial model file
            job_name  --> Description of simulation that will be in the
                          job_info/inputs file
            max_levs  --> The maximum number of levels the AMR will refine to
            coarse_res --> Resolution of the base (coarsest) mesh
            anelastic_cutoff --> Density below which Maestro's velocity
                                 constraint switches to the anelastic constraint instead of the fancy
                                 low Mach constraint.  This helps alleviate issues at the edge of the
                                 star as density plummets.
            octant    --> .true. or .false.
                          If true, model an octant of a star, a full star otherwise
            dim       --> Dimensionality of the problem (2D or 3D)
            physical_size --> Physical extent of the domain cube in cm
            plt_delta --> Periodic interval in seconds at which pltfiles should
                          be saved. A pltfile will be saved every plt_delta of
                          simulation time.
            miniplt_delta --> Same as plt_delta but for minipltfiles
            chk_int   --> Periodic interval in timesteps (integer! not seconds)
                          at which checkpoint files should be saved
        """
        #Define default dictionary values
        inputs_keywords = {}
        inputs_keywords['im_file'] = "sub_chandra.M_WD-1.00.M_He-0.045.hse.C.10240"
        inputs_keywords['job_name'] = "512^3 base grid, T_core = 10^7, T_base = 210 MK -- M_WD=1.0, M_He=0.04"
        inputs_keywords['max_levs'] = 4
        inputs_keywords['coarse_res'] = 512
        inputs_keywords['anelastic_cutoff'] = 64000.0
        inputs_keywords['octant'] = ".false."
        inputs_keywords['dim'] = 3
        inputs_keywords['physical_size'] = 1500000000.0
        inputs_keywords['plt_delta'] = 5.0
        inputs_keywords['miniplt_delta'] = 0.2
        inputs_keywords['chk_int'] = 10

        #Define derived dictionary values
        pltfile_base = sim_label + "_plt"
        inputs_keywords['pltfile_base'] = pltfile_base
        miniplt_base = sim_label + "_miniplt"
        inputs_keywords['miniplt_base'] = miniplt_base
        chkfile_base = sim_label + "_chk"
        inputs_keywords['chkfile_base'] = chkfile_base
       
        if inputs_keywords['octant'].lower().count('false') > 0:
            inputs_keywords['bc_lo'] = 12
            inputs_keywords['bc_hi'] = 12
        else:
            inputs_keywords['bc_lo'] = 13
            inputs_keywords['bc_hi'] = 12

        #Define the base template string
        inputs_template = """&PROBIN
         model_file = "{im_file:s}"
         drdxfac = 5
        
         job_name = "{job_name}"
        
         use_alt_energy_fix = T
         ppm_trace_forces = 0
        
         hg_bottom_solver = 4
         mg_bottom_solver = 4
         max_mg_bottom_nlevels = 2
        
         max_levs = {max_levs:d}
         regrid_int = 2
        
         n_cellx = {coarse_res:d}
         n_celly = {coarse_res:d}
         n_cellz = {coarse_res:d}
        
         stop_time = 30000.
        
         max_step = 100000000
        
         init_iter = 1
         init_divu_iter = 3
         do_initial_projection = T
        
         max_grid_size_1 = 32
         max_grid_size_2 = 64
         max_grid_size_3 = 128
        
         the_sfc_threshold = 32768
        
         anelastic_cutoff = {anelastic_cutoff:f}
         base_cutoff_density = 10000.0
         buoyancy_cutoff_factor = 2.d0
         sponge_center_density = {anelastic_cutoff:f}
         sponge_start_factor = 2.0
         sponge_kappa = 10.0d0
        
         spherical_in = 1
         octant = {octant:s}
         dm_in = {dim:d}
         do_sponge = .true.
        
         prob_hi_x = {physical_size:f}
         prob_hi_y = {physical_size:f}
         prob_hi_z = {physical_size:f}
        
         plot_base_name = "{pltfile_base:s}"
         plot_int = -1
         plot_deltat = {plt_delta:f}
        
         mini_plot_base_name = "{miniplt_base:s}"
         mini_plot_int = -1
         mini_plot_deltat = {miniplt_delta:f}
         mini_plot_var1 = "species"
         mini_plot_var2 = "velocity"
         mini_plot_var3 = "temperature"
         mini_plot_var4 = "radial_velocity"
        
         check_base_name = "{chkfile_base:s}"
         chk_int = {chk_int:d}
        
         cflfac = 0.7d0
         init_shrink = 0.1d0
         max_dt_growth = 1.1d0
         use_soundspeed_firstdt = T
         use_divu_firstdt = T
        
         bcx_lo = {bc_lo:d}
         bcx_hi = {bc_hi:d}
         bcy_lo = {bc_lo:d}
         bcy_hi = {bc_hi:d}
         bcz_lo = {bc_lo:d}
         bcz_hi = {bc_hi:d}
        
            verbose = 1
         mg_verbose = 1
         cg_verbose = 0
        
         do_burning = T
        
         enthalpy_pred_type = 1
         evolve_base_state = T
        
         dpdt_factor = 0.0d0
         species_pred_type = 1
         use_tfromp = T
        
         single_prec_plotfiles = T
        
         use_eos_coulomb = T
        
         plot_trac = F
         plot_base = T
        
         velpert_amplitude = 1.d5
         velpert_scale = 5.d7
         velpert_steep = 1.d7
        
        /"""

        #Trim leading spaces

        #Populate template with data from dictionary
        


class SCOutput(SimOutput):
    """Represents the products of a sub-Chandra, such as the diagnostics files,
    reduced data from pltfiles, output from the executable, etc..."""

    def __init__(self, simdir):
        """Constructs an SCOutput object using an existing configuration in the
        given directory."""
        super(simdir)

    def _initFromDir(self, simdir):
        """Initialize this object using an existing configuration."""

        raise NotImplementedError("A subclass of SimOutput did not implement
        this method or you're directly instantiating SimOutput.  Either way,
        NO!")

#TODO Make driver for tests here (or make an independent test driver?)
#   + Verify input info
#   + Verify all object construction, generation works

