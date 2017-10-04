# This code implements the various classes from the simmy framework to represent
# a sub-chandra grid of models.  The code is based in part on code I previously
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
    needed binary files.
    
    The class tries to handle the details and set reasonable defaults for things
    that don't change much from simulation to simulation.  The user-tunable
    properties of the configuration are stored in a set of dictionaries."""

    def __init__(self, simdir, config_dicts=None):
        """Constructs an SCConfig object using an existing configuration in the
        given directory."""
        super().__init__(simdir, config_dicts)

    def _initFromDir(self, simdir):
        """Initialize this object using an existing configuration."""
        from os.path import isfile, isdir, join, basename

        #Should include
        #   + inputs file data
        #   + initial model data
        #   + location of template files?
        #   + job configuration
        #
        #TODO Currently I have all dictionary values from files as strings (even
        #   things that are naturally floats, ints, or boolean).
        #   This makes it easy to pull them out of and put them back in file form.
        #   If proves annoying, change.
        self._config_dicts = {
                {'inputs_dict': {}},
                {'im_dict': {}}
                }
        self._initFiles(simdir) 
        self._initInputsDictFromDir()
        self._initIMDictFromDir()

    def _initFiles(self, simdir):
        """Initialize variables with paths to inputs and config files.
        
        Initializes self._inputs_file, self._params_file, self._imdata_files
        """
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

    def _initInputsDictFromDir(self):
        """Initialize a dictionary of inputs variables based on the inputs file.
        
        Populates self._config_dicts['inputs_dict']
        """
        #An inputs file consists of definitions of the form "var = value".
        #Here we convert this into a dictionary that will allow easy programmatic
        #access to these variables.
        #TODO As of now, this isn't consistent with init from scratch.  I use
        #   more human-readable keys like 'coarse_res' instead of n_cellx =
        #   n_celly = n_cellz.  But this will read in n_cell*.
        inputs_dict = self._config_dicts['inputs_dict']
        with open(self._inputs_file, 'r') as f:
            for line in f:
                tokens = line.partition('=')
                if tokens[1]: #Only do anything if a '=' was found
                    key = tokens[0].strip()
                    strval = tokens[2].strip()
                    inputs_dict[key] = strval

    def _initIMDictFromDir(self):
        """Initialize a dictionary of initial model data, config from the files
        found in the simulation directory.
        
        Populates self._config_dicts['im_dict'].  Contains all of the data from
        initial model files and the _params file used to generate this data.
        """
        from numpy import loadtxt, array
        #TODO Make sure loadtxt is robust for things like blank lines, bad lines, etc
        im_dict = self._config_dicts['im_dict']

        #Store initial model parameters
        #TODO I use this logic multiple times, move to helper function?
        with open(self._params_file, 'r') as f:
            for line in f:
                tokens = line.partition('=')
                if tokens[1]: #Only do anything if a '=' was found
                    key = tokens[0].strip()
                    strval = tokens[2].strip()
                    im_dict[key] = strval

        #Store initial model data
        hse_file = self._imdata_files[0]
        extras_file = self._imdata_files[1]
        rad, rho, temp, pressure, Xhe4, Xc12, Xo16, Xfe56 = loadtxt(
                hse_file, unpack=True)
        rad, cs, ent = loadtxt(extras_file, unpack=True)
        im_dict['radius'] = rad
        im_dict['density'] = rho
        im_dict['temperature'] = temp
        im_dict['pressure'] = pressure
        im_dict['soundspeed'] = cs
        im_dict['entropy'] = ent
        self._ihe4, self._ic12, self._io16 = 0, 1, 2
        im_dict['species'] = array([Xhe4, Xc12, Xo16])

    def _initFromDicts(self):
        """Initialize this object using the configuration dictionaries found in
        self._config_dicts.  The files associated with these dictionaries will
        be created.

        For sub-Chandra, the config dicts are labeled "im_dict" and "inputs_dict"

        im_dict:
            A dictionary of parameters used to generate the 1D initial model as
            well as the data from the generated model.
        inputs_dict:
            A dictionary of the variables in the inputs file that's passed to
            the Maestro executable.

        Note that these dictionaries need not be fully specified by the user and
        that some entries may be derived.
        """



        #Finally, fully initialize all dictionaries from the created files.
        #TODO I like to do this to make sure things are consistent for the two
        #   methods of initialization.  But maybe it's not needed?
        self._initInputsDictFromDir()

    def _initInputsDict(self):
        """Initialize the inputs dictionary of key properties describing the
        inputs parameters passed to the Maestro executable.
       
        Initialize self._config_dicts['inputs_dict'] with default values for any
        keys missing in the current dictionary.
        """
        #TODO As of now, users can let all values be default.  This doesn't make
        #   sense, choose which values are required to be passed.
        #Define default dictionary values
        inputs_dict = self._config_dicts['inputs_dict']
        inputs_defaults = {}
        inputs_defaults['im_file'] = "sub_chandra.M_WD-1.00.M_He-0.045.hse.C.10240"
        inputs_defaults['job_name'] = "512^3 base grid, T_core = 10^7, T_base = 210 MK -- M_WD=1.0, M_He=0.04"
        inputs_defaults['max_levs'] = '4'
        inputs_defaults['coarse_res'] = '512'
        inputs_defaults['anelastic_cutoff'] = '64000.0'
        inputs_defaults['octant'] = ".false."
        inputs_defaults['dim'] = '3'
        inputs_defaults['physical_size'] = '1500000000.0'
        inputs_defaults['plt_delta'] = '5.0'
        inputs_defaults['miniplt_delta'] = '0.2'
        inputs_defaults['chk_int'] = '10'

        for key in inputs_defaults:
            if key not in inputs_dict:
                inputs_dict[key] = inputs_default[key]

        #Define derived dictionary values
        pltfile_base = self._label + "_plt"
        inputs_dict['pltfile_base'] = pltfile_base
        miniplt_base = self._label + "_miniplt"
        inputs_dict['miniplt_base'] = miniplt_base
        chkfile_base = self._label + "_chk"
        inputs_dict['chkfile_base'] = chkfile_base
       
        if inputs_dict['octant'].lower().count('false') > 0:
            inputs_dict['bc_lo'] = '12'
            inputs_dict['bc_hi'] = '12'
        else:
            inputs_dict['bc_lo'] = '13'
            inputs_dict['bc_hi'] = '12'

    def _genInputsFile(self, savepath):
        """Generate and save an inputs file populated with data
        self._config_dicts['inputs_dict'].
        
        Data is inserted into a template based on the inputs dictionary, which
        should already be initialized before calling.
        If expected values are missing, default values defined in this method
        are used.

        Arguments:
            savepath  --> Full path to the location where this file should be saved

        inputs_dict keys used.  Those with * are typically unique to a
        simulation and should have been initialized in the dictionary:
            im_file*   --> Name of the "hse" initial model file
            job_name*  --> Description of simulation that will be in the
                           job_info/inputs file
            max_levs   --> The maximum number of levels the AMR will refine to
            coarse_res --> Resolution of the base (coarsest) mesh
            anelastic_cutoff* --> Density below which Maestro's velocity
                                  constraint switches to the anelastic constraint instead of the fancy
                                  low Mach constraint.  This helps alleviate issues at the edge of the
                                  star as density plummets.
            octant     --> .true. or .false.
                           If true, model an octant of a star, a full star otherwise
            dim        --> Dimensionality of the problem (2D or 3D)
            physical_size* --> Physical extent of the domain cube in cm
            plt_delta  --> Periodic interval in seconds at which pltfiles should
                           be saved. A pltfile will be saved every plt_delta of
                           simulation time.
            miniplt_delta --> Same as plt_delta but for minipltfiles
            chk_int    --> Periodic interval in timesteps (integer! not seconds)
                           at which checkpoint files should be saved
        """
        from simmy import TemplateFile
        #Define the base template string
        inputs_template = """&PROBIN
         model_file = "{im_file:s}"
         drdxfac = 5
        
         job_name = "{job_name:s}"
        
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

        #Now create and save the file
        inputs_file = TemplateFile(inputs_dict, inputs_template, 8)
        inputs_file.saveFile(self._inputs_file)

    def _initIMDict(self, gen_im=True):
        """Initialize the initial model dictionary of key properties describing the
        1D initial model for this sub-Chandra simulation.
       
        Initialize self._config_dicts['im_dict'] with default values for any
        keys missing in the current dictionary.  This will also run the 1d model
        generator using the given parameters, unless gen_im is False.
        """
        #TODO As of now, users can let all values be default.  This doesn't make
        #   sense, choose which values are required to be passed.
        #Define default dictionary values for the '_params' part of the IM dict.
        #This is needed for the 1d model generator
        im_dict = self._config_dicts['im_dict']
        im_defaults = {}
        im_defaults['M_tot'] = '1.0'                   #Mass of WD core in M_sol
        im_defaults['M_He']  = '0.0445'                #Mass of He envelope in M_sol
        im_defaults['delta'] = '2000000.0'             #Transition delta in cm
        im_defaults['temp_core'] = '90000000.0'        #Isothermal core temp in K
        im_defaults['temp_base'] = '210000000.0'       #Temp at base of He env in K
        im_defaults['mixed_co_wd'] = '.false.'         #Core is C/O or just C?
        im_defaults['low_density_cutoff'] = '1.d-4'    #Density floor for init model (not 3D Maestro simulation)
        im_defaults['temp_fluff'] = '7.5d7'            #Temp floor, temp also set to this when density floor hit
        im_defaults['smallt'] = '1.d6'                 #As far as I can tell, #this isn't used.  Maybe was synonym for temp_fluff

        #The initial model resolution should match Maestro's base state
        #resolution.  This is derived from inputs.  TODO Users are allowed to
        #override this, but shouldn't?
        indict = self._config_dicts['inputs_dict']
        max_levs = float(indict['max_levs'])
        drdxfac = float(indict['drdxfac'])
        coarse_res = float(indict['drdxfac'])
        fine_res = coarse_res*2*(max_levs-1)
        basestate_res = fine_res*drdxfac
        nx = basestate_res
        im_defaults['nx'] = str(nx)  #Resolution of the 1D model, number of cells

        #The physical size of initial model also can be derived from inputs.
        #for octant, it is same as domain size.  For full star, half.
        #TODO In practice, an initial model is tried to get an idea of the
        #   physical size of the domain, which is then put into inputs.
        #   However, below we get IM size from inputs.  Would be
        #   nice to formalize this algorithm here instead of the current method of
        #   doing it manually.  
        #   The basic algorithm is to do an initial model with the desired
        #   properties with a relatively huge xmax.  Then, redo the initial
        #   model with xmax set to be the radius of T_peak + 50% of that radius.
        #   This will also be the size of the domain of the 3D grid.
        #   This gives reasonable balance of buffer zone between surface of star
        #   and edge of domain without wasting too many resources on unimportant
        #   parts of the domain.
        octant =
        if octant:
            xmax = domain_size
        else:
            xmax = domain_size/2.0
        
        im_defaults['xmin'] = '0.0'
        im_defaults['xmax'] = str(xmax)

        for key in im_defaults:
            if key not in im_dict:
                im_dict[key] = im_defaults[key]

        #TODO Maybe separate out 1d model generation?  But, IM dict can't be
        #   fully init'd without 1d model file.
        if gen_im:
            #Use _params to generate the 1D model.

            #Add results to dict

        else:
            #IM will be incomplete! Blank out uninitializeable fields.


        im_dict['radius'] = rad
        im_dict['density'] = rho
        im_dict['temperature'] = temp
        im_dict['pressure'] = pressure
        im_dict['soundspeed'] = cs
        im_dict['entropy'] = ent
        self._ihe4, self._ic12, self._io16 = 0, 1, 2
        im_dict['species'] = array([Xhe4, Xc12, Xo16])

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

