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
    properties of the configuration are stored in ConfigRecords."""

    #TODO Organize methods: public first, static last, etc
    def __init__(self, simdir, config_recs=None):
        """Constructs an SCConfig object using an existing configuration in the
        given directory."""
        super().__init__(simdir, config_recs)

    def _initFromDir(self):
        """Initialize this object using an existing configuration."""

        #Should include
        #   + inputs file data
        #   + initial model data
        #   + Xlocation of template files?
        #   + Xjob configuration
        #
        #TODO Should fields all be strings or have appropriate type?
        #   Strings makes it easy to pull them out of and put them back in file form.
        config_recs = []
        self._initFilesFromDir() 
        #TODO For now I'm having SCConfig keep a local reference to the
        #   ConfigRecords.  This suggests to me I might want to have subclasses of
        #   ConfigRecord like InputsRecord and IMRecord.  I like the SimConfig super
        #   class being able to do generic operations by iterating over
        #   ConfigRecords, so if I do subclass I want to be careful to maintain this
        #   ability.
        self._inputs_rec = self._initInputsRecFromDir()
        self._im_rec     = self._initIMRecFromDir()
        config_recs.append(self._inputs_rec)
        config_recs.append(self._im_rec)
        return config_recs

    def _initFromRecs(self):
        """Initialize this object using the partially initialized ConfigRecords
        in self._config_recs.

        The minimum non-default field values that need to be defined for this are:

        Initial Model: M_tot, M_he, delta, temp_core, temp_base
        Everything else can be derived.
        """
        #TODO Other things can be derived as mentioned above, but maybe I want
        #to make sure they can be customized without derivation overwriting?
        #Design
        #  + Use the given parameters to build a first attempt at initial model
        #  with large radius
        #  + Based on result, rebuild initial model with a more reasonable
        #  radius
        #  + Store initial model and fully initialize im record
        #  + Use results to fully initialize inputs record
        #
        #  Reference _computeRmacAndCutoffs, _generateInitModel from subchandra.py
        self._initFiles()
        self._buildIM()
        #Fully initialize im record
        #With im record in hand, fully initialize inputs record

    def _buildIM(self):
        """Use the partially initialized self._im_rec to build an initial model.

        Fields needed from self._im_rec.
        Fields with * can reasonably be left to the default:
            + M_tot  = Mass of the WD core in M_sol.
            + M_He   = Mass of He envelope in M_sol.
            + temp_base = Temperature at the base of the He envelope in K.
            + *delta = Transition delta from core to envelope in cm.
            + *temp_core = Isothermal core temperature in K.

        Fields needed from self._inputs_rec:
            + max_levs = Number of levels of refinement.
            + coarse_res = Resolution of the base (coarsest) grid.
            + drdxfac  = Factor by which finest grid's resolution is multiplied
                         to get the base state resolution in spherical geometry.
            + octant   = Boolean, .true. means we model an octant, not the full star.
        """
        from subprocess import call, Popen, PIPE, STDOUT
        from os.path import join, isfile
        from os import remove
        from glob import glob
        from shlex import split
        #Design
        #  + Use the given parameters to build a first attempt at initial model
        #  with large radius
        #  + Based on result, rebuild initial model with a more reasonable
        #  radius

        #Calculate base state resolution, which should also be the initial
        #model's resolution
        max_levs = int(self._inputs_rec.getField('max_levs'))
        coarse_res = int(self._inputs_rec.getField('coarse_res'))
        fine_res = coarse_res*2**(max_levs-1) 
        drdxfac = int(self._inputs_rec.getField('drdxfac'))
        octant = self._inputs_rec.getField('octant')
        octant = octant.lower().count('true') > 0
        if octant:
            base_state_res = drdxfac*fine_res
        else:
            base_state_res = drdxfac*fine_res/2
        #TODO Here I'm forcing the correct im resolution, but users may expect
        #     any nx they pass to be used.  Should design to make it clear users
        #     can't set nx.
        self._im_rec.setField('nx', str(base_state_res))

        #For now, choose a pretty huge size.  We'll adjust down later.
        self._im_rec.setField('xmax', '1.1e9')

        #We should now have all the information we need to write a _params file
        #for the initial model builder.
        im_template_text, lead_space = SCConfig._getIMTempText()
        field_dict = im_rec.getFieldDict()
        im_tempfile = TemplateFile(field_dict, im_template_text, lead_space)
        self._im_rec.associateFile(im_tempfile)
        self._im_rec.saveFile(self._params_file)

        #Execute the initial model builder
        #TODO RESTART HERE
        #Make sure helmtable is linked
        if not isfile('helm_table.dat'):
            call(['ln', '-s', join(simconfig.Maestro_home, 'Microphysics', 'EOS', 'helmeos', 'helm_table.dat')])

        #Build the executable command
        init1d_exe = join(simconfig.Maestro_home, 'Util', 'initial_models', 'sub_chandra', 
              'init_1d.Linux.gfortran.debug.exe') + ' ' + pfilename

        #Execute, removing any old IM data files
        old_files = glob('sub_chandra.M_WD*')
        for f in old_files:
            remove(f)
        i1d_proc = Popen(split(init1d_exe), stdout=PIPE, stderr=PIPE)
        (i1d_out, i1d_err) = i1d_proc.communicate()
        if i1d_err:
            print('init1d error: ', i1d_err)
            print('init1d out, last 5 lines:')
            for line in i1d_out.split('\n')[-5:]:
                print(line)
 

    def _initFilesFromDir(self):
        """Initialize variables with paths to inputs and config files for an
        existing simulation.
        
        Initializes self._inputs_file, self._params_file, self._imdata_files
        """
        from os.path import join, basename
        from glob import glob
        #NOTE Convention established here: 
        #     A simulation directory contains a `model` directory with inputs
        #     file in it.  inputs file looks like
        #         inputs<dim>d.<simulation label>
        #TODO Handle multiple inputs files in model dir?
        #TODO Use regex instead of glob, it's safer and better-defined

        #Input file that is passed to executable as arg
        inputs_list = glob(join(self._simdir, 'model', 'inputs*'))
        if len(inputs_list) > 1:
            raise NotImplementedError("Having multiple inputs files in model dir not currently implemented")
        self._inputs_file = inputs_list[0] #NOTE full path
        #print('inputs file: {}'.format(self._inputs_file))
      
        #Parameters file describing initial model parameters
        params_list = glob(join(self._simdir, 'model', '_params*'))
        if len(params_list) > 1:
            raise NotImplementedError("Having multiple _params files in model dir not currently implemented")
        self._params_file = params_list[0]
        #print('_params file: {}'.format(self._params_file))

        #Initial model data, pointed to by inputs file
        hse_list = glob(join(self._simdir, 'model', 'sub_chandra.*.hse.*'))
        if len(hse_list) > 1:
            raise NotImplementedError("Having multiple hse initial model files in model dir not currently implemented")
        extras_list = glob(join(self._simdir, 'model', 'sub_chandra.*.extras.*'))
        if len(hse_list) > 1:
            raise NotImplementedError("Having multiple extras initial model files in model dir not currently implemented")
        self._imdata_files = (hse_list[0], extras_list[0])
        #print('im data files: {}'.format(self._imdata_files))

    def _initFiles(self):
        """Initialize variables with paths to inputs and config files to be
        written for a new simulation.
        
        Initializes self._inputs_file, self._params_file.
        self._imdata_files is initialized after the initial model is built.
        """
        from os.path import join, basename
        from glob import glob
        #NOTE Convention established here: 
        #     A simulation directory contains a `model` directory with inputs
        #     file in it.  inputs file looks like
        #         inputs<dim>d.<simulation label>
        #TODO Handle multiple inputs files in model dir?
        #TODO Use regex instead of glob, it's safer and better-defined

        #All config files for the model go here
        base_dir = join(self._simdir, 'model')
        model_label = self._label
       
        #Input file that is passed to executable as arg
        inputs_filename = "inputs3d.{}".format(model_label)
        self._inputs_file = join(base_dir, inputs_filename)
        #print('inputs file: {}'.format(self._inputs_file))
      
        #Parameters file describing initial model parameters
        params_filename = "_params.{}".format(model_label)
        self._params_file = join(base_dir, params_filename)
        #print('_params file: {}'.format(self._params_file))

    def _initInputsRecFromDir(self):
        """Initialize a ConfigRecord of inputs variables based on the inputs file."""
        from simmy import ConfigRecord, TemplateFile
        #An inputs file consists of definitions of the form "var = value".
        #Here we convert this into a ConfigRecord that will allow easy programmatic
        #access to and manipulation of these variables.

        #Get file variables
        file_vars = {}
        with open(self._inputs_file, 'r') as f:
            for line in f:
                tokens = line.partition('=')
                if tokens[1]: #Only do anything if a '=' was found
                    key = tokens[0].strip()
                    strval = tokens[2].strip()
                    file_vars[key] = strval

        #Define fields and initialize ConfigRecord
        inputs_rec = SCConfig.genInputsConfigRec()
        for key, val in file_vars.items():
            try:
                inputs_rec.setField(key, val)
            except KeyError:
                pass
                #print('{} is an extra key in the file'.format(key))

        #Define TemplateFile
        inputs_template_text, lead_space = SCConfig._getInputsTempText()
        field_dict = inputs_rec.getFieldDict()
        inputs_tempfile = TemplateFile(field_dict, inputs_template_text, lead_space)

        #Associate the file and return
        inputs_rec.associateFile(inputs_tempfile)
        return inputs_rec

    def _initIMRecFromDir(self):
        """Initialize an initial model ConfigRecord from the files
        found in the simulation directory.
        
        Returns a ConfigRecord representing initial model configuration.
        Contains all of the data from initial model files and the _params file
        used to generate this data.
        """
        from numpy import loadtxt, array
        from simmy import ConfigRecord, TemplateFile

        #Store initial model parameters from _params file
        #TODO I use this logic multiple times, move to helper function?
        file_vars = {}
        with open(self._params_file, 'r') as f:
            for line in f:
                tokens = line.partition('=')
                if tokens[1]: #Only do anything if a '=' was found
                    key = tokens[0].strip()
                    strval = tokens[2].strip()
                    file_vars[key] = strval

        #Define fields and initialize ConfigRecord
        im_rec = SCConfig.genIMConfigRec()
        for key, val in file_vars.items():
            try:
                im_rec.setField(key, val)
            except KeyError:
                pass
                #print('{} is an extra key in the file'.format(key))

        #Store initial model data
        hse_file = self._imdata_files[0]
        extras_file = self._imdata_files[1]
        #TODO Make sure loadtxt is robust for things like blank lines, bad lines, etc
        rad, rho, temp, pressure, Xhe4, Xc12, Xo16, Xfe56 = loadtxt(
                hse_file, unpack=True)
        rad, cs, ent = loadtxt(extras_file, unpack=True)
        im_rec.setField('radius',  rad)
        im_rec.setField('density', rho)
        im_rec.setField('temperature', temp)
        im_rec.setField('pressure', pressure)
        im_rec.setField('soundspeed', cs)
        im_rec.setField('entropy', ent)
        self._ihe4, self._ic12, self._io16 = 0, 1, 2
        im_rec.setField('species', array([Xhe4, Xc12, Xo16]))

        #Define TemplateFile
        im_template_text, lead_space = self._getIMTempText()
        field_dict = im_rec.getFieldDict()
        im_tempfile = TemplateFile(field_dict, im_template_text, lead_space)

        #Associate the file and return
        im_rec.associateFile(im_tempfile)
        return im_rec

    @staticmethod
    def genInputsConfigRec():
        """Return an inputs ConfigRecord with some default values set.
        
        This provides a baseline for users to fully initialize and then use to
        create new SCConfig objects from scratch.
        """
        from simmy import ConfigRecord, TemplateFile
        #Define fields and initialize ConfigRecord
        fields_dict, fieldmap = SCConfig._getInputsFields()
        rec_label = 'Inputs Configuration'
        rec_desc = """Configuration of the inputs file.  This is the file passed
        to the Maestro executable that sets various Maestro parameters,
        configures the simulation, and provides the location of initial model
        data."""
        inputs_rec = ConfigRecord(fields_dict, rec_label, rec_desc, fieldmap)
        inputs_defaults = SCConfig._getInputsDefaults()
        for key, val in inputs_defaults.items():
            inputs_rec.setField(key, val)

        return inputs_rec

    @staticmethod
    def genIMConfigRec():
        """Return an initial model ConfigRecord with some default values set.
        
        This provides a baseline for users to fully initialize and then use to
        create new SCConfig objects from scratch.
        """
        from simmy import ConfigRecord
        #Define fields and initialize ConfigRecord
        fields_dict = SCConfig._getIMFields()
        rec_label = "Initial Model Configuration"
        rec_desc  = """Configures the initial model for this simulation.  This
        corresponds to the _params file used by init1d to build an initial 1D
        model to be mapped to the 3D domain.  The data from this model are also
        stored."""
        im_rec = ConfigRecord(fields_dict, rec_label, rec_desc)
        im_defaults = SCConfig._getIMDefaults()
        for key, val in im_defaults.items():
            im_rec.setField(key, val)

        return im_rec

    @staticmethod
    def _getInputsDefaults():
        """Get a dictionary of default values for inputs fields."""
        #TODO I'm redundantly setting things that do not make sense to have a
        #default for to None.  Helps me keep track, but maybe should just delete.
        inputs_defaults = {}
        inputs_defaults['im_file'] = None
        inputs_defaults['job_name'] = None
        inputs_defaults['max_levs'] = '4'
        inputs_defaults['coarse_res'] = '512'
        inputs_defaults['anelastic_cutoff'] = None
        inputs_defaults['octant'] = ".false."
        inputs_defaults['dim'] = '3'
        inputs_defaults['physical_size'] = None
        inputs_defaults['plot_deltat'] = '5.0'
        inputs_defaults['mini_plot_deltat'] = '0.2'
        inputs_defaults['chk_int'] = '10'
        
        return inputs_defaults

    @staticmethod
    def _getIMDefaults():
        """Get a dictionary of default values for initial_model fields."""
        #TODO I'm redundantly setting things that do not make sense to have a
        #default for to None.  Helps me keep track, but maybe should just delete.
        im_defaults = {}
        im_defaults['M_tot'] = None
        im_defaults['M_He']  = None
        im_defaults['delta'] = None
        im_defaults['temp_core'] = None
        im_defaults['temp_base'] = None
        im_defaults['mixed_co_wd'] = '.false.'
        im_defaults['low_density_cutoff'] = '1.d-4'
        im_defaults['temp_fluff'] = '7.5d7'
        im_defaults['smallt'] = '1.d6'
        im_defaults['xmin'] = '0.0'

        #The initial model resolution should match Maestro's base state
        #resolution.  This is derived from inputs.  TODO Users are allowed to
        #override this, but shouldn't?
        #These should be derived:
        im_defaults['nx'] = None
        im_defaults['xmax'] = None

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

        im_defaults['radius'] = None
        im_defaults['density'] = None
        im_defaults['temperature'] = None
        im_defaults['pressure'] = None
        im_defaults['soundspeed'] = None
        im_defaults['entropy'] = None
        im_defaults['species'] = None

        return im_defaults

    def _deriveInputs(self, inputs_rec):
        """Derive inputs fields based on a partially initialized inputs ConfigRecord."""
        #TODO Derive job_name, im_file, anelastic cutoff, physical size
        pltfile_base = self._label + "_plt"
        inputs_rec.setField('plot_base_name', pltfile_base)
        miniplt_base = self._label + "_miniplt"
        inputs_rec.setField('mini_plot_base_name',miniplt_base)
        chkfile_base = self._label + "_chk"
        inputs_rec.setField('check_base_name', chkfile_base)
       
        if inputs_rec.getField('octant').lower().count('false') > 0:
            inputs_rec.setField('bc_lo', '12')
            inputs_rec.setField('bc_hi', '12')
        else:
            inputs_rec.setField('bc_lo', '13')
            inputs_rec.setField('bc_hi', '12')
    
    def _deriveIM(self, im_rec, inputs_rec):
        """Derive initial model fields based on a partially initialized inputs
        and im ConfigRecord."""
        #TODO Implement this
        pass

    #def _initInputsDict(self):
    @staticmethod
    def _getInputsFields():
        """Get a dictionary of inputs fields and their descriptions, as well as
        a mapping from file variable names to ConfigRecord field names."""
        inputs_fields = {}
        fieldmap = {}
        inputs_fields['im_file'] = 'Initial model file with data to be read into the Maestro basestate.'
        fieldmap['model_file'] = 'im_file'
        inputs_fields['job_name'] = 'Description of the simulation.'
        inputs_fields['max_levs'] = 'Number of levels the AMR will refine to.'
        inputs_fields['coarse_res'] = 'Resolution of the base (coarsest) level'
        fieldmap['n_cellx'] = 'coarse_res'
        fieldmap['n_celly'] = 'coarse_res'
        fieldmap['n_cellz'] = 'coarse_res'
        inputs_fields['anelastic_cutoff'] = 'Density cutoff below which the Maestro velocity constraint is simplified to the anelastic constraint.'
        inputs_fields['octant'] = "Boolean that sets if an octant or full star should be modeled."
        inputs_fields['dim'] = 'Dimensionality of the problem.'
        fieldmap['dm_in'] = 'dim'
        inputs_fields['physical_size'] = 'Sidelength in cm of the square domain.'
        fieldmap['prob_hi_x'] = 'physical_size'
        fieldmap['prob_hi_y'] = 'physical_size'
        fieldmap['prob_hi_z'] = 'physical_size'
        inputs_fields['plot_deltat'] = 'Time interval in s at which to save pltfiles.'
        inputs_fields['mini_plot_deltat'] = 'Time interval in s at which to save minipltfiles.'
        inputs_fields['chk_int'] = 'Timestep interval at which to save chkpoint files.'
        inputs_fields['plot_base_name'] = 'Basename for pltfiles. Pltfiles will be saved with this name plus their timestep.'
        inputs_fields['mini_plot_base_name'] = 'Basename for minipltfiles. Minipltfiles will be saved with this name plus their timestep.'
        inputs_fields['check_base_name'] = 'Basename for checkpoint files. Chkfiles will be saved with this name plus their timestep.'
        inputs_fields['bc_lo'] = 'Integer flag for the lower (x=y=z=0) boundary'
        fieldmap['bcx_lo'] = 'bc_lo'
        fieldmap['bcy_lo'] = 'bc_lo'
        fieldmap['bcz_lo'] = 'bc_lo'
        inputs_fields['bc_hi'] = 'Integer flag for the hi (x=y=z=max) boundary'
        fieldmap['bcx_hi'] = 'bc_hi'
        fieldmap['bcy_hi'] = 'bc_hi'
        fieldmap['bcz_hi'] = 'bc_hi'

        return inputs_fields, fieldmap

    @staticmethod
    def _getInputsTempText():
        """Returns the template text and leading space for an inputs file."""
        #TODO Currently, programmer should make sure fields here are the same as
        #in ConfigRecord.  Would be nice to automagically do this.
        #TODO Does this make sense as method?  Can I just define it as property
        #or some such?
        #TODO Should decide if it makes sense to have more specific format
        #specifiers.  For now, assume string.
        inputs_template = """&PROBIN
         model_file = "{im_file:s}"
         drdxfac = 5
        
         job_name = "{job_name:s}"
        
         use_alt_energy_fix = T
         ppm_trace_forces = 0
        
         hg_bottom_solver = 4
         mg_bottom_solver = 4
         max_mg_bottom_nlevels = 2
        
         max_levs = {max_levs:s}
         regrid_int = 2
        
         n_cellx = {coarse_res:s}
         n_celly = {coarse_res:s}
         n_cellz = {coarse_res:s}
        
         stop_time = 30000.
        
         max_step = 100000000
        
         init_iter = 1
         init_divu_iter = 3
         do_initial_projection = T
        
         max_grid_size_1 = 32
         max_grid_size_2 = 64
         max_grid_size_3 = 128
        
         the_sfc_threshold = 32768
        
         anelastic_cutoff = {anelastic_cutoff:s}
         base_cutoff_density = 10000.0
         buoyancy_cutoff_factor = 2.d0
         sponge_center_density = {anelastic_cutoff:s}
         sponge_start_factor = 2.0
         sponge_kappa = 10.0d0
        
         spherical_in = 1
         octant = {octant:s}
         dm_in = {dim:s}
         do_sponge = .true.
        
         prob_hi_x = {physical_size:s}
         prob_hi_y = {physical_size:s}
         prob_hi_z = {physical_size:s}
        
         plot_base_name = "{plot_base_name:s}"
         plot_int = -1
         plot_deltat = {plot_deltat:s}
        
         mini_plot_base_name = "{mini_plot_base_name:s}"
         mini_plot_int = -1
         mini_plot_deltat = {mini_plot_deltat:s}
         mini_plot_var1 = "species"
         mini_plot_var2 = "velocity"
         mini_plot_var3 = "temperature"
         mini_plot_var4 = "radial_velocity"
        
         check_base_name = "{check_base_name:s}"
         chk_int = {chk_int:s}
        
         cflfac = 0.7d0
         init_shrink = 0.1d0
         max_dt_growth = 1.1d0
         use_soundspeed_firstdt = T
         use_divu_firstdt = T
        
         bcx_lo = {bc_lo:s}
         bcx_hi = {bc_hi:s}
         bcy_lo = {bc_lo:s}
         bcy_hi = {bc_hi:s}
         bcz_lo = {bc_lo:s}
         bcz_hi = {bc_hi:s}
        
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
        lead_space = 8
        return inputs_template, lead_space

    @staticmethod
    def _getIMFields():
        """Get a dictionary of initial model fields and their descriptions."""
        im_fields = {}
        im_fields['nx'] = 'Resolution (number of cells) of the 1D model, should match Maestro base state resolution.'
        im_fields['M_tot'] = 'Mass of the WD core in M_sol.'
        im_fields['M_He']  = 'Mass of He envelope in M_sol.'
        im_fields['delta'] = 'Transition delta from core to envelope in cm.'
        im_fields['temp_core'] = 'Isothermal core temperature in K.'
        im_fields['temp_base'] = 'Temperature at the base of the He envelope in K.'
        im_fields['xmin'] = 'Spatial coordinate in cm the model starts at.'
        im_fields['xmax'] = 'Spatial coordinate in cm of the last cell, should match the sidelength of domain in octant simulation, half sidelength for full star.'
        im_fields['mixed_co_wd'] = 'Boolean that sets if core is C/O or just C.'
        im_fields['low_density_cutoff'] = 'Density floor in the initial model (NOT for the 3D Maestro domain).'
        im_fields['temp_fluff'] = 'Temperature floor, will also be temperature when below density floor.'
        im_fields['smallt'] = 'An unused parameter that used to be like temp_fluff.'

        im_fields['radius'] = 'NumPy array of initial model radius in cm.'
        im_fields['density'] = 'NumPy array of initial model density in g/cm^3.'
        im_fields['temperature'] = 'NumPy array of initial model temperature in K.'
        im_fields['pressure'] = 'NumPy array of initial model pressure in dyn/cm^2.'
        im_fields['soundspeed'] = 'NumPy array of initial model sound speed in cm/s.'
        im_fields['entropy'] = 'NumPy array of initial model specific entropy in erg/(g*K).'
        im_fields['species'] = 'NumPy 2D array of initial model species mass fractions.'
        return im_fields
 
    @staticmethod
    def _getIMTempText():
        """Returns the template text and leading space for an initial model _params file."""
        #TODO Currently, programmer should make sure fields here are the same as
        #in ConfigRecord.  Would be nice to automagically do this.
        #TODO Does this make sense as method?  Can I just define it as property
        #or some such?
        im_template = """&params
        
          nx = {nx:s}
        
          M_tot = {M_tot:s}
          M_He =  {M_He:s}
        
          delta = {delta:s}
        
          xmin = {xmin:s}
          xmax = {xmax:s}
        
          temp_core = {temp_core:s}
          temp_base = {temp_base:s}
        
          mixed_co_wd = {mixed_co_wd:s}
        
          low_density_cutoff = {low_density_cutoff:s}
          temp_fluff = {temp_fluff:s}
          smallt = 1.d6
        
        /"""
        lead_space = 8
        return im_template, lead_space

class SCOutput(SimOutput):
    """Represents the products of a sub-Chandra, such as the diagnostics files,
    reduced data from pltfiles, output from the executable, etc..."""

    def __init__(self, simdir):
        """Constructs an SCOutput object using an existing configuration in the
        given directory."""
        super(simdir)

    def _initFromDir(self, simdir):
        """Initialize this object using an existing configuration."""

        raise NotImplementedError("""A subclass of SimOutput did not implement
        this method or you're directly instantiating SimOutput.  Either way,
        NO!""")

#TODO Make driver for tests here (or make an independent test driver?)
#   + Verify input info
#   + Verify all object construction, generation works

