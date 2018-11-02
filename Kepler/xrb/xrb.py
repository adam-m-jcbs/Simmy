#This file contains implementations of SimulationGrid and Simulation for
#simulations of X-ray bursts using the Kepler stellar evolution code.

from simmy import SimulationGrid, Simulation, Machine

#class XRBGrid(SimulationGrid):
#    pass

class XRBSim(Simulation):
    """
    An implementation of Simulation for Kepler X-ray burst models.
    """

    def __init__(self, label, base_dir):
        """Initializes a simulation object using an existing setup in a
        directory. To generate a new simulation, use this class's static factory
        methods.
        
        Arguments:
            label      --> label for this simulation, will be used for naming
                           (e.g. directories and files)
            base_dir   --> path to the base directory this simulation is stored in
                           Simulation files will be in <base_dir>/<label>
        """
        #Initializes the general (not Kepler-/XRB-specific) components of
        #   self._simcon, self._execon, self._simout
        super().__init__(label, base_dir)


    def _initSimConfig(self, label, sim_root):
        """
        Setup a simulation configuration DefinedDict.

        This will init self._simcon
        """
        simcon_base = super()._initSimConfig(label, sim_root)
        simcon_desc = simcon_base.getDesc()
        simcon_dict = simcon_base.getData()

        #TODO: Add any xrb-specific sim config here

        return DefinedDict(simcon_desc, simcon_dict)

    def _initEXEConfig(self):
        """
        Setup a DefinedDict for data needed to specify execution instructions.

        This DefinedDict stores details used to compose shell scripts for
        executing the simulation.

        Assumption:
            self._simcon has been initialized
        """
        execon_base = super()._initEXEConfig()
        execon_desc = execon_base.getDesc()
        execon_dict = execon_base.getData()

        execon_desc['kepler_data'] = 'str, full path to the data directory Kepler expects in $KEPLER_DATA'
        execon_dict['kepler_data'] = None

        return DefinedDict(execon_desc, execon_dict)

    def _initSimOutput(self):
        """
        Setup a DefinedDict to assist with managing simulation output.

        Assumption:
            self._simcon has been initialized
        """
        simout_desc = {'output_dir': "str, name of the directory output will be stored in, defaults to 'output'"}
        simout_dict = {'output_dir': 'output'}

        return DefinedDict(execon_desc, execon_dict)

    def getEXETemplateText(batch=False):
        """
        Get the text for a TemplateFile of executable shell instructions.

        Arguments:
            batch --> bool, If True, it's assumed these instructions will be
                appended to a batch job script header (e.g. the header for a SLURM
                or TORQUE script).  Otherwise, this will be a stand-alone script
                that's intended to be executed directly
                Default: False
        """

        exe_text = """
EXE_PATH={exe_path[EXE_PATH=]:.*}
WORK_DIR={work_dir[WORK_DIR=]:.*}
export KEPLER_DATA={kepler_data[export KEPLER_DATA=]:.*}
cd $WORK_DIR

#TODO continue here


"""
        if batch:
            #TODO add batch-specific stuff
            pass

        return exe_text


