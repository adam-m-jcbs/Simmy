# This code implements the Machine classes from the simmy framework to represent
# the Titan supercomputer at Oak Ridge National Lab.  The code is based in part
# on code I previously wrote to carry out the models described in
# http://adsabs.harvard.edu/abs/2016ApJ...827...84J

# Author:        Adam Jacobs
# Creation date: September 28, 2017

class TitanConfig(RunConfig):
    """Represents the configuration and files needed to execute a simulation on
    a particular machine.  
    
    An instance of this class is stored by the Machine
    class.
    """

    def __init__(self, simdir, config_dict=None):
         """Constructs a TitanConfig object for a simulation found in
        simdir.
        
        If a simulation already exists in simdir, it will be used to create this
        class.  If not, you need to provide config_dict.  The run's
        label will be the name of the base directory of simdir.
        """
        from os.path import join
        #Will initialize self._label and self._simdir, will call subclass method
        #to initialize self._config_dict if config_dict not given.
        self._runscript = join(simdir, 'titan.pbs')
        super(simdir, config_dict)

    def _initFromDir(self):
        """Initialize this object using an existing configuration found in
        self._simdir."""
        #Read in runscript
        
        #Parse file for dictionary values
        
        #return dictionary
        pass

    def _initConfigDict(self):
        """Initialize self._config_dict from an existing runscript."""
        from re import match

        allocation_re = ' *#PBS *-A.*'
        job_re = ' *#PBS *-N.*'
        resources_re = ' *#PBS *-l.*'
        #TODO Make sure I don't define config_dict before now
        self._config_dict = {}
       
        #TODO For now, assume process.titan, but might be nice to infer from runscript
        self._config_dict['process_script'] = 'process.titan'
        with open(self._runscript, 'r') as f:
            for line in f:
                if match(allocation_re, line):
                    self._config_dict['allocation_label'] = line.partition('-A')[2].strip()
                if match(job_re, line):
                    self._config_dict['job_label'] = line.partition('-N')[2].strip()
                if match(resources_re, line):
                    resources_args = line.partition('-l')[2].strip()
                    for arg in resources_args.split(','):
                        if arg.count('walltime') > 0:
                            self._config_dict['walltime'] = arg.partition('=')[2].strip()
                        if arg.count('nodes') > 0:
                            self._config_dict['nodes'] = arg.partition('=')[2].strip()
                if line.count('OMP_NUM_THREADS') > 0:
                    #TODO Restart here
                    self._config_dict['threads_per_node'] = line
                if line.lstrip().startswith('aprun'):

                    
        run_dict['threads_per_node'] = 
        run_dict['process_script'] = 
        run_dict['exe_label'] = 




    def _initFromDict(self):
        """Initialize this object using the configuration dictionary found in
        self._config_dict.  The runscript associated with this dictionary will
        be created and  saved.
        
        Note that the config_dict will be fully initialized with any missing key
        values and any derived values."""
        #Define dictionary for runscript

        #Generate and save the runscript file
        
        #Return the fully initialized dictionary

    def _genRunscript(self, savepath):
        """Generate the runscript file and save to savepath.  self._config_dict
        keys will be used, with missing keys given default values.

        Arguments:
            savepath  --> Full path to the location where this file should be saved

        config_dict keys used.  Those with * are typically unique to a
        simulation and should have been initialized in the dictionary:
            job_label* --> The name this run will have in the queue
            walltime   --> The amount of wall time requested for the script to run
            nodes      --> Number of nodes to request
            threads_per_node --> The number of threads per node
            process_script --> The name of a script to run in parallel with the
                               job that will process output (e.g. archive to HPSS)
            exe_label  --> The name of the executable to run
            inputs_file --> The name of the inputs file to be passed to the
                            executable
        """



    def _initRunDict(self):
        """Initialize a dictionary of key properties describing the runscript
        submitted to the supercomputer queue."""
        #TODO Generalize this to run on either supercomputer/cluster or local machine
        #     Actually, this should be in Machine.  This handles the TODO, and
        #     makes it possible to configure running on various machines.
        
        run_dict = self._config_dicts['run_dict']
        run_dict['allocation_label'] = 
        run_dict['job_label'] = 
        run_dict['walltime'] = 
        run_dict['nodes'] = 
        run_dict['threads_per_node'] = 
        run_dict['process_script'] = 
        run_dict['exe_label'] = 





