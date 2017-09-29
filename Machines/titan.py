# This code implements the Machine classes from the simmy framework to represent
# the Titan supercomputer at Oak Ridge National Lab.  The code is based in part
# on code I previously wrote to carry out the models described in
# http://adsabs.harvard.edu/abs/2016ApJ...827...84J

# Author:        Adam Jacobs
# Creation date: September 28, 2017
from simmy import RunConfig, TemplateFile

class TitanConfig(RunConfig):
    """Represents the configuration and files needed to execute a simulation on
    a particular machine.  
    
    An instance of this class is stored by the Machine class.
    """

    def __init__(self, simdir, config_dict=None):
        """Constructs a TitanConfig object for a simulation found in
        simdir.
        
        If a simulation already exists in simdir, it will be used to create this
        class.  If not, you need to provide config_dict.  The run's
        label will be the name of the base directory of simdir.

        Below are the config_dict keys expected along with example defaults.
        Those with * are typically unique to a simulation and should have been
        initialized.  Others will be given reasonable defaults.
            allocation_label --> Name of the allocation to be charged
                                 Default: 'ast106'
            job_label*   --> The name this run will have in the queue
                             Default: '10044-090-210-4lev-full-512'
            walltime     --> The amount of wall time requested for the script to run
                             Default: '04:00:00'
            nodes        --> Number of nodes to request
                             Default: '256'
            threads_per_node --> The number of threads per node
                                 Default: '16'
            process_script --> The name of a script to run in parallel with the
                               job that will process output (e.g. archive to HPSS)
                               Default: 'process.titan'
            exe_label    --> The name of the executable to run
                             Default: 'main.Linux.Cray.mpi.omp.exe'
            inputs_file* --> The name of the inputs file to be passed to the executable
                             Default: 'inputs3d.10044-090-210-4lev-full-512'
        """
        from os.path import join
        #Will initialize self._label and self._simdir, will call subclass method
        #to initialize self._config_dict if config_dict not given.
        self._runscript = join(simdir, 'model', 'titan.pbs')
        super().__init__(simdir, config_dict)

    def _initFromDir(self):
        """Initialize this object using an existing configuration found in
        self._simdir."""
        #Read in runscript and parse file for dictionary values
        #ASSUMPTION: self._runscript has been defined
        #TODO if I don't end up adding more logic here, just replace with
        #initConfigDict
        self._initConfigDictFromFile()

    def _initConfigDictFromFile(self):
        """Initialize self._config_dict from an existing runscript.
        
        config_dict keywords are
            allocation_label --> Name of the allocation to be charged
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
        from re import match
        from argparse import ArgumentParser, REMAINDER
        #TODO Much of this makes assumptions that makes sense for sub-Chandra,
        #but I should modify to be more generic for Titan.  Should make defaults
        #more generic, for example, and probably for some things defaults
        #shouldn't be allowed -- user must provide them.

        #Define a parser for aprun. TODO add more possible args used
        aprun_parser = ArgumentParser(prog='aprun',
            description='Command for launching parallel jobs on compute nodes')
        aprun_parser.add_argument('-n')
        aprun_parser.add_argument('-N')
        aprun_parser.add_argument('-d')
        aprun_parser.add_argument('-ss', action='store_true')
        aprun_parser.add_argument('exe')
        aprun_parser.add_argument('exe_args', nargs=REMAINDER)

        #Define regex strings
        allocation_re = ' *#PBS *-A.*'
        job_re = ' *#PBS *-N.*'
        resources_re = ' *#PBS *-l.*'
        
        #TODO Make sure I don't define config_dict before now
        self._config_dict = {}
       
        #TODO For now, assume process.titan, but might be nice to infer from runscript
        self._config_dict['process_script'] = 'process.titan'
        #Process each line of the runscript
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
                    self._config_dict['threads_per_node'] = line.partition('=')[2].strip()
                if line.lstrip().startswith('aprun'):
                    #TODO Does parser break if something like '-N1' is done?
                    #   I've known such arg syntax to be permitted, so might need to
                    #   watch for this.
                    args_string = line.lstrip()[5:] #Everything but `aprun`
                    aprun_args = aprun_parser.parse_args(args_string.split())
                    self._config_dict['exe_label'] = aprun_args.exe[2:] #Strip leading `./`
                    if int(self._config_dict['nodes']) != int(aprun_args.n):
                        #TODO Maybe right my own exception?
                        raise UserWarning("Warning: Nodes requested doesn't match nodes used by aprun.")
                    self._config_dict['inputs_file'] = aprun_args.exe_args[0] 

    def _initFromDict(self):
        """Initialize this object using the configuration dictionary found in
        self._config_dict.  The runscript associated with this dictionary will
        be created and saved.
        
        Note that the config_dict will be fully initialized with any missing
        values and any derived values."""
        #ASSUMPTION: self._config_dict has been defined with the minimum number
        #            of keys needed to specify a configuration.  This can be
        #            done for example by passing a dict into the constructor.
        #Set reasonable defaults for any missing values
        self._initConfigDict()

        #Generate and save the runscript file
        self._genRunscript()
        
        #Finally, fully initialize the runscript by reading in the newly
        #created file.  This makes sure initFromDict is consistent with
        #initFromFile.  TODO Can I avoid this?
        self._initConfigDictFromFile()

    def _genRunscript(self):
        """Generate and save the runscript file.  self._config_dict
        keys will be used.

        config_dict keys used.  Those with * are typically unique to a
        simulation and should have been initialized in the dictionary:
            allocation_label --> Name of the allocation to be charged
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
        #Define the base template string
        runscript_template = """#!/bin/ksh
        #PBS -A {allocation_label:s}
        #PBS -N {job_label:s}
        #PBS -j oe
        #PBS -q batch
        #PBS -l walltime={walltime},nodes={nodes}
        
        # this script runs with {threads_per_node} threads, 1 MPI tasks/node,
        # and {nodes} nodes on titan
        
        export PSC_OMP_AFFINITY=FALSE
        export OMP_NUM_THREADS={threads_per_node}
        
        # These are needed for very large MAESTRO runs, such as
        # inputs_3d.2304.5dr.eq.dx_3levels with 6 threads
        #export MPICH_PTL_OTHER_EVENTS=16384
        #export MPICH_UNEX_BUFFER_SIZE=100000000
        #export MPICH_MAX_SHORT_MSG_SIZE=10000
        
        cd $PBS_O_WORKDIR
        
        # run the compression script to tar up the plot and checkpoint files
        # as they are created.
        ./{process_script} &
        PID=$!
        trap 'kill -s TERM $PID' EXIT TERM HUP XCPU KILL
        
        # find the latest restart file -- first look for one with 6 digits then fall
        # back to 5
        restartFile=$(find . -maxdepth 1 -type d -name "*chk??????" -print | sort | tail -1)
        
        # the Header is the last thing written -- check if it's there, otherwise,
        # fall back to the second-to-last check file written
        if [ ! -f ${{restartFile}}/Header ]; then
          # how many *chk?????? files are there? if only one, then skip
          nl=$(find . -maxdepth 1 -type d -name "*chk??????" -print | sort | wc -l)
          if [ $nl -gt 1 ]; then
        	  restartFile=$(find . -maxdepth 1 -type d -name "*chk??????" -print | sort | tail -2 | head -1)    
          else
        	  restartFile=""
          fi
        fi
        
        # if the above checks failed, then there are no valid 6-digit chk files, so
        # check the 5-digit ones
        if [ "${{restartFile}}" = "" ]; then
          restartFile=$(find . -maxdepth 1 -type d -name "*chk?????" -print | sort | tail -1)
        
          # make sure the Header was written, otherwise, check the second-to-last
          # file
          if [ ! -f ${{restartFile}}/Header ]; then
            # how many *chk????? files are there? if only one, then skip
            nl=$(find . -maxdepth 1 -type d -name "*chk?????" -print | sort | wc -l)
            if [ $nl -gt 1 ]; then
        	    restartFile=$(find . -maxdepth 1 -type d -name "*chk?????" -print | sort | tail -2 | head -1)    
            else
        	    restartFile=""
            fi
          fi
        fi
        
        # cut out the numerical part of the *chkXXXXX file, here we use the 'k' in 
        # 'chk' as the delimiter
        restartNum=`echo ${{restartFile}} | cut -d'k' -f2`  
        
        # restartString will be empty if no chk files are found -- i.e. new run
        if [ "${{restartNum}}" = "" ]; then
            restartString=""
        else
            restartString="--restart ${{restartNum}}"
            echo "Restarting with: " ${{restartString}}
        fi
        
        # Titan has 18688 physical nodes, each of which has 16 cores and 2 NUMA nodes
        # 
        # -n  is the total number of MPI tasks
        # -S  is the number of MPI tasks per NUMA node 
        # -N  is the number of MPI tasks per node
        # -d  is the number of OpenMP threads per MPI task (must match OMP_NUM_THREADS)
        # -ss forces MPI tasks to only allocate memory in their local NUMA node.
        #   This can boost performance by preventing costly remote memory I/O, though 
        #   it also restricts the amount of memory available to MPI tasks.
        
        aprun -n {nodes} -N 1 -d $OMP_NUM_THREADS -ss ./{exe_label} {inputs_file} ${{restartString}}
        
        rm -f process.pid
        """
        #TODO Need to derive N vs S, calculate MPI/node based on nodes and
        #threads per node

        #Now create and save the file
        runscript_file = TemplateFile(self._config_dict, runscript_template, 8)
        runscript_file.saveFile(self._runscript)

    def _initConfigDict(self):
        """Initialize a dictionary of key properties describing the runscript
        submitted to the supercomputer queue.
       
        Initialize self._config_dict with default values for any keys missing in
        the current self._config_dict.
        """
        #TODO Generalize this to run on either supercomputer/cluster or local machine
        #     Actually, this should be in Machine.  This handles the TODO, and
        #     makes it possible to configure running on various machines.
        config_defaults = {}
        config_defaults['allocation_label'] = 'ast106'
        config_defaults['job_label'] = self._label
        config_defaults['walltime'] = '04:00:00'
        config_defaults['nodes'] = '256'
        config_defaults['threads_per_node'] = '16'
        #TODO Make process_script optional
        config_defaults['process_script'] = 'process.titan'
        config_defaults['exe_label'] = 'main.Linux.Cray.mpi.omp.exe'
        config_defaults['inputs_file'] = 'inputs3d.10044-090-210-4lev-full-512'

        for key in config_defaults:
            if key not in self._config_dict:
                self._config_dict[key] = config_defaults[key]
