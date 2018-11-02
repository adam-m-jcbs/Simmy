# This code implements the Machine classes from the simmy framework to represent
# the iCER clusters at Michigan State University.

# Author:        Adam Jacobs
# Creation date: October 30, 2018
from simmy import Machine, TemplateFile, DefinedDict

class ICER(Machine):
    """Represents the iCER cluster."""

    def __init__(self):
        #This will serve to initialize many base variables we can then add to as
        #needed for iCER.  They are:
        #   self._meta_base, self._partitions, self._partition_base,
        #   self._filesys_base
        super().__init__()

        #Now go through and initialize the Machine that is currently a
        #template/null definition.
        #+meta
        self._meta = DefinedDict(self._meta_base._data, self._meta_base._desc)
        self._meta['name'] = 'iCER HPCC Clusters'
        #TODO Add all possible hosts
        self._meta['hosts'] = ['gateway-00', 'dev-intel14']
        self._meta['has_batch'] = True

        #+partitions
        #creates self._partitions['icer_gen'] and self._partitions['intel14']
        self._init_parts(self._partitions)

        #+filesys
        self._filesys = DefinedDict(self._filesys_base._data, self._filesys_base._desc)
        self._filesys['user_root'] = '/mnt/home/jacob308' #TODO Make arg
        self._filesys['scratch_root'] = None #don't use for now

        #TODO add defensive superclass method for checking init is successful
        #(all non-optional data is initialized)

    def _init_parts(self, partitions):
        """
        Initialize partition DefinedDicts for ICER instance.
        
        For now, initializes some psuedo-partitions that are actually
        collections of partitions.  I only declare the minimum available
        resources available to all nodes in the partition. This inits

        partitions['icer_gen'] and partitions['intel14']
        """
        #TODO iCER's a bit different.  It's a cluster of clusters, but all
        #   clusters are byte-compatible and are usually treated as a single machine
        #   by SLURM.  For now, defining an icer_generic to represent all of it, and
        #   one example of a particular cluster: intel-14.
        #+generic representation of the whole machine, which forces me to only
        #   declare the minimum resources available per node
        icer_gen = DefinedDict(self._partition_base._data, self._partition_base._desc)
        icer_gen['arch'] = "A mix of Intel Xeon E-series, with 10-18 cores"
        icer_gen['node_count'] = 595
        icer_gen['hw_cores_pn'] = 20 #Many 20, many 28, and a few even up to 96
        icer_gen['logical_core_fac'] = 1
        icer_gen['gpus_pn'] = 0 #many GPUS, but not on all nodes, so can't declare
        icer_gen['mem_domains_pn'] = 1 #I'm not actually sure, but for now I don't need to do NUMA optimization
        icer_gen['mem_per_domain'] = "64 GB" #This is the minimum guaranteed, many nodes have 256 GB and some nodes have up to 6 TB
        partitions['icer_gen'] = icer_gen

        #+intel-14 cluster
        #TODO: intel-14 is actually composed of 4 sub-clusters.  In the
        #   definition below, I defined resources as the minimum available to any
        #   sub-cluster.  This is hiding resources, which is OK for now but need to
        #   update after prototyping
        #   See this for description of iCER hardware: 
        #   https://wiki.hpcc.msu.edu/display/hpccdocs/Description+of+the+Processing+Hardware
        i14_part = DefinedDict(self._partition_base._data, self._partition_base._desc)
        i14_part['arch'] = "Intel Xeon E5-2670v2" #2.5 Ghz 10-core, 2 per node
        i14_part['node_count'] = 220
        i14_part['hw_cores_pn'] = 20
        i14_part['logical_core_fac'] = 1
        i14_part['gpus_pn'] = 0 #intel14 actually has some, but only on 40 nodes
        i14_part['mem_domains_pn'] = 1 #I'm not actually sure, but for now I don't need to do NUMA optimization
        i14_part['mem_per_domain'] = "64 GB" #This is the minimum guaranteed, some nodes have up to 256 GB
        partitions['intel14'] = i14_part 

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
                               script's pragmas.  These are the mostly
                               simulation/problem-specific shell statements that
                               will be executed, and that's why we separate them
                               from the machine-specific batch script template.
            batch_template --> TemplateFile, optional, the TemplateFile to be
                               used to generate the batch script.  By default will be generated
                               with static method getBatchTemplateFile().
        """
        from math import floor
        #NOTE: As recommended in Machine, this method will first derive
        #   iCER-specific data and then pass this along to Machine's genBatch

        #Get copies of the batch DefinedDict descriptions and data.
        icer_desc = batch_ddict.getDesc()
        icer_data = batch_ddict.getData()

        #Derive ntasks
        icer_desc['ntasks'] = "int, Total number of tasks requested"
        tasks_pn = batch_ddict['tasks_pn']
        nodes = batch_ddict['nodes']
        icer_data['ntasks'] = tasks_pn * nodes
        icer_batch_ddict = DefinedDict(icer_data, icer_desc)

        #If not given by user, derive max mem per task
        #if icer_data['mem_pt'] == '':
        #    mem_pt = self.partitions['icer_gen']icer_data['']

        #Derive memory per cpu
        icer_desc['mem_pc'] = "int, Memory per cpu in MB"
        mem_pt_mb = Machine.memStr2MB(batch_ddict['mem_pt'])
        icer_data['mem_pc'] = floor(
                float(mem_pt_mb) / 
                float(batch_ddict['cores_pt'])
            )

        #Create expanded DDict
        icer_batch_ddict = DefinedDict(icer_data, icer_desc)

        #Pass along to super, now that we've translated the data
        super().genBatch(batch_path, icer_batch_ddict, exe_text, batch_template)

    @staticmethod
    def getBatchTemplateFile():
        """Return the base TemplateFile for iCER SLURM batch scripts."""
        #TODO Make TemplateFile smart about removing indendation like python
        #   does for docstrings.
        #TODO Figure out a way to make this responsive to different
        #   configurations, e.g. a template with or without array jobs
        icer_slurm_template_text = """#!/bin/bash --login
###################################
#SBATCH --job-name {job_name[#SBATCH --job-name]:gs_.+}
#SBATCH --array={array_str[#SBATCH --array=]:[0-9,-]+}
#SBATCH --time={walltime[#SBATCH --time=]:[0-9]{2}.[0-9]{2}.[0-9]{2}}
#SBATCH --ntasks={ntasks[#SBATCH --ntasks=]:[0-9]+}
#SBATCH --cpus-per-task={cores_pt[#SBATCH --cpus-per-task=]:[0-9]+} 1
#SBATCH --mem-per-cpu={mem_pc[#SBATCH --mem-per-cpu=]:[0-9]+ *[kmgtibKMGTB]*}
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user={user_email[#SBATCH --mail-user=]:.*@.*\..*}
###################################

#Prepare the environment
source ${HOME}/PyVE/bin/activate
module load GNU/7.3.0-2.30
N=$SLURM_ARRAY_TASK_ID
"""

        return TemplateFile(icer_slurm_template_text)

    #@staticmethod
    #def getBaseBatchDDict():
    #    """Return the base DefinedDict for batch job data."""
    #    #TODO This is now redundant and inconsistent with what's in super's __init__
    #    batch_desc = {'job_name': 'str, Name to be used for the job, no spaces!',
    #                  'walltime': 'str, The requested walltime, in form [DD:]HH:MM:SS',
    #                  'array_str': 'str, optional, If you want an array job, give a valid array str here, e.g. "1-24", "1,4,8-9". "" if you do not want an array job (default)',
    #                  'nodes': 'int, Number of nodes requested',
    #                  'tasks_pn': 'int, Number of tasks (roughly, processes or MPI tasks) per node, defaults to 1',
    #                  'cores_pt': 'int, Number of cores (cpus) per task',
    #                  'mem_pt': 'str, memory requested per task in MB, defaults to most memory guaranteed to be available (system-dependent)',
    #                  'user_email': 'str, optional, email to send any alerts to',
    #                  'exe_script': 'str, filename for the script to be executed in the batch submission, put all simulation-specific execution instructions here.'}
    #    batch_dict = {'job_name': None,
    #                  'walltime': None,
    #                  'array_str': "",
    #                  'nodes': None,
    #                  'tasks_pn': 1,
    #                  'mem_pt': None,
    #                  'user_email': "",
    #                  'exe_script': None}
    #    
    #    return DefinedDict(batch_dict, batch_desc)

