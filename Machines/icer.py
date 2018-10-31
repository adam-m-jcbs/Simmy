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
        self._meta_base['name'] = 'iCER HPCC Clusters'
        #TODO Add all possible hosts
        self._meta_base['hosts'] = ['gateway-00', 'dev-intel14']

        #+partitions
        #creates self._partitions['icer_gen'] and self._partitions['intel14']
        self._init_parts(self._partitions)

        #+filesys
        self._filesys_base['user_root'] = '/mnt/home/jacob308' #TODO Make arg
        self._filesys_base['scratch_root'] = None #don't use for now

        filesys_dict = {'user_root': None, 'scratch_root': None}
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

    def genBatch(self, batch_path, batch_ddict, batch_template=None):
        """
        Generate a batch script.

        Takes the full path to the generated batch file, a TemplateFile for the batch
        scripts, and a DefinedDict of the batch data what will populate
        TemplateFile.  Often, the 
        """
        if batch_template is None:
            batch_template = ICER.getBatchTemplateFile()


    @staticmethod
    def getBatchTemplateFile():
        """Return the base TemplateFile for iCER SLURM batch scripts."""
        #TODO Make TemplateFile smart about removing indendation like python
        #    does for docstrings.
        icer_slurm_template_text = """#!/bin/bash --login
###################################
#SBATCH --job-name {job_name[#SBATCH --job-name]:gs_.+}
#SBATCH --array={array_str[#SBATCH --array=]:[0-9,-]+}
#SBATCH --time={walltime[#SBATCH --time=]:[0-9]{2}.[0-9]{2}.[0-9]{2}}
#SBATCH --nodes={nodes[#SBATCH --nodes=]:[0-9]+}
#SBATCH --ntasks={ntasks[#SBATCH --ntasks=]:[0-9]+}
#SBATCH --cpus-per-task={cpus_pt[#SBATCH --cpus-per-task=]:[0-9]+} 1
#SBATCH --mem-per-cpu={mem_pc[#SBATCH --mem-per-cpu=]:[0-9]+ *[kmgtibKMGTB]*}
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user={user_email[#SBATCH --mail-user=]:.*@.*\..*}
###################################

./{exe_script[\./]:.*}
"""

        return TemplateFile(icer_slurm_template_text)

    @staticmethod
    def getBaseBatchDDict():
        """Return the base DefinedDict for batch job data."""
        #TODO This is now redundant and inconsistent with what's in super's __init__
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

