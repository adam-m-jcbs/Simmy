# This code implements the Machine classes from the simmy framework to represent
# a generic Linux workstation or laptop.

# Author:        Adam Jacobs
# Creation date: October 13, 2017
from simmy import Machine, RunConfig, TemplateFile

class LinuxMachine(Machine):
    """Represents a generic linux setup.

    This is intended to work for most workstations or laptops running recent
    Linux distros.
    """
