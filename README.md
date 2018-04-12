# Simmy

*What is it, and why is it needed?*

`simmy` is a Python framework for managing simulations currently under
development.  Computational scientists must develop and deploy code on a variety
of machines (from laptops to clusters to supercomputers) which can have a
variety of configurations (e.g. different job submission systems, different
filesystem configurations, different ways for interacting with high performance
storage).  In practice, we end up writing scripts to facilitate managing all of
this.  Every computational scientist I know has written such "glue" code, often
implementing similar concepts and solving similar problems.

There are many issues with such scripts.  They can take a great deal of effort
to write and maintain, especially for scientists just starting to work with a
new code or new model/simulation.  They tend, out of expediency and an
expectation that no one else will be using them, to be specific to certain
machines and/or particular codes.  They get the job done, but what happens when
the problem needs to be run on a different machine, or if the computing facility
switches to a new job management system, or if you want to use your scripts to
manage a new problem, or if the scientist lands a sweet grant and wants to get a
graduate student quickly up to speed on carrying out a suite of simulations?
I've yet to write my own or see the scripts of another (and I've seen quite a
few) that can robustly and quickly adjust to such changes.

`simmy` is meant to be a solution to this problem.  The design goal is a
lightweight Python package providing useful abstractions that allow
computational scientists to easily and programmatically configure, build, and
execute a simulation or suite of related simulations on a variety of machines.
This includes being able to manage the data generated.

To facilitate usability, the lightweight component is key.  This framework
should be a light layer on top of the infrastructure involved in carrying out
simulations.  Only popular, widely-adopted third-party software should be
required (e.g. NumPy).  The simulation will in no way require `simmy` to run,
meaning none of the files or executables of the simulation should be aware of it
or contain any `simmy`-aware syntax or commands.  There exist many, many
"workflow" tools, commercial and otherwise, and I have never used them because
they frequently require third-party software, are expensive, cannot easily be
deployed on a range of different hardware, and/or are too complicated to use.
`simmy` should be none of these things.

Though this alone motivates `simmy`, it has another massive bit of utility: it
makes your science more reproducible.  Using `simmy` will serve to document your
exact configuration and make it readily shareable with others.  As an added
bonus, this facilitates getting people new to a code or supercomputer up to
speed quickly.

As of now, the project is focused on astrophysical codes running in a Linux
environment, such as `AMReX` codes, `Kepler`, `MESA`, `dStar`, and
`FLASH`.  Once an alpha is available with some working prototype
implementations, the code design should allow for application to other domains.
At this point, pull requests and issues will also be welcomed from the
community.
