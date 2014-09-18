Introduction
============

Overview
--------

Firedrake-Fluids is a collection of finite element-based numerical
models for the study of fluid dynamical systems. It uses the `Firedrake <http://firedrakeproject.org>`_
framework to automate the solution of the governing equations written in
their weak form using the high-level, compact, near-mathematical Unified
Form Language (UFL). The complexity of writing a numerical model is
hidden through abstraction, and model developers do not need to concern
themselves with hand-writing the low-level (e.g. C or Fortran) code
required to solve the equations; this is all derived and optimised
automatically from a high-level specification. Furthermore, model
developers do not need to be experts in parallel programming to enable
their code to be performance-portable across different hardware
architectures (e.g. a cluster of multi-core CPUs, or a single GPU); the
generated code is targetted, compiled and executed automatically on a
desired architecture using the `PyOP2 <https://github.com/OP2/PyOP2>`_ 
library with which Firedrake is coupled.

Some information briefly outlining Firedrake’s automated solution
technique and the setup of Firedrake-Fluids can be found in the sections
following this one. The remaining chapters provide details on the models
available within Firedrake-Fluids, along with any auxiliary
parameterisations that the user may wish to include. In addition,
information regarding how to set up a model is also given.

Automated solution technique
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

When a model in Firedrake-Fluids is executed by the Python interpreter,
the model’s UFL (along with the computational mesh used to discretise
the domain) is first passed to the Firedrake framework. Within this
framework, the UFL is first converted to an abstract syntax tree (AST)
by a modified version of the FEniCS Form Compiler (FFC). Additionally,
the topology of the mesh is described by a `PETSc DMPlex object
<http://www.mcs.anl.gov/petsc/petsc-current/docs/manualpages/DM/DMPLEX.html>`_
to allow the efficient execution of the generated code over the whole
mesh. The DMPlex object and the AST are then passed to the PyOP2 library
which, after the AST has been optimised by the COFFEE compiler and
converted into low-level generated C code, targets and compiles the
generated code towards a specific hardware architecture and executes it
on that hardware. 

As an example, consider the UFL statement in Figure
ufl_expression_.

.. _ufl_expression:
.. figure::  images/ufl_expression.png
   :align:   center

   An example of a UFL expression.

This one single line of UFL is converted
to a kernel comprising many lines of generated C code, which perform the
evaluation of the expression, as shown in Figure c_kernel_.

.. _c_kernel:
.. figure::  images/c_kernel.png
   :align:   center

   An example of C code, generated automatically, for the purpose
   of evaluating an expression defined by a high-level, near-mathematical
   UFL statement.

Directory structure
-------------------

The directory structure of the Firedrake-Fluids codebase is as follows:

-  ``/``: The Firedrake-Fluids base directory contains general
   information in the README file, information about the license in the
   COPYING file, and a full list of authors in the AUTHORS file.

-  ``/doc``: Contains the source code and images for this
   documentation.

-  ``/firedrake_fluids``: Comprises a collection of Python files containing the
   implementation of the different models and auxiliary functionality.

-  ``/schema``: Contains a set of schema files used to define the
   different options a simulation configuration file can take (see
   Section [sect:configuring\ :sub:`as`\ imulation] for more details).

-  ``/tests``: A set of test cases to help ensure the correctness of the
   models.

Setup
-----

Dependencies
~~~~~~~~~~~~

Before running the models in Firedrake-Fluids, please ensure that all
the dependencies specified in the README file are satisfied.
Installations for Firedrake (and its dependencies) can be found `here <http://www.firedrakeproject.org/download.html>`_.
Firedrake-Fluids also
relies on the `libspud <https://launchpad.net/spud>`_ library (and the Python bindings) to retrieve
simulation-related options (e.g. the time-step size and initial
conditions) from a configuration/setup file. Following the steps below
at the command line will download and build libspud, and install the
Python bindings::

   1. ``bzr checkout lp:~spud/spud/trunk libspud``
   2. ``cd libspud``
   3. ``./configure``
   4. ``make``
   5. ``cd python``
   6. ``sudo python setup.py install``

Installation
~~~~~~~~~~~~

The Firedrake-Fluids Python module can be installed with::

   sudo python setup.py install
   
Alternatively, the ``firedrake_fluids`` directory may be added to the
``PYTHONPATH`` environment variable in order to use the module. This can be done at the command
line, e.g.::

   export PYTHONPATH=$PYTHONPATH:/home/christian/firedrake-fluids/firedrake_fluids

Following this, it is recommended that you run ``make test`` (from the
Firedrake-Fluids base directory) to ensure that the setup and models are
working correctly.

