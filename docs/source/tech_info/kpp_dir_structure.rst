.. _directory-structure:

#######################
KPP directory structure
#######################

The KPP distribution will unfold a directory :envvar:`$KPP_HOME` with the
following subdirectories:

====
src/
====

Contains the KPP source code files:

.. _table-kpp-dirs:

.. list-table:: KPP source code files
   :align: center
   :widths: 30 70
   :header-rows: 1

   * - File
     - Description
   * - :file:`kpp.c`
     - Main program
   * - :file:`code.c`
     - generic code generation functions
   * - :file:`code.h`
     - Header file
   * - :file:`code_c.c`
     - Generation of C code
   * - :file:`code_f90.c`
     - Generation of F90 code
   * - :file:`code_matlab.c`
     - Generation of Matlab code
   * - :file:`debug.c`
     - Debugging output
   * - :file:`gdata.h`
     - Header file
   * - :file:`gdef.h`
     - Header file
   * - :file:`gen.c`
     - Generic code generation functions
   * - :file:`lex.yy.c`
     - Flex generated file
   * - :file:`scan.h`
     - Input for Flex and Bison
   * - :file:`scan.l`
     - Input for Flex
   * - :file:`scan.y`
     - Input for Bison
   * - :file:`scanner.c`
     - Evaluate parsed input
   * - :file:`scanutil.c`
     - Evaluate parsed input
   * - :file:`y.tab.c`
     - Bison generated file
   * - :file:`y.tab.h`
     - Bison generated header file

====       
bin/
====

Contains the KPP executable. This directory should be added to the
:envvar:`PATH` environment variable.

=====	
util/
=====

Contains different function templates useful for the simulation. Each
template file has a suffix that matches the appropriate target
language (Fortran90, C, or Matlab). KPP will run the template files
through the substitution preprocessor (cf.
:ref:`list-of-symbols-replaced`). The user can define their own
auxiliary functions by inserting them into the files.

=======
models/
=======

Contains the description of the chemical models. Users
can define their own models by placing the model description files in
this directory. The KPP distribution contains several models from
atmospheric chemistry which can be used as templates for model
definitions.

====
drv/
====

Contains driver templates for chemical simulations. Each driver has a
suffix that matches the appropriate target language (Fortran90, C, or
Matlab). KPP will run the appropriate driver through the substitution
preprocessor (cf. :ref:`list-of-symbols-replaced`). Users can also
define their own driver templates here.

====
int/
====

Contains numerical solvers (integrators). The :command:`#INTEGRATOR`
command will force KPP to look into this directory for a definition
file with suffix :code:`.def`. This file selects the numerical solver
etc. Each integrator template is found in a file that ends with the
appropriate suffix (:code:`.f90`, :code:`.c`, or :code:`.m`). The
selected template is processed by the substitution preprocessor (cf.
:ref:`list-of-symbols-replaced`). Users can define their own
numerical integration routines in the :code:`user_contributed`
subdirectory.

=========
examples/
=========

Contains several model description examples (:file:`.kpp` files)
which can be used as templates for building simulations with KPP.

==========
site-lisp/
==========
Contains the file :file:`kpp.el` which provides a KPP mode for emacs
with color highlighting.

=========
ci-tests/
=========

Contains directories defining several :ref:`ci-tests`.

==============
.ci-pipelines/
==============

Hidden directory containing a YAML file with settings for automatically
running the continuous integration tests on `Azure DevOps Pipelines
<https://azure.microsoft.com/en-us/services/devops/pipelines/>`_

Also contains bash scripts (ending in :file:`.sh`) for running the
continuous integration tests either automatically in Azure Dev
Pipelines, or manually from the command line.  For more
information, please see :ref:`ci-tests`.

=================
.github/workflows
=================

Contains configuration files for `GitHub Actions
<https://github.com/KineticPreProcessor/KPP/actions>`_ that will run
automatically when commits are pushed or when pull requests are opened.
