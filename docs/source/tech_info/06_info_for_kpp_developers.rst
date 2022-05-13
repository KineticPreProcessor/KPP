
.. _developer-info:

##############################
Information for KPP Developers
##############################

This chapter is meant for KPP Developers. It describes the internal
architecture of the KPP preprocessor, the basic modules and their
functionalities, and the preprocessing analysis performed on the input
files. KPP can be very easily configured to suit a broad class of users.

.. _directory-structure:

=======================
KPP directory structure
=======================

.. container:: center

   .. container::
      :name: tab:source

      .. table:: Source code files

         ==== =================================
         File Description
         \    main program
         \    generic code generation functions
         \    header file
         \    generation of C code
         \    generation of Fortran77 code
         \    generation of Fortran90 code
         \    generation of matlab code
         \    debugging output
         \    header file
         \    header file
         \    generic code generation functions
         \    flex/bison-generated file
         \    input for flex and bison
         \    input for flex
         \    input for bison
         \    evaluate parsed input
         \    evaluate parsed input
         \    flex/bison-generated file
         \    flex/bison-generated header file
         \
         ==== =================================

The KPP distribution will unfold a directory ``$KPP_HOME`` with the
following subdirectories:

-  **src/** Contains the KPP source code files, as listed in
   Table `9 <#tab:source>`__.

-  **bin/** Contains the KPP executable. The path to this directory
   needs to be added to the environment variable.

-  **util/** Contains different function templates useful for the
   simulation. Each template file has a suffix that matches the
   appropriate target language (, , , or ). KPP will run the template
   files through the substitution preprocessor. The user can define
   their own auxiliary functions by inserting them into the files.

-  **models/** Contains the description of the chemical models. Users
   can define their own models by placing the model description files in
   this directory. The KPP distribution contains several models from
   atmospheric chemistry which can be used as templates for model
   definitions.

-  **drv/** Contains driver templates for chemical simulations. Each
   driver has a suffix that matches the appropriate target language (, ,
   , or ). KPP will run the appropriate driver through the substitution
   preprocessor. The driver template provided with the distribution
   works with any example. Users can define here their own driver
   templates.

-  **int/** Contains numerical time stepping (integrator) routines. The
   command “*integrator*” will force KPP to look into this directory for
   a definition file *integrator*. This file selects the numerical
   routine (with the command) and sets the function type, the Jacobian
   sparsity type, the target language, etc. Each integrator template is
   found in a file that ends with the appropriate suffix (, , , or ).
   The selected template is processed by the substitution preprocessor.
   Users can define here their own numerical integration routines.

-  **examples/** Contains several model description examples ( files)
   which can be used as templates for building simulations with KPP.

-  **site-lisp/** Contains the file which provides a KPP mode for emacs
   with color highlighting.

.. _kpp-env-vars:

=========================
KPP environment variables
=========================

In order for KPP to find its components, it has to know the path to the
location where the KPP distribution is installed. This is achieved by
requiring the :envvar:`$KPP_HOME` environment variable to be set to the path
where KPP is installed.

The :envvar:`PATH` variable should be updated to contain the
:file:`$KPP_HOME/bin` directory.

There are also several optional environment variable that control the places
where KPP looks for module files, integrators, and drivers.  All KPP
environment variables are summarized in the subsections below.

KPP_HOME
--------
Required, stores the absolute path to the KPP distribution.

Default setting: none

KPP_MODEL
---------
Optional, specifies additional places where KPP will look for model
files before searching the default location.

Default setting: :file:`$KPP_HOME/models`.

KPP_INT
-------
Optional, specifies additional places where KPP will look for integrator files before searching the default.

Default setting: :file:`$KPP_HOME/int`.

KPP_DRV
-------

Optional specifies additional places where KPP will look for driver
files before searching the default folder.

Default setting: :file:`$KPP_HOME/drv`

.. _kpp-internal-modules:

====================
KPP internal modules
====================

Scanner and Parser
------------------

This module is responsible for reading the kinetic description files and
extracting the information necessary in the code generation phase. We
make use of the flex and bison generic tools in implementing our own
scanner and parser. Using these tools this module gathers information
from the input files and fills in the following data structures in
memory:

-  The atom list

-  The species list

-  The left hand side matrix of coefficients

-  The right hand side matrix of coefficients

-  The equation rates

-  The option list

Error checking is performed at each step in the scanner and the parser.
For each syntax error the exact line and input file, along with an
appropriate error message are produced. In most of the cases the exact
cause of the error can be identified, therefore the error messages are
very precise. Some other errors like mass balance, and equation
duplicates, are tested at the end of this phase.

Species reordering
------------------

When parsing the input files, the species list is updated as soon as a
new species is encountered in a chemical equation. Therefore the
ordering of the species is the order in which they appear in the
equation description section. This is not a useful order for subsequent
operations. The species have to be first sorted such that all variable
species and all fixed species are put together. Then if a sparsity
structure of the Jacobian is required, it might be better to reorder the
species in such a way that the factorization of the Jacobian will
preserve the sparsity. This reordering is done using a Markovitz type of
algorithm.

Expression trees computation
----------------------------

This is the core of the preprocessor. This module has to generate the
production/destruction functions the Jacobian and all the data structure
nedeed by these functions. This module has to build a language
independent structure of each function and statement in the target
source file. Instead of using an intermediate format for this as some
other compilers do, KPP generates the intermediate format for just one
statement at a time. The vast majority of the statements in the target
source file are assignments. The expression tree for each assignment is
incrementally build by scanning the coefficient matrices and the rate
constant vector. At the end these expression trees are simplified.
Similar approaches are applied to function declaration and prototypes,
data declaration and initialization.

Code generation
---------------

There are basically two modules, each dealing with the syntax
particularities of the target language. For example, the C module
includes a function that generates a valid C assignment when given an
expression tree. Similarly there are functions for data declaration,
initializations, comments, function prototypes, etc. Each of these
functions produce the code into an output buffer. A language specific
routine reads from this buffer and splits the statements into lines to
improve readability of the generated code.

Adding new KPP commands
-----------------------

To add a new KPP command, the source code has to be edited at several
locations. A short summary is presented here, using the new command as
an example:

-  Add to several files in the directory:

.. code-block:: C

   void CmdNEWCMD( char *cmd );
   -  : ``{ "NEWCMD", PRM_STATE, NEWCMD },``

   -  : ``void CmdNEWCMD( char *cmd )``

   -  :

      -  ``%token NEWCMD``

      -  ``NEWCMD PARAMETER``

      -  ``{ CmdNEWCMD( $2 ); }``

-  Maybe add a CI-test:

   -  Create a new directory

   -  Add new CI-test to and in the directory

-  Other:

   -  Explain in user manual :

      -  Add to table

      -  Add a section

      -  Add to BNF description table
