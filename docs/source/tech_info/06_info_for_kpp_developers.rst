.. _developer-info:

##############################
Information for KPP developers
##############################

This chapter describes the internal architecture of the KPP
preprocessor, the basic modules and their functionalities, and the
preprocessing analysis performed on the input files. KPP can be very
easily configured to suit a broad class of users.

.. _directory-structure:

=======================
KPP directory structure
=======================

The KPP distribution will unfold a directory :envvar:`$KPP_HOME` with the
following subdirectories:

.. option:: src/

   Contains the KPP source code files:

.. _table-kpp-dirs:

.. table:: KPP source code files
   :align: center

   +-----------------------+-------------------------------------+
   | File                  | Description                         |
   +=======================+=====================================+
   | :file:`kpp.c`         | Main program                        |
   +-----------------------+-------------------------------------+
   | :file:`code.c`        | generic code generation functions   |
   +-----------------------+-------------------------------------+
   | :file:`code.h`        | Header file                         |
   +-----------------------+-------------------------------------+
   | :file:`code_c.c`      | Generation of C code                |
   +-----------------------+-------------------------------------+
   | :file:`code_f90.c`    | Generation of F90 code              |
   +-----------------------+-------------------------------------+
   | :file:`code_matlab.c` | Generation of Matlab code           |
   +-----------------------+-------------------------------------+
   | :file:`debug.c`       | Debugging output                    |
   +-----------------------+-------------------------------------+
   | :file:`gdata.h`       | Header file                         |
   +-----------------------+-------------------------------------+
   | :file:`gdef.h`        | Header file                         |
   +-----------------------+-------------------------------------+
   | :file:`gen.c`         | Generic code generation functions   |
   +-----------------------+-------------------------------------+
   | :file:`lex.yy.c`      | Flex generated file                 |
   +-----------------------+-------------------------------------+
   | :file:`scan.h`        | Input for Flex and Bison            |
   +-----------------------+-------------------------------------+
   | :file:`scan.l`        | Input for Flex                      |
   +-----------------------+-------------------------------------+
   | :file:`scan.y`        | Input for Bison                     |
   +-----------------------+-------------------------------------+
   | :file:`scanner.c`     | Evaluate parsed input               |
   +-----------------------+-------------------------------------+
   | :file:`scanutil.c`    | Evaluate parsed input               |
   +-----------------------+-------------------------------------+
   | :file:`y.tab.c`       | Bison generated file                |
   +-----------------------+-------------------------------------+
   | :file:`y.tab.h`       | Bison generated header file         |
   +-----------------------+-------------------------------------+

.. option:: bin/

   Contains the KPP executable. This directory should be added to the
   :envvar:`PATH` environment variable.

.. option:: util/

   Contains different function templates useful for the simulation. Each
   template file has a suffix that matches the appropriate target
   language (Fortran90, C, or Matlab). KPP will run the template files
   through the substitution preprocessor (cf.
   :ref:`list-of-symbols-replaced`). The user can define their own
   auxiliary functions by inserting them into the files.

.. option:: models/

   Contains the description of the chemical models. Users
   can define their own models by placing the model description files in
   this directory. The KPP distribution contains several models from
   atmospheric chemistry which can be used as templates for model
   definitions.

.. option:: drv/

   Contains driver templates for chemical simulations. Each driver has a
   suffix that matches the appropriate target language (Fortran90, C, or
   Matlab). KPP will run the appropriate driver through the substitution
   preprocessor (cf. :ref:`list-of-symbols-replaced`). Users can also
   define their own driver templates here.

.. option:: int/

   Contains numerical solvers (integrators). The :command:`#INTEGRATOR`
   command will force KPP to look into this directory for a definition
   file with suffix :code:`.def`. This file selects the numerical solver
   etc. Each integrator template is found in a file that ends with the
   appropriate suffix (:code:`.f90`, :code:`.c`, or :code:`.m`). The
   selected template is processed by the substitution preprocessor (cf.
   :ref:`list-of-symbols-replaced`). Users can define their own
   numerical integration routines in the :code:`user_contributed`
   subdirectory.

.. option:: examples/

   Contains several model description examples (:file:`.kpp` files)
   which can be used as templates for building simulations with KPP.

.. option:: site-lisp/

   Contains the file :file:`kpp.el` which provides a KPP mode for emacs
   with color highlighting.

.. option:: ci-tests/

   Contains directories defining several :ref:`ci-tests`.

.. option:: .ci-pipelines/

   Hidden directory containing a YAML file with settings for automatically
   running the continuous integration tests on `Azure DevOps Pipelines
   <https://azure.microsoft.com/en-us/services/devops/pipelines/>`_

   Also contains bash scripts (ending in :file:`.sh`) for running the
   continuous integration tests either automatically in Azure Dev
   Pipelines, or manually from the command line.  For more
   information, please see :ref:`ci-tests`.

.. _kpp-env-vars:

=========================
KPP environment variables
=========================

In order for KPP to find its components, it has to know the path to the
location where the KPP distribution is installed. This is achieved by
setting the :envvar:`$KPP_HOME` environment variable to the path where
KPP is installed.

The :file:`$KPP_HOME/bin` directory. should be added to the
:envvar:`PATH` variable.

There are also several optional environment variable that control the places
where KPP looks for module files, integrators, and drivers:

.. option:: KPP_HOME

   Required, stores the absolute path to the KPP distribution.

   Default setting: none.

.. option:: KPP_FLEX_LIB_DIR

   Optional. Use this to specify the path to the :ref:`flex library
   file <flex-dep>` (:file:`libfl.so` or :file:`libfl.a`) that are
   needed to :ref:`build the KPP executable <build-kpp-exec>`. The KPP
   build sequence will use the path contained in
   :envvar:`KPP_FLEX_LIB_DIR` if the flex library file cannot be found
   in  :file:`/usr/lib`, :file:`/usr/lib64`, and similar standard
   library paths.

.. option:: KPP_MODEL

   Optional, specifies additional places where KPP will look for model
   files before searching the default location.

   Default setting: :file:`$KPP_HOME/models`.

.. option:: KPP_INT

   Optional, specifies additional places where KPP will look for
   integrator files before searching the default.

   Default setting: :file:`$KPP_HOME/int`.

.. option:: KPP_DRV

   Optional specifies additional places where KPP will look for driver
   files before searching the default directory.

   Default setting: :file:`$KPP_HOME/drv`.

.. _kpp-internal-modules:

====================
KPP internal modules
====================

.. _scanner-parser:

Scanner and parser
------------------

This module is responsible for reading the kinetic description files and
extracting the information necessary in the code generation phase. We
make use of the flex and bison generic tools in implementing our own
scanner and parser. Using these tools, this module gathers information
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
appropriate error message are produced. Some other errors like mass
balance, and equation duplicates, are tested at the end of this phase.

.. _species-reordering:

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
preserve the sparsity. This reordering is done using a Markovitz type
algorithm.

.. _expression-trees:

Expression trees computation
----------------------------

This is the core of the preprocessor. This module generates the
production/destruction functions, the Jacobian and all the data
structure nedeed by these functions. It builds a language-independent
structure of each function and statement in the target source file.
Instead of using an intermediate format for this as some other compilers
do, KPP generates the intermediate format for just one statement at a
time. The vast majority of the statements in the target source file are
assignments. The expression tree for each assignment is incrementally
built by scanning the coefficient matrices and the rate constant vector.
At the end, these expression trees are simplified. Similar approaches are
applied to function declaration and prototypes, data declaration and
initialization.

.. _code-generation:

Code generation
---------------

There are basically two modules, each dealing with the syntax
particularities of the target language. For example, the C module
includes a function that generates a valid C assignment when given an
expression tree. Similarly there are functions for data declaration,
initializations, comments, function prototypes, etc. Each of these
functions produce the code into an output buffer. A language-specific
routine reads from this buffer and splits the statements into lines to
improve readability of the generated code.

.. _adding-new-commands:

=======================
Adding new KPP commands
=======================

To add a new KPP command, the source code has to be edited at several
locations. A short summary is presented here, using :code:`NEWCMD` as an
example:

* Add the new command to several files in the :file:`src/` directory:

  - :file:`scan.h`: add :code:`void CmdNEWCMD( char *cmd );`

  - :file:`scan.l`: add :code:`{ "NEWCMD", PRM_STATE, NEWCMD },`

  - :file:`scanner.c`: add :code:`void CmdNEWCMD( char *cmd )`

  - :file:`scan.y`:

    - Add :code:`%token NEWCMD`

    - Add :code:`| NEWCMD PARAMETER`

    - Add :code:`{ CmdNEWCMD( $2 ); }`

* Add :ref:`ci-tests`:

  - Create a new directory :file:`ci-tests/ros_newcmd/ros_newcmd.kpp`

  - Add new :ref:`ci-tests` to the :file:`ci-tests` directory and
    update the scripts in the :file:`.ci-pipelines` directory.

* Other:

  - Explain in user manual :file:`docs/source/*/*.rst`:

    - Add to Table :ref:`table-cmd-defaults`

    - Add a new subsection to :ref:`kpp-commands`

    - Add to the Table :ref:`bnf-description`

  - Add to :file:`site-lisp/kpp.el`

.. _ci-tests:

============================
Continuous integration tests
============================

KPP contains several continuous integration (aka C-I) tests. Each C-I
test calls KPP to generate source code for a given
:ref:`chemical mechanism <model-cmd>`, :ref:`integrator
<integrator-cmd>`, and :ref:`target language <language-cmd>`, and
then runs a short "box model" simulation with the generated code. C-I
tests help to ensure that new features and updates added to KPP will
not break any existing functionality.

The continuous integration tests will run automatically on `Azure
DevOps Pipelines
<https://azure.microsoft.com/en-us/services/devops/pipelines/>`_ each time a
commit is pushed to the `KPP Github repository
<https://github.com/KineticPreProcessor/KPP>`_.  You can also run the
integration tests :ref:`locally on your own computer
<running-ci-tests-locally>`.

.. _list-of-ci-tests:

List of continuous integration tests
------------------------------------

.. _table-ci-tests:

.. table:: Continuous integration tests
   :align: center

   +------------------------+-----------+-----------------+-----------------------+
   | C-I test               | Language  | Model           | Integrator            |
   +========================+===========+=================+=======================+
   | ``C_rk``               | C         | small_strato    | runge_kutta           |
   +------------------------+-----------+-----------------+-----------------------+
   | ``C_rosadj``           | C         | small_strato    | rosenbrock_adj        |
   +------------------------+-----------+-----------------+-----------------------+
   | ``C_sd``               | C         | small_strato    | sdirk                 |
   +------------------------+-----------+-----------------+-----------------------+
   | ``C_sdadj``            | C         | small_strato    | sdirk_adj             |
   +------------------------+-----------+-----------------+-----------------------+
   | ``C_small_strato``     | C         | small_strato    | rosenbrock            |
   +------------------------+-----------+-----------------+-----------------------+
   | ``F90_lsode``          | Fortran90 | small_strato    | lsode                 |
   +------------------------+-----------+-----------------+-----------------------+
   | ``F90_radau``          | Fortran90 | saprc99         | radau5                |
   +------------------------+-----------+-----------------+-----------------------+
   | ``F90_rk``             | Fortran90 | small_strato    | runge_kutta           |
   +------------------------+-----------+-----------------+-----------------------+
   | ``F90_rktlm``          | Fortran90 | small_strato    | runge_kutta_tlm       |
   +------------------------+-----------+-----------------+-----------------------+
   | ``F90_ros``            | Fortran90 | small_strato    | rosenbrock            |
   +------------------------+-----------+-----------------+-----------------------+
   | ``F90_ros_autoreduce`` | Fortran90 | saprc99         | rosenbrock_autoreduce |
   +------------------------+-----------+-----------------+-----------------------+
   | ``F90_ros_split``      | Fortran90 | small_strato    | rosenbrock            |
   +------------------------+-----------+-----------------+-----------------------+
   | ``F90_ros_upcase``     | Fortran90 | saprc99         | rosenbrock            |
   +------------------------+-----------+-----------------+-----------------------+
   | ``F90_rosadj``         | Fortran90 | small_strato    | rosenbrock_adj        |
   +------------------------+-----------+-----------------+-----------------------+
   | ``F90_rosenbrock``     | Fortran90 | saprc99         | rosenbrock            |
   +------------------------+-----------+-----------------+-----------------------+
   | ``F90_rostlm``         | Fortran90 | small_strato    | rosenbrock_tlm        |
   +------------------------+-----------+-----------------+-----------------------+
   | ``F90_saprc_2006``     | Fortran90 | saprcnov        | rosenbrock            |
   +------------------------+-----------+-----------------+-----------------------+
   | ``F90_sd``             | Fortran90 | small_strato    | sdirk                 |
   +------------------------+-----------+-----------------+-----------------------+
   | ``F90_sdadj``          | Fortran90 | small_strato    | sdirk_adj             |
   +------------------------+-----------+-----------------+-----------------------+
   | ``F90_seulex``         | Fortran90 | saprcnov        | seulex                |
   +------------------------+-----------+-----------------+-----------------------+
   | ``F90_small_strato``   | Fortran90 | small_strato    | rosenbrock            |
   +------------------------+-----------+-----------------+-----------------------+
   | ``X_minver``           | Fortran90 | small_strato    | runge_kutta           |
   +------------------------+-----------+-----------------+-----------------------+

Notes about C-I tests:

#. :file:`F90_ros_split` also uses :command:`#FUNCTION SPLIT`.
#. :file:`F90_ros_upcase` also uses :command:`#UPPERCASEF90 ON`.
#. :file:`F90_small_strato` is the example from
   :ref:`running-kpp-with-an-example-mechanism`.
#. :file:`X_minver` tests if the :ref:`minversion-cmd` command works
   properly.

Each continuous integration test is contained in a subdirectory of
:file:`$KPP_HOME/ci-tests`.  In each subdirectory is a KPP definition
file (ending in :file:`.kpp`).

.. _running-ci-tests-on-azure:

Running continuous integration tests on Azure DevOps Pipelines
--------------------------------------------------------------

The files needed to run the C-I tests are located in the
:file:`$KPP_HOME/.ci-pipelines` directory:

.. _table-ci-pipelines:

.. table:: Files needed to execute C-I tests
   :align: center

   +-------------------------------+------------------------------------------+
   | File                          | Description                              |
   +===============================+==========================================+
   | :file:`Dockerfile`            | File containing specifications for the   |
   |                               | Docker container that will be used to    |
   |                               | run C-I tests on Azure DevOps Pipelines. |
   |                               | Also contains commands needed to run     |
   |                               | the C-I scripts in the Docker container. |
   +-------------------------------+------------------------------------------+
   | :file:`build_testing.yml`     | Contains options for triggering C-I      |
   |                               | tests on Azure DevOps Pipelines.         |
   +-------------------------------+------------------------------------------+
   | :file:`ci-testing-script.sh`  | Driver script for running C-I tests.     |
   |                               | Can be used on Azure DevOps Pipelines    |
   |                               | or on a local computer.                  |
   +-------------------------------+------------------------------------------+
   | :file:`ci-cleanup-script.sh`  | Script to remove compiler-generated      |
   |                               | files (e.g. ``*.o``, ``.mod``, and       |
   |                               | ``.exe``) from C-I test folders.         |
   +-------------------------------+------------------------------------------+
   | :file:`ci-common-defs.sh`     | Script with common variable and function |
   |                               | definitions needed by                    |
   |                               | :file:`ci-testing-script.sh` and         |
   |                               | :file:`ci-cleanup-script.sh`.            |
   +-------------------------------+------------------------------------------+

The :file:`Dockerfile` contains the software environment for `Azure
DevOps Pipelines
<https://azure.microsoft.com/en-us/services/devops/pipelines/>`_.  You
should not have to update this file.

File :file:`build_testing.yml` defines the runtime options for Azure
DevOps Pipelines.  The following settings determine which branches
will trigger C-I tests:

.. code-block:: yaml

   # Run a C-I test when a push to any branch is made.
   trigger:
     branches:
       include:
          - '*'
   pr:
     branches:
       include:
         - '*'

Currently this is set to trigger the C-I tests when a commit or pull
request is made to any branch of
`https://github.com/KineticPreProcessor/KPP
<https://github.com/KineticPreProcessor/KPP>`_. This is the recommended
setting, but you can restrict this so that only pushes or pull requests
to certain branches will trigger the C-I tests.

The script :file:`ci-testing-script.sh` executes all of the C-I tests
whenever a push or a pull request is made to the selected branches in
the KPP Github repository.

.. _running-ci-tests-locally:

Running continuous integration tests locally
--------------------------------------------

To run the C-I tests on a local computer system, use this command:

.. code-block:: console

   $ $KPP_HOME/.ci-pipelines/ci-testing-script.sh | tee ci-tests.log

This will run all C-I tests on your own computer system and pipe the
results to a log file. This will easily allow you to check if the
results of the C-I tests are identical to C-I tests that were run on a
prior commit or pull request.

To remove the files generated by the continuous integration tests, use
this command:

.. code-block :: console

   $ $KPP_HOME/.ci-pipelines/ci-cleanup-script.sh

If you add new C-I tests, be sure to add the name of the new tests to
the variable :code:`GENERAL_TESTS` in :file:`ci-common-defs.sh`.
