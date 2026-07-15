.. _ci-tests:

############################
Continuous integration tests
############################

KPP contains several continuous integration (aka C-I) tests. Each C-I
test calls KPP to generate source code for a given
:ref:`chemical mechanism <model-cmd>`, :ref:`integrator
<integrator-cmd>`, and :ref:`target language <language-cmd>`, and
then runs a short "box model" simulation with the generated code. C-I
tests help to ensure that new features and updates added to KPP will
not break any existing functionality.

C-I tests will run automatically as a `GitHub Action
<https://github.com/KineticPreProcessor/KPP/actions>`_ when commits
are pushed to the `KPP Github repository
<https://github.com/KineticPreProcessor/KPP>`_, or when a new pull
requests are opened.  You may also run the integration tests
:ref:`locally on your own computer <running-ci-tests-locally>`.

.. _list-of-ci-tests:

====================================
List of continuous integration tests
====================================

.. list-table:: Continuous integration tests
   :header-rows: 1
   :align: center

   * - C-I test
     - Language
     - Model
     - Integrator
     - What it tests
   * - :code:`C_rk`
     - C
     - small_strato
     - runge_kutta
     - :ref:`Runge-Kutta integrator <rk-methods-3stage>`
   * - :code:`C_rosadj`
     - C
     - small_strato
     - rosenbrock_adj
     - :ref:`Rosenbrock discrete adjoint <rosenbrock-adjoint>`
   * - :code:`C_sd`
     - C
     - small_strato
     - sdirk
     - :ref:`SDIRK integrator <rk-methods-sdirk>`
   * - :code:`C_sdadj`
     - C
     - small_strato
     - sdirk_adj
     - :ref:`SDIRK discrete adjoint <rk-methods-sdirk>`
   * - :code:`C_small_strato`
     - C
     - small_strato
     - rosenbrock
     - :ref:`Rosenbrock integrator <rosenbrock-methods>`
   * - :code:`F90_feuler`
     - Fortran90
     - carbon
     - feuler
     - :ref:`Forward Euler integrator <other-methods-feuler>`
   * - :code:`F90_graph`
     - Fortran90
     - small_strato
     - rosenbrock
     - :ref:`#GRAPH command <graph-cmd>`
   * - :code:`F90_lsode`
     - Fortran90
     - small_strato
     - lsode
     - :ref:`LSODE integrator <back-diff-lsode>`
   * - :code:`F90_mcm`
     - Fortran90
     - mcm
     - rosenbrock
     - Master Chemical Mechanism
   * - :code:`F90_mcm_h211b`
     - Fortran90
     - mcm
     - rosenbrock_h211b_qssa
     - :ref:`Rosenbrock int. w/ H211b time stepping <rosenbrock-h211b-qssa>`
   * - :code:`F90_radau`
     - Fortran90
     - saprc99
     - radau5
     - :ref:`RADAU5 integrator <rk-methods-radau5>`
   * - :code:`F90_rk`
     - Fortran90
     - small_strato
     - runge_kutta
     - :ref:`Runge-Kutta integrator <rk-methods-3stage>`
   * - :code:`F90_rkadj`
     - Fortran90
     - small_strato
     - runge_kutta_adj
     - :ref:`Runge-Kutta discrete adjoint <rk-methods-adj>`
   * - :code:`F90_rktlm`
     - Fortran90
     - small_strato
     - runge_kutta_tlm
     - :ref:`Runge-Kutta tangent linear model <rk-methods-tlm>`
   * - :code:`F90_ros`
     - Fortran90
     - small_strato
     - rosenbrock
     - :ref:`Rosenbrock integrator <rosenbrock-methods>`
   * - :code:`F90_rosadj`
     - Fortran90
     - small_strato
     - rosenbrock_adj
     - :ref:`Rosenbrock discrete adjoint <rosenbrock-adjoint>`
   * - :code:`F90_ros_autoreduce`
     - Fortran90
     - saprc99
     - rosenbrock_autoreduce
     - :ref:`Rosenbrock integrator w/ auto-reduction <rosenbrock-autoreduce>`
   * - :code:`F90_rosenbrock`
     - Fortran90
     - saprc99
     - rosenbrock
     - :ref:`Rosenbrock integrator <rosenbrock-methods>`
   * - :code:`F90_ros_h211b`
     - Fortran90
     - saprc99
     - rosenbrock_h211b_qssa
     - :ref:`Rosenbrock int. w/ H211b time stepping <rosenbrock-h211b-qssa>`
   * - :code:`F90_ros_passivespc`
     - Fortran90
     - small_strato
     - rosenbrock
     - Excluding passive species
   * - :code:`F90_ros_split`
     - Fortran90
     - small_strato
     - rosenbrock
     - :ref:`#FUNCTION SPLIT <function-cmd>`
   * - :code:`F90_rostlm`
     - Fortran90
     - small_strato
     - rosenbrock_tlm
     - :ref:`Rosenbrock tangent linear model <rosenbrock-tlm>`
   * - :code:`F90_ros_upcase`
     - Fortran90
     - saprc99
     - rosenbrock
     - :ref:`#UPPERCASEF90 ON <uppercasef90-cmd>`
   * - :code:`F90_saprc_2006`
     - Fortran90
     - saprcnov
     - rosenbrock
     - :ref:`Rosenbrock integrator <rosenbrock-methods>`
   * - :code:`F90_sd4`
     - Fortran90
     - small_strato
     - sdirk4
     - :ref:`Rosenbrock integrator <rosenbrock-methods>`
   * - :code:`F90_sd`
     - Fortran90
     - small_strato
     - sdirk
     - :ref:`SDIRK integrator <rk-methods-sdirk>`
   * - :code:`F90_sdadj`
     - Fortran90
     - small_strato
     - sdirk_adj
     - :ref:`SDIRK discrete adjoint <rk-methods-sdirk>`
   * - :code:`F90_sdtlm`
     - Fortran90
     - small_strato
     - sdirk_tlm
     - :ref:`SDIRK tangent linear model <rk-methods-sdirk>`
   * - :code:`F90_seulex`
     - Fortran90
     - saprcnov
     - seulex
     - :ref:`SEULEX integrator <rk-methods-seulex>`
   * - :code:`F90_small_strato`
     - Fortran90
     - small_strato
     - rosenbrock
     - :ref:`Rosenbrock integrator <rosenbrock-methods>`
   * - :code:`X_minver`
     - Fortran90
     - small_strato
     - runge_kutta
     - :ref:`#MINVERSION command <minversion-cmd>`

Notes about C-I tests:

#. :file:`F90_ros_split` uses :command:`#FUNCTION SPLIT`.
#. :file:`F90_ros_upcase` uses :command:`#UPPERCASEF90 ON`.
#. :file:`F90_small_strato` is the example from
   :ref:`running-kpp-with-an-example-mechanism`.
#. :file:`X_minver` tests if the :ref:`minversion-cmd` command works
   properly.

Each continuous integration test is contained in a subdirectory of
:file:`$KPP_HOME/ci-tests`.  In each subdirectory is a KPP definition
file (ending in :file:`.kpp`).

.. _running-ci-tests:

=======================================================
Running continuous integration tests as a GitHub Action
=======================================================

The files needed to run the C-I tests are described below.

.. _running-ci-tests-action:

run-ci-tests.yml
----------------

**Path**: :file:`$KPP_HOME/.github/workflows/run-ci-tests.yml`

**Description:** Configuration file with commands to download KPP, load
libraries, and run the C-I tests as a GitHub Action.

C-I tests will run automatically when a commit is pushed to any branch
at `https://github.com/KineticPreProcessor/KPP
<https://github.com/KineticPreProcessor/KPP>`_, or when a new pull
request is opened there.  This is the recommended setting, but you can
restrict this so that only pushes or pull requests to certain branches
will trigger the C-I tests.

.. _running-ci-tests-testing:

ci-testing-script.sh
--------------------

**Path:** :file:`$KPP_HOME/.ci-pipelines/ci-testing-script.sh`

**Description:**  Runs the KPP C-I tests as a GitHub Action, or on a
local computer system.

.. _running-ci-tests-cleanup:

ci-cleanup-script.sh
--------------------

**Path:** :file:`$KPP_HOME/.ci-pipelines/ci-cleanup-script.sh`

**Description:** Removes compiler-generated files (e.g. :file:`*.o`,
:file:`.mod` , and  :file:`.exe`) from C-I test folders.

.. _running-ci-tests-defs:

ci-common-defs.sh
-----------------

**Path:** :file:`$KPP_HOME/.ci-pipelines/ci-common-defs.sh`

**Description** Contains common variable and function definitions needed by
:ref:`running-ci-tests-testing` and :ref:`running-ci-tests-cleanup`.

.. _running-ci-tests-locally:

============================================
Running continuous integration tests locally
============================================

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
