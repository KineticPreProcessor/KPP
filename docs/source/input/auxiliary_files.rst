.. _auxiliary-files-and-the-substitution-preprocessor:

#################################################
Auxiliary files and the substitution preprocessor
#################################################

The `auxiliary files <auxiliary-files-for-fortran-90_>`_ in the
:file:`$KPP_HOME/util` subdirectory are templates for integrators,
drivers, and utilities. They are inserted into the KPP output after
being run through the substitution preprocessor. This preprocessor
replaces `several placeholder symbols <list-of-symbols-replaced_>`_ in
the template files with their particular values in the model at hand.
Usually, only :command:`KPP_ROOT` and :command:`KPP_REAL` are needed
because the other values can also be obtained via the variables listed
in :ref:`table-inl-type`.

:command:`KPP_REAL` is replaced by the appropriate single or double
precision declaration  type. Depending on the target language KPP will
select the correct declaration type. For example if one needs to
declare an array BIG of size 1000, a declaration like the following
must be used:

.. code-block:: fortran

   KPP_REAL :: BIG(1000)

When used with the command :command:`#DOUBLE ON`, the above line will be
automatically translated into:

.. code-block:: fortran

   REAL(kind=dp) :: BIG(1000)

and when used with the command :command:`#DOUBLE OFF`, the same line will
become:

.. code-block:: fortran

   REAL(kind=sp) :: BIG(1000)

in the resulting Fortran90 output file.

:command:`KPP_ROOT` is replaced by the root file name of the main kinetic
description file.  In our example where we are processing
:file:`small_strato.kpp`, a line in an auxiliary Fortran90 file like

.. code-block:: fortran

   USE KPP_ROOT_Monitor

will be translated into

.. code-block:: fortran

   USE small_strato_Monitor

in the generated Fortran90 output file.

.. _auxiliary-files-for-fortran-90:

=====================================
List of auxiliary files for Fortran90
=====================================

.. _table-aux-files:

.. list-table:: Auxiliary files for Fortran90
   :align: center
   :header-rows: 1

   * - File
     - Contents
   * - :code:`dFun_dRcoeff.f90`
     - Derivatives with respect to reaction rates
   * - :code:`dJac_dRcoeff.f90`
     - Derivatives with respect to reaction rates
   * - :code:`Makefile_f90` and :code:`Makefile_upper_F90`
     - Makefiles to build Fortran-90 code
   * - :code:`Mex_Fun.f90` and :code:`Mex_Jac_SP.f90`
     - Mex files.
   * - :code:`Mex_Hessian.f90`
     - Mex files.
   * - :code:`sutil.f90`
     - Sparse utility functions.
   * - :code:`tag2num.f90`
     - Function related to equation tags.
   * - :code:`UpdateSun.f90`
     - Function related to solar zenith angle.
   * - :code:`UserRateLaws.f90` and :code:`UserRateLawsInterfaces.f90`
     - User-defined rate-law functions.
   * - :code:`util.f90`
     - Input/output utilities.

.. _list-of-symbols-replaced:

=========================================================
List of symbols replaced by the substitution preprocessor
=========================================================

.. _table-sym-repl:

.. list-table:: Symbols and their replacements
   :align: center
   :header-rows: 1

   * - Symbol
     - Replacement
     - Example
   * - **KPP_ROOT**
     - The :literal:`ROOT` name
     - :literal:`small_strato`
   * - **KPP_REAL**
     - The real data type
     - :code:`REAL(kind=dp)`
   * - **KPP_NSPEC**
     - Number of species
     - 7
   * - **KPP_NVAR**
     - Number of variable species
     - 5
   * - **KPP_NFIX**
     - Number of fixed species
     - 2
   * - **KPP_NREACT**
     - Number of chemical reactopms
     - 10
   * - **KPP_NONZERO**
     - Number of Jacobian nonzero elements
     - 18
   * - **KPP_LU_NONZERO**
     - Number of Jacobian nonzero elements, with LU fill0in
     - 19
   * - **KPP_LU_NHESS**
     - Number of Hessian nonzero elements
     - 19
   * - **KPP_FUN_OR_FUN_SPLIT**
     - Name of the function to be called
     - ``FUN(Y,FIX,RCONST,Ydot)``
