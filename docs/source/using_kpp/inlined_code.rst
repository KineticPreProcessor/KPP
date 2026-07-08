.. |br| raw:: html

   <br />

.. _inlined-code:

############
Inlined code
############

In order to offer maximum flexibility, KPP allows the user to include
pieces of code in the kinetic description file. Inlined code begins on a
new line with :command:`#INLINE` and the *inline_type*. Next, one or
more lines of code follow, written in the target language (Fortran90, C,
or Matlab) as specified by the *inline_type*. The inlined code ends with
:command:`#ENDINLINE`. The code is inserted into the KPP output at a
position which is also determined by *inline_type* as shown in
:ref:`table-inl-type`. If two inline commands with the same inline type
are declared, then the contents of the second is appended to the first
one.

.. _list-of-inlined-types:

=====================
List of inlined types
=====================

In this manual, we show the inline types for Fortran90. The inline
types for the other languages are produced by replacing :code:`F90`
by :code:`C`, or :code:`matlab`, respectively.

.. _table-inl-type:

.. list-table:: KPP inlined types
   :align: center
   :header-rows: 1

   * - Inline_type
     - File
     - Placement
     - Usage
   * - **F90_DATA**
     - :ref:`Monitor`
     - specification section
     - (obsolete)
   * - **F90_GLOBAL**
     - :ref:`Global`
     - specification section
     - global variables
   * - **F90_INIT**
     - :ref:`Initialize`
     - subroutine
     - integration parameters
   * - **F90_RATES**
     - :ref:`Rates`
     - executable section
     - rate law functions
   * - **F90_RCONST**
     - :ref:`Rates`
     - subroutine
     - rate coefficient definitions
   * - **F90_RCONST_USE**
     - :ref:`Rates`
     - subroutine
     - rate coefficient definitions
   * - **F90_UTIL**
     - :ref:`Util`
     - executable section
     - utility functions

.. _f90-data:

========
F90_DATA
========

This inline type was introduced in a previous version of KPP to
initialize variables. It is now obsolete but kept for compatibility. For
Fortran90, :command:`F90_GLOBAL` should be used instead.

.. _f90-global:

==========
F90_GLOBAL
==========

This inline type can be used to declare global variables, e.g. for a
special rate coefficient:

.. code-block:: fortran

   #INLINE F90_GLOBAL
     REAL(dp) :: k_DMS_OH
   #ENDINLINE

Inlining code can be useful to introduce additional state variables
(such as temperature, humidity, etc.) for use by KPP routines, such as
for calculating rate coefficients.

If a large number of state variables needs to be held in inline code, or
require intermediate computation that may be repeated for many rate
coefficients, a derived type object should be used for efficiency, e.g.:

.. code-block:: fortran

   #INLINE F90_GLOBAL
     TYPE, PUBLIC :: ObjGlobal_t
        ! ... add variable fields to this type ...
     END TYPE ObjGlobal_t
     TYPE(ObjGlobal_t), TARGET, PUBLIC :: ObjGlobal
   #ENDINLINE

This global variable :code:`ObjGlobal` can then be used globally in KPP.

Another way to avoid cluttering up the KPP input file is to
:code:`#include` a header file with global variables:

.. code-block:: fortran

   #INLINE F90_GLOBAL
   ! Inline common variables into KPP_ROOT_Global.f90
   #include "commonIncludeVars.f90"
   #ENDINLINE

In future versions of KPP, the global state will be reorganized into
derived type objects as well.

.. _inline-type-f90-init:

========
F90_INIT
========

This inline type can be used to define initial values before the start of the
integration, e.g.:

.. code-block:: fortran

   #INLINE F90_INIT
     TSTART = (12.*3600.)
     TEND = TSTART + (3.*24.*3600.)
     DT = 0.25*3600.
     TEMP = 270.
   #ENDINLINE

.. _f90-rates:

=========
F90_RATES
=========

This inline type can be used to add new subroutines to calculate rate
coefficients, e.g.:

.. code-block:: fortran

   #INLINE F90_RATES
     REAL FUNCTION k_SIV_H2O2(k_298,tdep,cHp,temp)
       ! special rate function for S(IV) + H2O2
       REAL, INTENT(IN) :: k_298, tdep, cHp, temp
       k_SIV_H2O2 = k_298 &
         * EXP(tdep*(1./temp-3.3540E-3)) &
         * cHp / (cHp+0.1)
     END FUNCTION k_SIV_H2O2
   #ENDINLINE

.. _f90-rconst:

==========
F90_RCONST
==========

This inline type can be used to define time-dependent values of rate
coefficients.  You may inline variables directly, e.g.:

.. code-block:: fortran

   #INLINE F90_RCONST
     k_DMS_OH = 1.E-9*EXP(5820./temp)*C(ind_O2)/ &
       (1.E30+5.*EXP(6280./temp)*C(ind_O2))
   #ENDINLINE

The inlined code will be placed directly into the subroutines
:code:`UPDATE_RCONST` and :code:`UPDATE_PHOTO` in the :ref:`Rates` file.

.. _f90-rconst-use:

==============
F90_RCONST_USE
==============

Similar to :ref:`f90-rconst`, but allows you to inline Fortran-90
:code:`USE` statements referencing modules where rate coefficients are
computed, such as:

.. code-block:: fortran

   #INLINE F90_RCONST_USE
     USE MyRateFunctionModule
   #ENDINLINE

The inlined code will be placed directly into the subroutines
:code:`UPDATE_RCONST` and :code:`UPDATE_PHOTO` in the :ref:`Rates`
file.  :code:`USE` statements will be placed before Fortran variable
definitions and executable statements, as is required by the
Fortran-90 language standard.

.. _f90-util:

========
F90_UTIL
========

This inline type can be used to define utility subroutines.
