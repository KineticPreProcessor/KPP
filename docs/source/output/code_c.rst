.. |br| raw:: html

   <br />

.. _C-code:

##########
The C code
##########

.. important::

   Some run-time options for C-language integrators (specified in
   the :ref:`ICNTRL and RCNTRL arrays <icntrl-rcntrl>`) do not exactly
   correspond to the Fortran90 run-time options.  We will standardize
   run-time integrator options across all target languages in a future
   KPP release.

The driver file :file:`ROOT.c` contains the main (driver) program and
numerical integrator functions, as well as declarations and
initializations of global variables.

The generated C code includes three header files which are
:code:`#include`-d in other files as appropriate.

#. :ref:`Global parameters <table-par>` are :code:`#include`-d
   in the header file :file:`ROOT_Parameters.h`. |br|
   |br|

#. :ref:`Global variables <table-glob>` are  extern-declared in
   :file:`ROOT_Global.h` and declared in the driver file
   :file:`ROOT.c`. |br|
   |br|

#. Extern declarations of sparse data structures for the
   :ref:`Jacobian <table-jac>`, :ref:`Hessian <table-hess>` and
   :ref:`stoichiometric matrix  <table-sto>`, and the :ref:`Jacobian
   of reaction products <table-jvrp>` are in the header file
   :file:`ROOT_Sparse.h`.  The actual declarations of each
   data structure is done in the corresponding files.

The code for the :ref:`ODE function <Function>` is in
:file:`ROOT_Function.c`.

The :ref:`ODE Jacobian and sparse multiplications
<Jacobian-and-JacobianSP>` is in :file:`ROOT_Jacobian.c`, and the
declaration and initialization of the Jacobian sparse data structures
is in the file :file:`ROOT_JacobianSP.c`.

Similarly, the :ref:`Hessian function and associated sparse
multiplications <Hessian-and-HessianSP>`) are in
:file:`ROOT_Hessian.c`, and the declaration and initialization of
Hessian sparse data structures are in :file:`ROOT_HessianSP.c`.

Functions for :ref:`reactant products and its Jacobian, and
derivatives with respect to rate coefficients
<Stoichiom-and-StoichiomSP>` are in the file :file:`ROOT_Stoichiom.c`.

The declaration and initialization of the :ref:`stoichiometric matrix
and the associated sparse data structures <table-sto>`) is done in
:file:`ROOT_StoichiomSP.c`.

:ref:`Sparse linear algebra routines <LinearAlgebra>` are in the file
:file:`ROOT_LinearAlgebra.c`.  The code to update the rate constants
and user defined code for rate laws is in :file:`ROOT_Rates.c`.

:ref:`Various utility and input/output functions <Util>` are in
:file:`ROOT_Util.c` and :file:`ROOT_Monitor.c`.

Finally, :ref:`mex gateway routines <mex-code>` that allow the C
implementation of the ODE function, Jacobian, and Hessian to be called
directly from Matlab are also generated in the files
:file:`ROOT_mex_Fun.c`, :file:`ROOT_mex_Jac_SP.c`, and
:file:`ROOT_mex_Hessian.c`.
