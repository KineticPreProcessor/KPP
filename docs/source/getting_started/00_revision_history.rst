####################
KPP Revision History
####################

Only the major new features are listed here. For a detailed description
of the changes, read the file :file:`$KPP_HOME/CHANGELOG.md`.

=========
KPP-2.5.0
=========

- A new forward-Euler method integrator (:program:`feuler.f90`)

- New KPP options :code:`#MINVERSION` and :code:`#UPPERCASEF90`
  (along with corresponding continuous integration tests)

- Optional variables :code:`Aout` and :code:`Vdotout` have been added
  to subroutine Fun()

- Replaced Fortran :code:`EQUIVALENCE` statements with thread-safe pointer
  assignments (Fortran90 only)

- New documentation in Sphinx/ReadTheDocs format

- Updates to allow KPP to be built on MacOS X systems

=========
KPP-2.4.0
=========

- New integrators: :file:`beuler.f90`, :file:`rosenbrock_mz.f90`,
  :file:`rosenbrock_posdef.f90,  :file:`rosenbrock_posdef_h211b_qssa.f90`
  
- Several memory sizes (:code:`MAX_EQN`, ...) increased to allow large
  chemical mechanisms 
  
- new Makefile target: :code:`list`
  
- LaTeX User Manual added

- Now use :code:`ICNTRL(15)` to decide whether or not to toggle calling the
  :code:`Update_SUN`, :code:`Update_RCONST`, and :code:`Update_PHOTO`
  routines from within the integrator

=========
KPP 2.2.3
=========

-  A new function called :code:`k_3rd_iupac` is available, calculating
   third-order rate coefficients using the formula used by IUPAC
   :raw-latex:`\citep{1610}`.

-  While previous versions of KPP were using :program:`yacc` (yet another
   compiler compiler), the current version has been modified to be
   compatible with the parser generator :program:`bison`, which is the
   successor of :program:`yacc`.

-  The new Runge-Kutta integrators :code:`kpp_sdirk4`, :code`runge_kutta`, and
   :code:`sdirk` were added.

-  The new KPP command :code:`#DECLARE` was added (see
   Sect. `3.2.1 <#sec:command-declare>`__).

=======
KPP 2.1
=======

This user manual describes recently added features of KPP as well as
those which have been available for a longer period. Here we give an
overview about the recent changes:

-  Fortran90 output has been available since the preliminary version
   “1.1-f90-alpha12” provided in :raw-latex:`\citet{1666}`.

-  Matlab is a new target language (see Sect. `4.4 <#sec:matlab>`__).

-  The set of integrators has been extended with a general Rosenbrock
   integrator, and the corresponding tangent linear and adjoint methods.

-  The KPP-generated Fortran90 code has a different file structure than
   the C or Fortran77 output (see Sects. `4.2 <#sec:c>`__ and
   `4.3 <f77>`__).

-  An automatically generated Makefile facilitates the compilation of
   the KPP-generated code (see Sect. `4.1.18 <#sec:output-makefile>`__).

-  Equation tags provide a convenient way to refer to specific chemical
   reactions (see Sect. `4.1.5 <#sec:output-monitor>`__).

-  The dummy index allows to test if a certain species occurs in the
   current chemistry mechanism. (see
   Sect. `3.2.4 <#sec:command-dummyindex>`__).

-  Lines starting with :code:`//` are comment lines.
