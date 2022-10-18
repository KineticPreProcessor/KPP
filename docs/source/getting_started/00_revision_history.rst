.. _kpp-revision-history:

####################
KPP revision history
####################

Only the major new features are listed here. For a detailed description
of the changes, read `CHANGELOG.md
<https://github.com/KineticPreProcessor/KPP/blob/main/CHANGELOG.md>`_.

.. _kpp300:

=========
KPP 3.0.0
=========

.. attention::

   When you are upgrading from an older KPP version to KPP 3.0.0 or later
   versions, a few minor changes in your code may be necessary:

   - The :file:`atoms` file is now called :file:`atoms.kpp`. Thus, you have
     to change :code:`#INCLUDE atoms` to :code:`#INCLUDE atoms.kpp` in your
     KPP input file.

   - The utility functions :code:`ARR`, :code:`ARR2`, :code:`k_3rd` and
     :code:`k_arr` have been replaced by the new set of the consistent
     functions :code:`ARR_abc`, :code:`ARR_ab`, :code:`ARR_ac`,
     :code:`k3rd_jpl`, :code:`k3rd_jpl_activation`, and
     :code:`k3rd_iupac`. We recommend to upgrade to the new functions,
     which all use the temperature from the :code:`temp` variable in
     :file:`ROOT_Global.f90`. Alternatively, it is possible to copy the
     old functions into a separate file and make them available via
     :ref:`f90-rconst`.

   - If you have been using :code:`ICNTRL(5)` for maximal order in the
     :code:`lsode` integrator, you now have to use :code:`ICNTRL(10)`
     instead. The index 5 in the :code:`ICNTRL` array is now used
     consistently for the maximum number of Newton iterations in all
     integrators.


- Updated the search for the :program:`flex` library in
  :file:`src/Makefile.defs`.  The build process will look for the
  :program:`flex:` library file (either :file:`libfl.so` or
  :file:`libfl.a`  file in several standard locations first.  If not
  found, the build process will look in the path specfiied by
  environment variable :envvar:`KPP_FLEX_LIB_DIR`.

- Added content to ReadTheDocs pages and fixed several formatting issues.

- Fixed various minor issues in generating C-language code.

- Fixed various minor issues in generating Matlab-language code.

- C-I tests folders have been renamed for clarity.  Also refactored
  the scripts used to submit C-I tests.  Updated the Dockerfile to
  always request Ubuntu 20.04 and an AMD64 platform, so that the same
  libraries will always be used when running C-I tests on Azure
  DevOps.

- Fortran type :code:`DOUBLE_COMPLEX` is now replaced by
  :code:`COMPLEX(kind=dp)`.

- Fixed incorrect license metadata in :file:`.zenodo.json`, which is
  used to auto-generate a DOI with each KPP release on Github.

- Added extra :code:`free()` statements in :file:`src/gen.c` to avoid
  memory leaks.

- :code:`Fun()` no longer uses Vdotout since it can be retrieved from
  :code:`Vdot`.

- Fixed a bug in :file:`int/feuler.f90`, where the wrong argument was
  being passed to routine :code:`Fun`.

.. _kpp260:

=========
KPP 2.6.0
=========

- Added the **rosenbrock_autoreduce** integrator :cite:t:`Lin_et_al._2022`.

.. _kpp250:

=========
KPP 2.5.0
=========

- Merged updates from the GEOS-Chem development stream (versions
  :ref:`kpp224gc`, :ref:`kpp225gc`, :ref:`kpp230gc`, :ref:`kpp231gc`,
  :ref:`kpp232gc` ) into the mainline KPP development
  stream.  Previously hardwired code has been removed and replaced
  with code selectable via KPP commands.

- Added a new forward-Euler method integrator (:file:`feuler.f90`).

- Added KPP commands :command:`#MINVERSION` and :command:`#UPPERCASEF90`
  (along with corresponding continuous integration tests).

- Added optional variables :code:`Aout` and :code:`Vdotout`
  to subroutine Fun().

- Replaced Fortran :code:`EQUIVALENCE` statements with thread-safe pointer
  assignments (Fortran90 only).

- Converted the KPP user manual to Sphinx/ReadTheDocs format (this now
  replaces the prior ReadTheDocs documentaton).

- Added updates to allow KPP to be built on MacOS X systems.

- Added :program:`small_strato` C-I test that uses the exact same
  options as is described in :ref:`running-kpp-with-an-example-mechanism`.

.. _kpp240:

=========
KPP 2.4.0
=========

- Added new integrators: :file:`beuler.f90`, :file:`rosenbrock_mz.f90`,
  :file:`rosenbrock_posdef.f90`,  :file:`rosenbrock_posdef_h211b_qssa.f90`.

- Several memory sizes (:code:`MAX_EQN`, ...) have been increased to
  allow large chemical mechanisms.

- Added new Makefile target: :code:`list`.

- Added LaTeX User Manual.

- Now use :code:`ICNTRL(15)` to decide whether or not to toggle calling the
  :code:`Update_SUN`, :code:`Update_RCONST`, and :code:`Update_PHOTO`
  routines from within the integrator.

.. _kpp232gc:

============
KPP 2.3.2_gc
============

NOTE: Contains KPP Modifications specific to GEOS-Chem.

- Added workaround for F90 derived-type objects in inlined code
  (i.e. properly parse :code:`State_Het%xArea`, etc).

- Updated Github issue templates.

- :code:`MAX_INLINE` (max # of inlined code lines to read) has been
  increased to 200000.

- Commented out the :code:`Update_Sun()` functions in :code:`update_sun.F90`,
  :code:`update_sun.F`. (NOTE: These have been restored in
  :ref:`kpp250`).

- Default rate law functions are no longer written to :code:`gckpp_Rates.F90`.
  (NOTE: These have been restored in :ref:`kpp250`).

.. _kpp231gc:

============
KPP 2.3.1_gc
============

NOTE: KPP modifications specific to GEOS-Chem.

ALSO NOTE: ReadTheDocs documentation has been updated in :ref:`kpp250`
to remove GEOS-Chem specific information.

- Added documentation for ReadTheDocs.

- Added Github issue templates.

- README.md now contains the ReadTheDocs badge.

- README.md now points to kpp.readthedocs.io for documentation.

.. _kpp230gc:

============
KPP 2.3.0_gc
============

NOTE: Contains KPP modifications specific to GEOS-Chem.

- Added :file:`README.md` for the GC_updates branch.

- Added MIT license for the GC_updates branch.

- Add :code:`Aout` argument to return reaction rates from
  :code:`SUBROUTINE Fun`.

- Rename :file:`KPP/kpp_2.2.3_01` directory to :file:`KPP/kpp-code`.

- Now write :file:`gckpp_Model.F90` and :file:`gckpp_Precision.F90`
  from :code:`gen.c`.

- Do not write file creation & time to KPP-generated files (as this
  will cause Git to interpret each file as a new file to be added).

- Now create Fortran-90 source code files with :file:`*.F90` instead
  of :file:`*.f90`. (NOTE: In :ref:`kpp250`, this can specified with
  the :ref:`uppercasef90-cmd` command.)

- Remove calls to UPDATE_SUN and UPDATE_RCONST from all :code:`*.f90`
  integrators. (NOTE: This has been restored in :ref:`kpp250`.)

.. _kpp225gc:

============
KPP 2.2.5_gc
============

NOTE: Contains KPP modifications specific to GEOS-Chem.

- Increase :code:`MAX_INLINE` from 20000 to 50000

.. _kpp224gc:

============
KPP 2.2.4_gc
============

NOTE: Contains KPP modifications specific to GEOS-Chem.

- Add MIT license files for GC_updates branch and update
  :file:`README.md` accordingly

- Create :file:`README.md` for main branch

- Set :envvar:`FLEX_LIB_DIR` using :envvar:`FLEX_HOME` env variable if
  it is defined.

- Added an exponential integrator.

- Added array to :file:`*_Monitor` for family names
  (:code:`FAM_NAMES`) string vector.

- Added functionality for Prod/Loss families using :ref:`families` token.

- Add scripts necessary to build a new mechanism for GEOS-Chem

- Completed the prod/loss option (token: :code:`#FLUX [on/off]`)

- Added :code:`OMP THREADPRIVATE` to LinearAlgebra.F90

- Added :file:`rosenbrock_split.def` integrator definition

- Added :code:`OMPThreadPrivate` function for F77.

- Added declaration of :code:`A` in :ref:`Function`

- Added :code:`OMP THREADPRIVATE` Functionality to F90 output.

- Completed the split-form Function for F90.

- Increase maximum number of equations.

- Increase :code:`MAX_FAMILIES` parameter from 50 to 300

- Extend equation length limit to 200 characters.

- Also changed the species name for a family to the family name itself.

- Modified Families to minimize the number of additional species created

- Renamed and change indexing convention

- Removed unnecessary files

.. _kpp223:

=========
KPP 2.2.3
=========

- A new function called :code:`k_3rd_iupac` is available, calculating
  third-order rate coefficients using the formula used by IUPAC
  :cite:`Atkinson_et_al._2004`.

- While previous versions of KPP were using :program:`yacc` (yet another
  compiler compiler), the current version has been modified to be
  compatible with the parser generator :program:`bison`, which is the
  successor of :program:`yacc`.

- New Runge-Kutta integrators were added: :file:`kpp_dvode.f90`,
  :file:`runge_kutta.f90`, :file:`runge_kutta_tlm.f90`,
  :file:`sdirk_adj.f90`, :file:`sdirk_tlm.f90`.

- New Rosebrock method :code:`Rang3` was added.

- The new KPP command :command:`#DECLARE` was added (see:
  :ref:`declare-cmd`).

- Several vector and array functions from :program:`BLAS` (:code:`WCOPY`,
  :code:`WAXPY`, etc.) were replaced by Fortran90 expressions.

.. _kpp21:

=======
KPP 2.1
=======

- Described by :cite:t:`Sandu_and_Sander_2006`.

- Matlab is a new target language (see: :ref:`matlab-code`).

- The set of integrators has been extended with a general Rosenbrock
  integrator, and the corresponding tangent linear and adjoint methods.

- The KPP-generated Fortran90 code has a different file structure than
  the C or Fortran77 output (see: :ref:`f90-code`).

- An automatically generated Makefile facilitates the compilation of
  the KPP-generated code (see: :ref:`Makefile`).

- Equation tags provide a convenient way to refer to specific chemical
  reactions (see: :ref:`lookat-and-monitor`.

- The dummy index allows to test if a certain species occurs in the
  current chemistry mechanism. (see: :ref:`dummyindex-cmd`)

- Lines starting with :code:`//` are comment lines.

===================
KPP 1.1-f90-alpha12
===================

- First KPP version with Fortran90 :cite:p:`Sander_et_al._2005`.
