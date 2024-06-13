# Changelog

<!-- Github markdown syntax: -->
<!-- https://docs.github.com/en/get-started/writing-on-github/getting-started-with-writing-and-formatting-on-github/basic-writing-and-formatting-syntax -->

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased] - TBD
### Added
- New integrator `implicit12`, with corresponding C-I test

## [3.1.1] - 2024-04-30
### Changed
- Updated Python package versions for ReadTheDocs in `docs/requirements.txt`
- Now request Python 3.12 for ReadTheDocs builds in `.readthedocs.yaml`
- Updated `MAX_NO_OF_LINES` to 2000 to parse the MCM mechanism
- Updaeed `MAX_EQN` to 18000 to parse the MCM mechanism

### Fixed
- Now only add tha extra `Aout` argument to `Fun` and `Fun_Split` for F90 (see issues #56, #96)

<!-- Version numbers must be synchronized in CHANGELOG.md, -->
<!-- src/gdata.h, and docs/source/conf.py-->
## [3.1.0] - 2023-12-20
### Added
- `#AUTOREDUCE` has been added to the list of KPP commands in the ReadTheDocs documentaton
- Added `examples/mcm` folder with minimal example for the Master Chemical Mechanism
- Added C-I test for MCM, based on the minimal example

### Removed
- `TRANSPORT` and `TRANSPORTALL` input options; these were obsolete
- `LUMP` input option; this was obsolete
- `DEFRAD`, `SETRAD`, `INITIALIZE`, `XGRID`, `YGRID`, `ZGRID`, `WRITE_OPT`, `RUN`, `USE`, `USES`; these were obsolete

## [3.0.2] - 2023-06-02
### Added
- Added `.readthedocs.yaml` for configuring ReadTheDocs builds
- Added a ReadTheDocs badge in `README.md`
- State that `#INTEGRATOR none` statements should be removed in RTD documentation

### Changed
- Added pybtex and pybtex-docutils to the `docs/requirements.txt`

### Fixed
- Minor issues in `site-lisp/kpp.el` file for Emacs

## [3.0.1] - 2023-03-21
### Added
- `CITATION.cff` file which will activate the "Cite this repository" option.

### Fixed
- Fixed a segmentation fault when using #STOICHMAT by using dynamically-sized variables `EqnNr` and `VarNr` rather than static variables `MAX_EQN` and `MAX_SPECIES`. 

## [3.0.0] - 2022-11-09

### Added

- C-language updates
  - Restored driver programs `general.c` and `general_adj.c`
  - Updated rate-law functions in `util/UserRateLaws.c`
  - Added C-I tests for C-language integrators (using the
    `small_strato` mechanism
  - Now use `//` instead of `/* */` comment strings
- ReadTheDocs documentation updates:
  - Authors in the reference page are now listed alphabetically
  - In-text citations now use the :cite:t style (e.g. "Smith et al [2000]")
  - Corrected several omissions
  - Removed table numbers from tables (to reduce confusion)
  - Now document additional installation steps for MacOS X
  - Added documentation about `KPP_FLEX_LIB_DIR`
  - Added badges for release, release date, and Azure pipelines to index.rst
  - All authors from the Lin et al [2022] paper about KPP 3.0.0 have
    been added to the front page of ReadTheDocs and .zenodo.json
- C-I test additions
  - Added a C-I test for the `rosenbrock_autoreduce` integrator

### Changed

- C-I test updates
  - Renamed CI-test folders for clarity (`C_*`, `F90_*`, 'X_*')
  - Abstracted reusable code into the `ci-common-defs.sh` script
  - Replaced `ci-manual-testing-script.sh` with `ci-testing-script.sh`
    (which now works both locally and on Azure Dev Pipelines)
  - Renamed `ci-manual-cleanup-script.sh` to `ci-cleanup-script.sh`
  - Updated docs accordingly
- Makefile changes
  - `Makefile.defs` now uses environment variable `KPP_FLEX_LIB_DIR`
    when it cannot find the flex library file in standard locations v
  - Removed host-specific if blocks from `Makefile.defs`
- Other changes
  - DOUBLE_COMPLEX is now replaced by COMPLEX(kind=dp)
  - Fixed incorrect license string in .zenodo.json
  - Added extra `free()` statements in `src/gen.c` to avoid memory
    leaks
- Bug fix in int/feuler.f90: FIX must be the 2nd argument passed to
  routine Fun (from routine FunTemplate)

### Removed
- bibtex.json (no longer needed w/ Sphinx 3.5.4)
- Fun() no longer uses Vdotout since it can be retrieved from `Vdot`

## [2.6.0] - 2022-06-24

### Added

- New integrator: rosenbrock_autoreduce.f90
  - Also added corresponding documentation for ReadTheDocs

### Changed

- The `ICNTRL(5)` option in the LSODE integrator is now `ICNTRL(10)`.

## [2.5.0] - 2022-05-18

### Added
- New integrators
  - feuler.f90
- New C-I tests
  - ros_minver
  - ros_upcase
  - Added scripts to run C-I tests manually (for development/debugging)
- Brought updates for GEOS-Chem into the main line of development
  - Added #MINVERSION switch to force KPP to stop unless a minimum
    version is used
  - Added #UPPERCASEF90 to generate F90 code with the .F90 suffixes
  - Removed EQUIVALENCE statements from F90-generated code; VAR and
    FIX now point to C witihin integrators.  This is to ensure
    thread-safe operation when using KPP-generated code in an OpenMP
    parallel environment.
- Bug fixes:
  - Make sure to inline parameter "sp" into the _Global.F90 file
    when the "#DOUBLE off" option is used.


### Changed
- Code updates
   - Rewrote code to remove compiler warnings
   - Fortran-90 makefiles now use GFORTRAN as the default compiler option
   - Added ReadTheDocs output
   - Updated top-of-file comment headers to point to the KPP Github
     site and to acknowledge new authors
   - Routine Fun() now returns optional arguments Aout and Vdotout
- Updates for building on MacOS
  - Reduce size of MAX_EQN and MAX_SPECIES to get KPP to run within
    65532 kb of stack memory

## [2.4.0] - 2022-04-25

### Added

- Brought updates from the MECCA branch into the main line of
  development:
  - new integrators: beuler.f90, rosenbrock_mz.f90, rosenbrock_posdef.f90,
    rosenbrock_posdef_h211b_qssa.f90
  - several memory sizes (MAX_EQN, ...) increased to allow large chemical mechanisms
  - new Makefile target: list
  - LaTeX User Manual added
- Now use ICNTRL(15) to decide whether or not to toggle calling the
  Update\_SUN, Update\_RCONST, and Update\_PHOTO routines from within
  the integrator

## [2.3.2_gc]

### Added

- Workaround for F90 derived-type objects in inlined code (i.e. properly parse State_Het%xArea, etc).
- Write global variables NUMDEN, MW, SR\_MW, SR\_TEMP, TEMP\_OVER\_K300, K300\_OVER\_TEMP to gckpp_Global.F90
- Documentation for ReadTheDocs (in the docs folder)
- Github issue templates

### Changed

- MAX_INLINE (max # of inlined code lines to read) is now 200000
- Version number in gdata.h is now 2.3.2
- README.md now contains the ReadTheDocs badge
- README.md now points to kpp.readthedocs.io for documentation

### Removed

- Comment out the Update\_Sun() functions in update\_sun.F90, update\_sun.F
- Default rate law functions are no longer written to gckpp_Rates.F90

## [2.3.1_gc]

### Added

- Documentation for ReadTheDocs (in the docs folder)
- Github issue templates

### Changed

- Version number in gdata.h is now 2.3.1
- README.md now contains the ReadTheDocs badge
- README.md now points to kpp.readthedocs.io for documentation

## [2.3.0_gc]

### Added

- Added README.md for the GC_updates branch
- Added MIT license for the GC_updates branch
- Add Aout argument to return reaction rates from SUBROUTINE Fun

### Changed

- Rename KPP/kpp\_2.2.3_01 folder to KPP/kpp-code
- Now write gckpp\_Model.F90 and gckpp\_Precision.F90 from gen.c
- Do not write file creation & time to KPP-generated files
- Now create Fortran-90 source code files with *.F90 instead of *.f90

### Removed

- Remove calls to UPDATE\_SUN and UPDATE\_RCONST from all *.f90 integrators

## [2.2.5_gc]

### Changed

- Increase MAX_INLINE from 20000 to 50000

## [2.2.4_gc]

### Added

- Add MIT license files and update README.md accordingly
- Create README.md for main branch
- Set FLEX_LIB_DIR using FLEX_HOME env variable if it is defined
- Added an exponential integrator
- Added array to *_Monitor for family names (FAM_NAMES) string vector
- Added functionality for Prod/Loss families using FAMILY token
- Add scripts necessary to build a new mechanism for GEOS-Chem
- Completed the prod/loss option (token: #FLUX [on/off])
- Added OMP THREADPRIVATE to LinearAlgebra.F90
- Added rosenbrock_split.def integrator definition
- Added OMPThreadPrivate function for F77
- Added declaration of "A" in *_Function.F90
- Added OMPThreadPrivate Functionality to F90 output
- Completed the split-form Function for f90

### Changed

- Increase maximum number of equations
- Increase MAX_FAMILIES parameter from 50 to 300
- Extend equation length limit to 200 characters.
- Also changed the species name for a family to the family name itself.
- Modified Families to minimize the number of additional species created
- Rename and change indexing convention

### Removed

- Remove unnecessary files
- Remove files for GEOS-Chem from top-level directory; these are now in the GEOS-Chem source code
- Delete the old name directory

### Fixed

- Fix to add coefficients to *_Monitor.f90
- Fixed a pad processing of terms leading to LOSS fam also generating PROD terms.
- Fix input argument to ComputeFamilies() from VAR to V - MSL
- Parameter NFAM was not defined. This is fixed. (MSL)
- Added fix to camculate NET rather than GROSS prod/loss families - MSL
- Add fix to calculate NET rather than GROSS prod/loss families - MSL
- Re-fixed the equivalence statement needed for GEOS-Chem to work. MSL
- Corrected miniterpretation of the split-form Function for f90.

## [2.2.3_rs3] - 2018-11-30

### CHANGES by Rolf Sander:
- several int/*.f90: TABs converted to spaces
- several dvode files converted to unix file format
- kpp/int/beuler.f90:
  - factor 101 added to Hstart
  - description of the error numbers IERR added
- kpp/int/rosenbrock_posdef_h211b_qssa.f90:
  USE KPP_ROOT_Function, ONLY: Fun_Split

### CHANGES by Sergey Gromov and Rolf Sander:
- several limits increased:
  - src/gdata.h:
    #define MAX_K 1000
  - src/scan.l:
    char crtToken[1000];
    char nextToken[1000];
    char crtFile[1000];
    char crt_rate[1000];
  - src/scan.y:
    %union{char str[1000];};

### CHANGES by Adrian Sandu:
- new integrator beuler.f90 (based on sdirk.f90)

## [2.2.3_rs2] - 2016-11-21

### CHANGES by Rolf Sander:
- src/gen.c: bug fix in GenerateParamHeader(): char name[MAX_SPNAME]
  instead of name[20]
- kpp.el updated: #REPLACE, #ENDREPLACE, uncertainty of rate coefficient
- Makefile.defs.Linux: compiler option -Wno-implicit-function-declaration

## [2.2.3_rs1] - 2016-08-28

### CHANGES by Rolf Sander:
- upgrade to the official KPP version 2.2.3

## [2.2.1_rs7] - 2016-02-09

### CHANGES by Rolf Sander:
- rosenbrock_mz.f90-orig deleted
- new function StoichNum to get stoichiometric numbers:
  - src/gen.c: new function GenerateFun_Split()
  - kp4.sh: GenerateStoichNum added to KPP_SUBROUTINE_LIST

## [2.2.1_rs6] - 2016-01-07

### CHANGES by Rolf Sander:
- src/code_f90.c: Maximum number of continuation lines MAX_NO_OF_LINES
  increased again, now to 250
- src/gdata.h: info adjusted
- src/gen.c: "USE kpp_Parameters" added to SUBROUTINE Initialize to
  allow species-specific code in #INLINE F90_INIT
- several int/*.f90: TABs converted to spaces

### CHANGES by Patrick Joeckel:
- int/rosenbrock_mz.f90: loop fusion according to BULL
- new file Makefile.m

Changes by Domenico Taraborrelli:
- option for faster chemistry calculation by alternative time stepping
  (integrator rosenbrock_posdef_h211b_qssa)
  - modified files:
    - src/gen.c: new function GenerateFun_Split
    - kp4.sh: GenerateFun_Split added to KPP_SUBROUTINE_LIST
  - new integrator files:
    rosenbrock_posdef_h211b_qssa.def, rosenbrock_posdef_h211b_qssa.f90

## [2.2.1_rs5] - since caaba_2.5n, 2010-02-22

### CHANGES by Rolf Sander:
- src/gdata.h: several limits increased to allow larger mechanisms:
  #define MAX_EQN       11000
  #define MAX_SPECIES    3500
  #define MAX_EQNTAG       32
  #define MAX_K           300
- src/gen.c:
  - crow and diag deleted because they are not used
  - Fortran90 double precision changed to SELECTED_REAL_KIND(12,307)
- src/Makefile: PHONY target list added to list the configuration
- src/scan.h: MAX_INLINE increased to 100000
- src/scan.l: several limits increased to allow larger mechanisms:
  char crtToken[300];
  char nextToken[300];
  char crtFile[300];
  char crt_rate[300];
- src/scan.y:
  - "char" increased to 300 (must be the same as MAX_K)
- util/sutil.f90: TABs converted to spaces

## [2.2.1_rs4] - since caaba_2.5m, 2009-11-27

### CHANGES by Patrick Joeckel:
- Makefile, Makefile.defs, Makefile.defs.*, src/Makefile: changes to
  allow compilation on different machines

### CHANGES by Rolf Sander:
- src/code_matlab.c: "#include <time.h>" added because it is necessary
  for time_t, see: http://en.wikipedia.org/wiki/Time_t
- src/scanutil.c:
  - "#include <malloc.h>" removed because it comes from <stdlib.h> (this
    had caused a problem on MAC-OSX)
- src/scan.y:
  - "#include <malloc.h>" removed because it comes from <stdlib.h> (this
    had caused a problem on MAC-OSX)
- using GNU Bison 2.3 now instead of 2.1

## [2.2.1_rs3] - since caaba_2.5c, 2008-07-17

### CHANGES by Patrick Joeckel and Rolf Sander:
- src/code_f90.c and src/gdata.h: changes in MAX_EQNTAG etc. to avoid
  problems with long equation tags

## [2.2.1_rs2] - since caaba_2.4, 2008-02-18

### CHANGES by Rolf Sander:
- kpp_lsode.f90, kpp_seulex.f90, rosenbrock.f90, rosenbrock_posdef.f90,
  runge_kutta.f90, and sdirk.f90: IERR_NAMES (for error messages)
  adapted to correct syntax
- src/code_f90.c:
  - Maximum number of continuation lines increased to 100. If
    MAX_NO_OF_LINES is too small, kpp may split lines incorrectly.

## [2.2.1_rs] - 2007-08-22

### CHANGES by Rolf Sander:
- code_f90.c: added FlushBuf() to F90_DeclareData (otherwise MAX_OUTBUF
  would have to be very large for large reaction mechanisms)
- gen.c: The declaration of RTOLS was deleted because it is not needed
  by the integrators. If the driver programs need it, they can define it
  themselves.
- Makefile: new option maintainer-clean is now consistent with
  src/Makefile

### CHANGES by Astrid Kerkweg:
- src/code_f90.c: string length of F90_types and in subroutine
  F90_DeclareData increased from 12 to 32 to avoid problems with long
  species names and long equation tags

### CHANGES by Adrian Sandu and Rolf Sander:
- It looks like the model is hypersensitive to negative concentrations;
  many times when small negative concentrations are produced the entire
  future trajectory is put in jeopardy. In the new integrator
  rosenbrock_posdef.f90, this is fixed by changing "CALL
  WCOPY(N,Ynew,1,Y,1)" to "Y = MAX(Ynew,ZERO)".

## [2.2.1] - 2006-12-01

### CHANGES by Adrian Sandu:
- rosenbrock_soa deleted
- util/sutil.c: new subroutine KppDecompCmplxR
- new files: examples/cell.kpp, examples/saprc2006.kpp,
  int/runge_kutta.c, and int/sdirk.c

### CHANGES by Rolf Sander:
- int/kpp_lsode.f90: like all other integrators, kpp_lsode now also
  returns IERR==1 after successful completion

### CHANGES by Adrian Sandu and Rolf Sander:
- src/Makefile, scan.l, scan.y: yacc replaced by bison, implementing the
  bug fixes suggested by Jason Lander

## [2.2_July2006] - 2006-07-23

### CHANGES by Adrian Sandu:
- models/CMAQ added
- changes in int/runge_kutta_adj.f90 and int/runge_kutta_tlm.f90
- util/blas.c: new subroutines Set2Zero and WADD
- util/blas.f90: new subroutines WGEFA and WGESL
- util/sutil.f90: new subroutines KppDecompCmplxR, KppSolveCmplxR, and
  KppSolveTRCmplxR

## [2.2.June2006] - 2006-06-06

### CHANGES by Philipp Miehe and Adrian Sandu:
- new integrators kpp_sdirk4, rosenbrock_soa, runge_kutta, and sdirk
- integrators rosenbrock, rosenbrock_tlm, and rosenbrock_adj: completely
  revised
- new kpp command #DECLARE [SYMBOL|VALUE]
- util/blas.f90: new subroutine WADD
- util/sutil.f90: new subroutines KppSolveTRIndirect and KppSolveTRCmplx
- several files examples/* added
- changes in drv/general_adj.f90 and drv/general_tlm.f90
- examples/mimi* deleted

## [2.1] - 2005-07-19
