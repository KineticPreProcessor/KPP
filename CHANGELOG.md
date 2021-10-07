All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [2.3.2_gc]

### Added

- Workaround for F90 derived-type objects in inlined code (i.e. properly parse State_Het%xArea, etc).
- Write global variables NUMDEN, MW, SR_MW, SR_TEMP, TEMP_OVER_K300, K300_OVER_TEMP to gckpp_Global.F90
- Documentation for ReadTheDocs (in the docs folder)
- Github issue templates

### Changed

- MAX_INLINE (max # of inlined code lines to read) is now 200000
- Version number in gdata.h is now 2.3.2
- README.md now contains the ReadTheDocs badge
- README.md now points to kpp.readthedocs.io for documentation

### Removed

- Comment out the Update_Sun() functions in update_sun.F90, update_sun.F
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

- Rename KPP/kpp_2.2.3_01 folder to KPP/kpp-code
- Now write gckpp_Model.F90 and gckpp_Precision.F90 from gen.c
- Do not write file creation & time to KPP-generated files
- Now create Fortran-90 source code files with *.F90 instead of *.f90

### Removed

- Remove calls to UPDATE_SUN and UPDATE_RCONST from all *.f90 integrators


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
