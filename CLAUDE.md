# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this repository is

KPP (Kinetic PreProcessor) is a source-code generator (written in C, built
with flex/bison) for chemical kinetics simulations. Given a chemical
mechanism description written in the KPP domain-specific language
(`.kpp`, `.spc`, `.eqn`, `.def` files), the `kpp` executable generates
standalone simulation code in Fortran90, C, or Matlab: ODE
function/Jacobian/Hessian evaluation, sparse-matrix data structures, a
numerical integrator, and a driver program.

KPP itself is not a chemistry solver — it is a compiler whose target
languages are Fortran90/C/Matlab and whose output is a complete,
compilable box-model simulation.

## Build

```bash
export KPP_HOME=/home/bob/repos/KPP3
export PATH=$PATH:$KPP_HOME/bin
export KPP_FLEX_LIB_DIR=/usr/lib   # wherever libfl.a/libfl.so lives

cd $KPP_HOME/src
make clean       # remove *.o (optional, only needed after a prior build)
make distclean   # also remove bin/kpp (optional)
make              # builds bin/kpp via gcc + flex + bison
```

Requires `gcc`, `flex`, `bison`, and `sed` on the PATH. Compiler/flags are
configured in `src/Makefile.defs` (defaults to `gcc`; honors `CC`/`FC`
from the shell environment if set). On macOS, `gcc` is aliased to
`clang` by default — see `docs/source/getting_started/installation.rst`
for the Homebrew GCC workaround, the flex-library-path workaround, and
the stack-size (`ulimit -s`) requirement.

`kpp` must be invoked with `KPP_HOME` set, since it looks up
`$KPP_HOME/models`, `$KPP_HOME/int`, and `$KPP_HOME/drv` by default
(overridable via `KPP_MODEL`, `KPP_INT`, `KPP_DRV`).

## Running a single mechanism through KPP

```bash
cd ci-tests/F90_small_strato
$KPP_HOME/bin/kpp F90_small_strato.kpp
make -f Makefile_F90_small_strato COMPILER=GFORTRAN
./F90_small_strato.exe
```

This generates `<ROOT>_Function.f90`, `<ROOT>_Jacobian*.f90`,
`<ROOT>_Global.f90`, `<ROOT>_Parameters.f90`, `<ROOT>_Rates.f90`,
`<ROOT>_Initialize.f90`, `<ROOT>_Integrator.f90`, `<ROOT>_Main.f90`,
etc., named after the `ROOT` in `#DRIVER`/model name.

## Tests (C-I tests)

There is no unit test framework — correctness is verified end-to-end by
generating code for real mechanisms, compiling, and running a box-model
simulation. Requires `bin/kpp` to be built first.

```bash
# Run the full C-I test suite (all mechanism/integrator/language combos)
$KPP_HOME/.ci-pipelines/ci-testing-script.sh | tee ci-tests.log

# Remove compiler-generated files from all C-I test folders afterward
$KPP_HOME/.ci-pipelines/ci-cleanup-script.sh
```

To run/inspect a single C-I test, see "Running a single mechanism"
above — each subdirectory of `ci-tests/` is one test, using the naming
convention `<LANG>_<description>` (e.g. `F90_rosadj`, `C_sd`). The list
of test directory names to run is maintained in
`GENERAL_TESTS`/`MCM_1`/`MCM_2`/`MINVERSION_TEST` in
`.ci-pipelines/ci-common-defs.sh` — new C-I tests must be added there
or they won't run in CI. GitHub Actions (`.github/workflows/run-ci-tests.yml`)
runs this same script across a matrix of GCC versions (9–16) on every
push/PR.

## Documentation

User docs are Sphinx/ReadTheDocs sources under `docs/source/`, built
from `docs/source/conf.py` (see `.readthedocs.yaml`). When changing
behavior (new command, new integrator, new output file), update the
corresponding `.rst` file — `docs/source/input/kpp_commands.rst` for
`#COMMAND` options, `docs/source/output/` for generated-file formats,
`docs/source/num_methods/` for integrators.

## Architecture

### Directory roles (each maps to a `KPP_*` search path or a fixed role)

- `src/` — the KPP program itself (C, flex/bison). Builds `bin/kpp`.
- `bin/` — holds the built `kpp` executable.
- `models/` — chemical mechanism definitions (`.def`/`.eqn`/`.spc`)
  usable via `#MODEL`; searched by default, overridable with `KPP_MODEL`.
- `int/` — integrator (numerical solver) templates, one `.def` +
  language-specific source (`.f90`/`.c`/`.m`) per integrator, selected
  via `#INTEGRATOR`; overridable with `KPP_INT`. `int/user_contributed/`
  holds community-contributed integrators.
- `drv/` — driver program templates (the `main`/box-model driver),
  selected via `#DRIVER`; overridable with `KPP_DRV`.
- `util/` — per-language function templates (Fortran90/C/Matlab
  variants of the same utility) substituted into generated code.
- `examples/` — example `.kpp` mechanism files.
- `ci-tests/` — one subdirectory per C-I test (see Tests above).
- `.ci-pipelines/` — shell scripts that run/clean the C-I tests.
- `site-lisp/` — `kpp.el`, an Emacs major mode for KPP files.
- `docs/` — Sphinx source for the ReadTheDocs user manual.

### How KPP itself works (`src/`)

1. **Scanner/parser** (`scan.l` → `lex.yy.c`, `scan.y` → `y.tab.c`,
   `scanner.c`, `scanutil.c`): flex/bison-based lexer+parser reads the
   `.kpp`/`.def`/`.eqn`/`.spc` input language and fills in-memory tables
   (atom list, species list, stoichiometry matrices, rate expressions,
   option list), declared in `gdata.h`. `kpp.c` is the driver/`main`.
2. **Species reordering**: variable vs. fixed species are separated and
   variable species are reordered with a Markowitz-type algorithm to
   minimize Jacobian LU fill-in (unless `#REORDER OFF`).
3. **Expression-tree / code generation core** (`gen.c`, `code.c`,
   `code.h`): builds a language-independent expression-tree
   representation of each generated statement (production/destruction
   terms, Jacobian, Hessian, stoichiometric matrix, etc.) directly from
   the coefficient matrices — one statement at a time, not one
   whole-program IR.
4. **Per-target-language emitters** (`code_f90.c`, `code_c.c`,
   `code_f77.c` [deprecated], `code_matlab.c`): each expression tree is
   rendered into valid syntax for the selected `#LANGUAGE`, then
   line-wrapped for readability.
5. Driver/integrator/util *template* files (in `drv/`, `int/`, `util/`)
   are passed through a substitution preprocessor (symbol substitution,
   not full macro expansion) rather than being generated from scratch.

Adding a new `#COMMAND` touches `scan.h`, `scan.l`, `scan.y`,
`scanner.c` in `src/`, plus a new C-I test and new `.rst` docs and
`site-lisp/kpp.el` — see
`docs/source/tech_info/adding_new_commands.rst` for the exact checklist.

### Versioning

The version string (`X.Y.Z`, semver, checked against `#MINVERSION` in
mechanism files) must be kept in sync in four places whenever it
changes: `src/gdata.h` (`KPP_VERSION`), `CHANGELOG.md`,
`docs/source/conf.py` (`release`), and the KPP Wikipedia page.

### macOS-specific constraint

macOS's default 64KB stack limit means `MAX_EQN`/`MAX_SPECIES` in
`src/gdata.h` are conditionally compiled smaller on macOS (`#ifdef
MACOS`) than on Linux; if a mechanism is too large to parse, these
constants (and possibly `ulimit -s`) need adjusting.
