.. _input-for-kpp:

#############
Input for KPP
#############

KPP basically handles two types of input files: **Kinetic description
files** and **auxiliary files**. Kinetic description files are in KPP
syntax and described in the following sections. Auxiliary files are
described in the section entitled
:ref:`auxiliary-files-and-the-substitution-preprocessor`.

KPP kinetic description files specify the chemical equations, the
initial values of each of the species involved, the integration
parameters, and many other options. The KPP preprocessor parses the
kinetic description files and generates several output files. Files
that are written in KPP syntax have one of the suffixes :file:`.kpp`,
:file:`.spc`, :file:`.eqn`, or :file:`def`.

The following general rules define the structure of a kinetic
description file:

-  A KPP program is composed of :ref:`kpp-sections`,
   :ref:`kpp-commands`, and :ref:`inlined-code`. Their syntax is
   presented in :ref:`bnf-description`.

-  Comments are either enclosed between the curly braces ":code:`{`"
   and ":code:`}`", or written in a line starting with two slashes and
   a space "// ".

-  Any name given by the user to denote an atom or a species is
   restricted to be less than 32 character in length and can only
   contain letters, numbers, or the underscore character. The first
   character cannot be a number. All names are case insensitive.

The kinetic description files contain a detailed specification of the
chemical model, information about the integration method and the desired
type of results. KPP accepts only one of these files as input, but using
the :ref:`include-cmd` command, code from separate files can be
combined. The include files can be nested up to 10 levels. KPP will
parse these files as if they were a single big file. By carefully
splitting the chemical description, KPP can be configured for a broad
range of users. In this way the users can have direct access to that
part of the model that they are interested in, and all the other details
can be hidden inside several include files. Often, the atom definitions
(:file:`atoms.kpp`) are included first, then species definitions
(:file:`*.spc`), and finally the equations of the chemical mechanism
(:file:`*.eqn`).

.. _kpp-sections:

============
KPP sections
============

A :literal:`#` sign at the beginning of a line followed by a section
name starts a new KPP section. Then a list of items separated by
semicolons follows. A section ends when another KPP section or command
occurs, i.e. when another :literal:`#` sign occurs at the beginning of
a line. The syntax of an item definition is different for each
particular section.

.. _atoms:

#ATOMS
------

The atoms that will be further used to specify the components of a
species must be declared in an :command:`#ATOMS` section, e.g.:

.. code-block:: console

   #ATOMS N; O; Na; Br;

Usually, the names of the atoms are the ones specified in the periodic
table of elements. For this table there is a predefined file containing
all definitions that can be used by the command:

.. code-block:: console

   #INCLUDE atoms.kpp

This should be the first line in a KPP input file, because it allows to
use any atom in the periodic table of elements throughout the kinetic
description file.

.. _check:

#CHECK
------

KPP is able to do mass balance checks for all equations. Some chemical
equations are not balanced for all atoms, and this might still be
correct from a chemical point of view. To accommodate for this, KPP can
perform mass balance checking only for the list of atoms specified in
the :command:`#CHECK` section, e.g.:

.. code-block:: console

   #CHECK N; C; O;

The balance checking for all atoms can be enabled by using the
:command:`#CHECKALL` command. Without :command:`#CHECK` or
:command:`#CHECKALL`, no checking is performed. The :literal:`IGNORE`
atom can also be used to control mass balance checking.

.. _defvar-and-deffix:

#DEFVAR and #DEFFIX
-------------------

There are two ways to declare new species together with their atom
composition: :command:`#DEFVAR` and :command:`#DEFFIX`. These sections
define all the species that will be used in the chemical mechanism.
Species can be variable or fixed. The type is implicitly specified by
defining the species in the appropriate sections. A fixed species does
not vary through chemical reactions.

For each species the user has to declare the atom composition. This
information is used for mass balance checking.  To ignore mass balance
checking for a given species, one can declare the predefined atom
:command:`IGNORE` as being part of the species composition. Examples
for these sections are:

.. code-block:: console

   #DEFVAR
     NO2 = N + 2O;
     CH3OOH = C + 4H + 2O;
     HSO4m = IGNORE;
     RCHO = IGNORE;
   #DEFFIX
     CO2 = C + 2O;

.. _equations:

#EQUATIONS
----------

The chemical mechanism is specified in the :command:`#EQUATIONS`
section. Each equation is written in the natural way in which a
chemist would write it:

.. code-block:: console

   #EQUATIONS

   <R1> NO2 + hv = NO + O3P :  6.69e-1*(SUN/60.0);
   <R2> O3P + O2 + AIR = O3 :  ARR_ac(5.68e-34,  -2.80);
   <R3> O3P + O3 = 2O2 :       ARR_ab(8.00e-12, 2060.0);
   <R4> O3P + NO + AIR = NO2 : ARR_ac(1.00e-31,  -1.60);
   //... etc ...

.. note::

   The above example is taken from the :command:`saprc99` mechanism
   (see :file:`models/saprc99.eqn`), with some whitespace deleted for
   clarity.  Optional :ref:`equation tags <eqntags-cmd>` are specified
   by text within :code:`< >` angle brackets.  Functions that compute
   **saprc99** equation rates (e.g. :code:`ARR_ac`,
   :code:`ARR_ab`) are defined in :file:`util/UserRateLaws.f90`
   and :file:`util/UserRateLawsInterfaces.f90`.

Only the names of already defined species can be used. The rate
coefficient has to be placed at the end of each equation, separated by a
colon. The rate coefficient does not necessarily need to be a numerical
value. Instead, it can be a valid expression (or a call to an
:ref:`inlined rate law function <inlined-code>`) in the :ref:`target
language <language-cmd>`.  If there are several :command:`#EQUATIONS`
sections in the input, their contents will be concatenated.

A minus sign in an equation shows that a species is consumed in a
reaction but it does not affect the reaction rate. For example, the
oxidation of methane can be written as:

.. code-block:: console

   CH4 + OH = CH3OO + H2O - O2 : k_CH4_OH;

However, it should be noted that using negative products may lead to
numerical instabilities.

Often, the stoichiometric factors are integers. However, it is also
possible to have non-integer yields, which is very useful to
parameterize organic reactions that branch into several side reactions:

.. code-block:: console

   CH4 + O1D = .75 CH3O2 + .75 OH + .25 HCHO + 0.4 H + .05 H2 : k_CH4_O1D;

KPP provides two pre-defined dummy species: :literal:`hv` and
:literal:`PROD`. Using dummy species does not affect the numerics of
the integrators. It only serves to improve the readability of the
equations. For photolysis reactions, :literal:`hv` can be specified as
one of the reagents to indicate that light (:math:`h\nu`) is needed for this
reaction, e.g.:

.. code-block:: console

   NO2 + hv = NO + O : J_NO2;

When the products of a reaction are not known or not important, the
dummy species :literal:`PROD` should be used as a product. This is
necessary because the KPP syntax does not allow an empty list of
products. For example, the dry deposition of atmospheric ozone to the
surface can be written as:

.. code-block:: console

   O3 = PROD : v_d_O3;

The same equation must not occur twice in the :command:`#EQUATIONS`
section. For example, you may have both the gas-phase reaction of
:literal:`N2O5` with water in your mechanism and also the
heterogeneous reaction on aerosols:

.. code-block:: console

   N2O5 + H2O = 2 HNO3 : k_gas;
   N2O5 + H2O = 2 HNO3 : k_aerosol;

These reactions must be merged by adding the rate coefficients:

.. code-block:: console

   N2O5 + H2O = 2 HNO3 : k_gas + k_aerosol;

.. _families:

#FAMILIES
---------

Chemical families (for diagnostic purposes) may be specified in the
:command:`#FAMILIES` section as shown below.  Family names beginning
with a :code:`P` denote production, and those beginning with an
:code:`L` denote loss.

.. code-block:: console

   #FAMILIES
     POx : O3 + NO2 + 2NO3 + HNO3 + ... etc. add more species as needed ...
     LOx : O3 + NO2 + 2NO3 + HNO3 + ... etc. add more species as needed ...
     PCO : CO;
     LCO : CO;
     PSO4 : SO4;
     LCH4 : CH4;
     PH2O2 : H2O2;

KPP will examine the chemical mechanism and create a dummy species for
each defined family.  Each dummy species will archive the production
and loss for the family.  For example, each molecule of CO that is
produced will be added to the :code:`PCO` dummy species.  Likewise,
each molecule of CO that is consumed will be added to the :code:`LCO`
dummy species. This will allow the :code:`PCO` and :code:`LCO` species
to be later archived for diagnostic purposes. Dummy species for chemical
families will not be included as active species in the mechanism.

.. _initvalues:

#INITVALUES
-----------

The initial concentration values for all species can be defined in the
:command:`#INITVALUES` section, e.g.:

.. code-block:: console

   #INITVALUES
     CFACTOR = 2.5E+19;
     NO2 = 1.4E-9;
     CO2 = MyCO2Func();
     ALL_SPEC = 0.0;

If no value is specified for a particular species, the default value
zero is used. One can set the default values using the generic species
names: :code:`VAR_SPEC`, :code:`FIX_SPEC`, and :code:`ALL_SPEC`. In order
to use coherent units for concentration and rate coefficients, it is
sometimes necessary to multiply each value by a constant factor. This
factor can be set by using the generic name :code:`CFACTOR`. Each of
the initial values will be multiplied by this factor before being
used. If :code:`CFACTOR` is omitted, it defaults to one.

The information gathered in this section is used to generate the
:code:`Initialize` subroutine (cf  :ref:`Initialize`). In more complex 3D
models, the initial values are usually taken from some input files or
some global data structures. In this case, :command:`#INITVALUES` may
not be needed.

.. note::

   If you are building your mechanism in Fortran 90, note that using
   floating point constants such as :literal:`2.5E+19` will only give
   you a single precision value (which will incur roundoff error after
   the 6th or 7th decimal place.

   If you are using double precision, you should use Fortran double
   precision exponents (e.g. :literal:`2.5D+19`) instead.  This will
   give you a true double precision value.

.. _lookat-and-monitor:

#LOOKAT and #MONITOR
--------------------

There are two sections in this category: :command:`#LOOKAT` and
:command:`#MONITOR`.

The section instructs the preprocessor what are the species for which
the evolution of the concentration, should be saved in a data file. By
default, if no :command:`#LOOKAT` section is present, all the species
are saved. If an atom is specified in the :command:`#LOOKAT` list then
the total mass of the particular atom is reported. This allows to
check how the mass of a specific atom was conserved by the integration
method. The :command:`#LOOKATALL` command can be used to specify all
the species. Output of :command:`#LOOKAT` can be directed to the file
:file:`ROOT.dat` using the utility subroutines described in the
section entitled :ref:`Util`.

The :command:`#MONITOR` section defines a different list of species
and atoms. This list is used by the driver to display the
concentration of the elements in the list during the integration. This
may give us a feedback of the evolution in time of the selected
species during the integration. The syntax is similar to the
:command:`#LOOKAT` section. With the driver :code:`general`,
output of :command:`#MONITOR` goes to the screen (STDOUT). The order
of the output is: first variable species, then fixed species, finally
atoms. It is not the order in the :command:`MONITOR` command.

Examples for these sections are:

.. code-block:: console

   #LOOKAT NO2; CO2; O3; N;
   #MONITOR O3; N;

.. _setvar-and-setfix:

#SETVAR and #SETFIX
-------------------

The commands :command:`#SETVAR` and :command:`#SETFIX` change the type of an
already defined species. Then, depending on the integration method,
one may or may not use the initial classification, or can easily move
one species from one category to another. The use of the generic
species :code:`VAR_SPEC`, :code:`FIX_SPEC`, and :code:`ALL_SPEC` is
also allowed. Examples for these sections are:

.. code-block:: console

   #SETVAR ALL_SPEC;
   #SETFIX H2O; CO2;

.. _kpp-commands:

============
KPP commands
============

A KPP command begins on a new line with a :code:`#` sign, followed by a
command name and one or more parameters.  Details about each command
are given in the following subsections.

.. _table-cmd-defaults:

.. list-table:: Default values for KPP commands
   :align: center
   :header-rows: 1

   * - KPP command
     - default value
     - KPP command
     - default value
   * - :command:`#AUTOREDUCE`
     - :code:`OFF`
     - :command:`#INTFILE`
     -
   * - :command:`#CHECKALL`
     -
     - :command:`#JACOBIAN`
     - :code:`SPARSE_LU_ROW`
   * - :command:`#DECLARE`
     - :code:`SYMBOL`
     - :command:`#LANGUAGE`
     -
   * - :command:`#DOUBLE`
     - :code:`ON`
     - :command:`#LOOKATALL`
     -
   * - :command:`#DRIVER`
     - :code:`none`
     - :command:`#MEX`
     - :code:`ON`
   * - :command:`#DUMMYINDEX`
     - :code:`OFF`
     - :command:`#MINVERSION`
     -
   * - :command:`#EQNTAGS`
     - :code:`OFF`
     - :command:`#MODEL`
     -
   * - :command:`#FUNCTION`
     - :code:`AGGREGATE`
     - :command:`#REORDER`
     - :code:`ON`
   * - :command:`#HESSIAN`
     - :code:`ON`
     - :command:`#STOCHASTIC`
     - :code:`OFF`
   * - :command:`#INCLUDE`
     -
     - :command:`#STOICMAT`
     - :code:`ON`
   * - :command:`#INTEGRATOR`
     -
     - :command:`#UPPERCASEF90`
     - :code:`OFF`

.. _autoreduce-cmd:

#AUTOREDUCE
-----------

The :command:`#AUTOREDUCE ON` command can be used with
:command:`#INTEGRATOR rosenbrock` to enable
:ref:`automatic  mechanism reduction <rosenbrock-autoreduce>` as
described in :cite:t:`Lin_et_al._2022`.  Automatic mechanism reduction
is disabled by default.

.. _declare-cmd:

#DECLARE
--------

The :command:`#DECLARE` command determines how constants like
:code:`dp`, :code:`NSPEC`, :code:`NVAR`, :code:`NFIX`, and
:code:`NREACT` are inserted into the KPP-generated code.
:command:`#DECLARE SYMBOL` (the default) will declare array variables
using parameters from the :ref:`Parameters` file. :command:`#DECLARE VALUE`
will replace each parameter with its value.

For example, the global array variable :code:`C` is declared in the
:ref:`Global` file generated by KPP.  In the :program:`small_strato`
example (described in :ref:`running-kpp-with-an-example-mechanism`),
:code:`C` has dimension :code:`NSPEC=7`. Using  :command:`#DECLARE
SYMBOL` will generate the following code in :ref:`Global`:

.. code-block:: fortran

   ! C - Concentration of all species
     REAL(kind=dp), TARGET :: C(NSPEC)

Whereas :command:`#DECLARE VALUE` will generate this code instead:

.. code-block:: fortran

   ! C - Concentration of all species
     REAL(kind=dp), TARGET :: C(7)

We recommend using :command:`#DECLARE SYMBOL`, as most modern compilers
will automatically replace each parameter (e.g. :code:`NSPEC`) with its
value (e.g :code:`7`). However, if you are using a very old compiler
that is not as sophisticated, :command:`#DECLARE VALUE` might result in
better-optmized code.

.. _double-cmd:

#DOUBLE
-------

The :command:`#DOUBLE` command selects single or double precision
arithmetic. :command:`ON` (the default) means use double precision,
:command:`OFF` means use single precision (see the section entitled
:ref:`Precision`).

.. important::

   We recommend using double precision whenever possible.  Using
   single precision may lead to integration non-convergence errors
   caused by roundoff and/or underflow.

.. _driver-cmd:

#DRIVER
-------

The :command:`#DRIVER` command selects the driver, i.e., the file from
which the main function is to be taken. The parameter is a file name,
without suffix. The appropriate suffix (:code:`.f90`, :code:`.F90`,
:code:`.c`, or :code:`.m`) is automatically appended.

Normally, KPP tries to find the selected driver file in the directory
:file:`$KPP_HOME/drv/`. However, if the supplied file name contains a slash,
it is assumed to be absolute. To access a driver in the current
directory, the prefix :file:`./` can be used, e.g.:

.. code-block:: console

   #DRIVER ./mydriver

It is possible to choose the empty dummy driver :command:`none`, if the
user wants to include the KPP generated modules into a larger model
(e.g. a general circulation or a chemical transport model) instead of
creating a stand-alone version of the chemical integrator. The driver
:command:`none` is also selected when the :command:`#DRIVER` command
is missing. If the command occurs twice, the second replaces the first.

.. _dummyindex-cmd:

#DUMMYINDEX
-----------

It is possible to declare species in the :ref:`defvar-and-deffix`
sections that are not used in the :ref:`equations` section. If your
model needs to check at run-time if a certain species is included in
the current mechanism, you can set to :command:`#DUMMYINDEX ON`. Then,
KPP will set the indices to zero for all species that do not occur in
any reaction. With :command:`#DUMMYINDEX OFF` (the default), those are
undefined variables. For example, if you frequently switch between
mechanisms with and without sulfuric acid, you can use this code:

.. code-block:: fortran

   IF (ind_H2SO4=0) THEN
     PRINT *, 'no H2SO4 in current mechanism'
   ELSE
     PRINT *, 'c(H2SO4) =', C(ind_H2SO4)
   ENDIF

.. _eqntags-cmd:

#EQNTAGS
--------

Each reaction in the :ref:`equations` section may start with an
equation tag which is enclosed in angle brackets, e.g.:

.. code-block:: console

   <R1> NO2 + hv = NO + O3P :  6.69e-1*(SUN/60.0);

With :command:`#EQNTAGS` set to :command:`ON`, this equation tag can be
used to refer to a specific equation (cf. :ref:`monitor`). The default
for :command:`#EQNTAGS` is :command:`OFF`.

.. _function-cmd:

#FUNCTION
---------

The :command:`#FUNCTION` command controls which functions are generated
to compute the production/destruction terms for variable
species. :command:`AGGREGATE` generates one function that computes the
normal derivatives. :command:`SPLIT` generates two functions
for the derivatives in production and destruction forms.

.. _hessian-cmd:

#HESSIAN
--------

The option :command:`ON` (the default) of the :command:`#HESSIAN`
command turns the Hessian generation on (see section
:ref:`Hessian-and-HessianSP`). With :command:`OFF` it is switched off.

.. _include-cmd:

#INCLUDE
--------

The :command:`#INCLUDE` command instructs KPP to look for the file
specified as a parameter and parse the content of this file before
proceeding to the next line. This allows the atoms definition, the
species definition and the equation definition to be shared between
several models. Moreover this allows for custom configuration of KPP to
accommodate various classes of users. Include files can be either in one
of the KPP directories or in the current directory.

.. _integrator-cmd:

#INTEGRATOR
-----------

The :command:`#INTEGRATOR` command selects the integrator definition
file. The parameter is the file name of an integrator, without
suffix. The effect of

.. code-block:: console

   #INTEGRATOR integrator_name

is similar to:

.. code-block:: console

   #INCLUDE $KPP_HOME/int/integrator_name.def

The :command:`#INTEGRATOR` command allows the use of different
integration techniques on the same model. If it occurs twice, the second
replaces the first. Normally, KPP tries to find the selected integrator
files in the directory :file:`$KPP_HOME/int/`. However, if the supplied
file name contains a slash, it is assumed to be absolute. To access an
integrator in the current directory, the prefix :file:`./` can be used,
e.g.:

.. code-block:: console

   #INTEGRATOR ./mydeffile

.. _intfile-cmd:

#INTFILE
--------

.. attention::

   :command:`#INTFILE` is used internally by KPP but should not be used
   by the KPP user. Using :ref:`integrator-cmd` alone suffices to
   specify an integrator.

The integrator definition file selects an integrator file with
:command:`#INTFILE` and also defines some suitable options for it. The
:command:`#INTFILE` command selects the file that contains the integrator
routine. The parameter of the
command is a file name, without suffix. The appropriate suffix
(:code:`.f90`, :code:`.F90`, :code:`.c`, or :code:`.m` is appended and
the result selects the file from which the integrator
is taken. This file will be copied into the code file in the appropriate
place.

.. _jacobian-cmd:

#JACOBIAN
---------

The :command:`#JACOBIAN` command controls which functions are generated
to compute the Jacobian. The option :command:`OFF` inhibits the
generation of the Jacobian routine. The option :command:`FULL` generates
the Jacobian as a square :code:`NVAR x NVAR` matrix. It should only be
used if the integrator needs the whole Jacobians. The options
:command:`SPARSE_ROW` and :command:`SPARSE_LU_ROW` (the default) both
generate the Jacobian in sparse (compressed on rows) format. They should
be used if the integrator needs the whole Jacobian, but in a sparse
form. The format used is compressed on rows. With
:command:`SPARSE_LU_ROW`, KPP extends the number of nonzeros to account
for the fill-in due to the LU decomposition.

.. _language-cmd:

#LANGUAGE
---------

.. attention::

   The :command:`Fortran77` language option is deprecated in
   :ref:`kpp250` and  later versions. All further KPP development will
   only support Fortran90.

The :command:`#LANGUAGE` command selects the target language in which the
code file is to be generated. Available options are :command:`Fortran90`,
:command:`C`, or :command:`matlab`.

You can select the suffix (:code:`.F90` or :code:`.f90`) to use for
Fortran90 source code generated by KPP (cf. :ref:`uppercasef90-cmd`).

.. _mex-cmd:

#MEX
----

:program:`Mex` is a Matlab extension that allows
to call functions written in Fortran and C directly from within the
Matlab environment. KPP generates the mex interface routines for the
ODE function, Jacobian, and Hessian, for the target languages C,
Fortran77, and Fortran90. The default is :command:`#MEX ON`. With
:command:`#MEX OFF`, no Mex files are generated.

.. _minversion-cmd:

#MINVERSION
-----------

You may restrict a chemical mechanism to use a given version of KPP or
later. To do this, add

.. code-block:: console

   #MINVERSION X.Y.Z

to the definition file.

The version number (:code:`X.Y.Z`) adheres to the Semantic
Versioning style (https://semver.org), where :code:`X` is the major
version number, :code:`Y` is the minor version number, and :code:`Z` is the
bugfix (aka “patch”) version number.

For example, if :command:`#MINVERSION 2.4.0` is specified, then KPP will
quit with an error message unless you are using KPP 2.4.0 or later.

.. _model-cmd:

#MODEL
------

The chemical model contains the description of the atoms, species, and
chemical equations. It also contains default initial values for the
species and default options including a suitable integrator for the
model. In the simplest case, the main kinetic description file, i.e. the
one passed as parameter to KPP, can contain just a single line selecting
the model. KPP tries to find a file with the name of the model and the
suffix :file:`.def` in the :file:`$KPP_HOME/models` subdirectory. This
file is then parsed. The content of the model definition file is written
in the KPP language. The model definition file points to a species file
and an equation file. The species file includes further the atom
definition file. All default values regarding the model are
automatically selected. For convenience, the best integrator and driver
for the given model are also automatically selected.

The :command:`#MODEL` command is optional, and intended for using a
predefined model. Users who supply their own reaction mechanism do not
need it.

.. _reorder-cmd:

#REORDER
--------

Reordering of the species is performed in order to minimize the fill-in
during the LU factorization, and therefore preserve the sparsity
structure and increase efficiency. The reordering is done using a
diagonal Markowitz algorithm. The details are explained in
:cite:t:`Sandu_et_al._1996`. The default is :command:`ON`.
:command:`OFF` means that KPP does not reorder the species. The order
of the variables is the order in which the species are
declared in the :command:`#DEFVAR` section.

.. _stochastic-cmd:

#STOCHASTIC
-----------

The option :command:`ON` of the :command:`#STOCHASTIC` command turns
on the generation of code for stochastic kinetic simulations (see the
section entitled :ref:`Stochastic`.  The default option is :command:`OFF`.

.. _stoicmat-cmd:

#STOICMAT
---------

Unless the :command:`#STOICMAT` command is set to :command:`OFF`, KPP
generates code for the stoichiometric matrix, the vector of reactant
products in each reaction, and the partial derivative of the time
derivative function with respect to rate coefficients
(cf. :ref:`Stoichiom-and-StoichiomSP`).

.. _checkall-lookatall-cmd:

#CHECKALL, #LOOKATALL
---------------------

The shorthand commands :command:`#CHECKALL` and :command:`#LOOKATALL`
apply :command:`#CHECK` and :command:`#LOOKAT`, respectively, to all
species in the mechanism.

.. _uppercasef90-cmd:

#UPPERCASEF90
-------------

If you have selected :command:`#LANGUAGE Fortran90` option, KPP will
generate source code ending in :code:`.f90` by default. Setting
:command:`#UPPERCASEF90 ON` will tell KPP to generate Fortran90 code
ending in :code:`.F90` instead.

.. _inlined-code:

============
Inlined Code
============

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

List of inlined types
---------------------

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

F90_DATA
--------

This inline type was introduced in a previous version of KPP to
initialize variables. It is now obsolete but kept for compatibility. For
Fortran90, :command:`F90_GLOBAL` should be used instead.

.. _f90-global:

F90_GLOBAL
----------

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

F90_INIT
--------

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

F90_RATES
---------

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

F90_RCONST
----------

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

F90_RCONST_USE
--------------

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

F90_UTIL
--------

This inline type can be used to define utility subroutines.

.. _auxiliary-files-and-the-substitution-preprocessor:

=================================================
Auxiliary files and the substitution preprocessor
=================================================

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

List of auxiliary files for Fortran90
--------------------------------------

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

List of symbols replaced by the substitution preprocessor
---------------------------------------------------------

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

.. _icntrl-rcntrl:

=================================================================
Controlling the Integrator with :code:`ICNTRL` and :code:`RCNTRL`
=================================================================

In order to offer more control over the integrator, KPP provides the
arrays :code:`ICNTRL` (integer) and :code:`RCNTRL` (real). Each of them
is an array of 20 elements that allow the fine-tuning of the integrator.
All integrators (except for :code:`tau_leap` and :code:`gillespie`) use
:code:`ICNTRL` and :code:`RCNTRL`. Details can be found in the comment
lines of the individual integrator files in :code:`$KPP_HOME/int/`.

ICNTRL
------

.. _table-icntrl:

.. table:: Summary of ICNTRL usage in the f90 integrators.
           Here, Y = used, and s = solver-specific usage.
   :align: center

   +------------------------+---+---+---+---+---+---+---+---+---+----+----+----+----+----+----+----+----+
   | ICNTRL                 | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9 | 10 | 11 | 12 | 13 | 14 | 15 | 16 | 17 |
   +========================+===+===+===+===+===+===+===+===+===+====+====+====+====+====+====+====+====+
   | beuler                 |   | Y | Y | Y | Y | Y | s |   |   |    |    |    |    |    | Y  |    |    |
   +------------------------+---+---+---+---+---+---+---+---+---+----+----+----+----+----+----+----+----+
   | dvode                  |   |   |   |   |   |   |   |   |   |    |    |    |    |    | Y  |    |    |
   +------------------------+---+---+---+---+---+---+---+---+---+----+----+----+----+----+----+----+----+
   | exponential            |   |   |   |   |   |   |   |   |   |    |    |    |    |    |    |    |    |
   +------------------------+---+---+---+---+---+---+---+---+---+----+----+----+----+----+----+----+----+
   | feuler                 |   |   |   |   |   |   |   |   |   |    |    |    |    |    | Y  | Y  | Y  |
   +------------------------+---+---+---+---+---+---+---+---+---+----+----+----+----+----+----+----+----+
   | gillespie              |   |   |   |   |   |   |   |   |   |    |    |    |    |    |    |    |    |
   +------------------------+---+---+---+---+---+---+---+---+---+----+----+----+----+----+----+----+----+
   | lsode                  |   | Y |   | Y |   |   |   |   |   | s  |    |    |    |    | Y  |    |    |
   +------------------------+---+---+---+---+---+---+---+---+---+----+----+----+----+----+----+----+----+
   | radau5                 |   | Y |   | Y | Y | Y |   |   |   |    | Y  |    |    |    | Y  |    |    |
   +------------------------+---+---+---+---+---+---+---+---+---+----+----+----+----+----+----+----+----+
   | rosenbrock_adj         | Y | Y | Y | Y |   | s | s | s |   |    |    |    |    |    | Y  |    |    |
   +------------------------+---+---+---+---+---+---+---+---+---+----+----+----+----+----+----+----+----+
   | rosenbrock             | Y | Y | Y | Y |   |   |   |   |   |    |    |    |    |    | Y  | Y  |    |
   +------------------------+---+---+---+---+---+---+---+---+---+----+----+----+----+----+----+----+----+
   | rosenbrock_tlm         | Y | Y | Y | Y |   |   |   |   |   |    |    | s  |    |    | Y  |    |    |
   +------------------------+---+---+---+---+---+---+---+---+---+----+----+----+----+----+----+----+----+
   | rosenbrock_autoreduce  | Y | Y | Y | Y |   |   |   |   |   |    |    | s  | s  | s  | Y  | Y  |    |
   +------------------------+---+---+---+---+---+---+---+---+---+----+----+----+----+----+----+----+----+
   | runge_kutta_adj        |   | Y | Y | Y | Y | s | s | s | s | s  | Y  |    |    |    | Y  |    |    |
   +------------------------+---+---+---+---+---+---+---+---+---+----+----+----+----+----+----+----+----+
   | runge_kutta            |   | Y | Y | Y | Y | Y |   |   |   | s  | Y  |    |    |    | Y  |    |    |
   +------------------------+---+---+---+---+---+---+---+---+---+----+----+----+----+----+----+----+----+
   | runge_kutta_tlm        |   | Y | Y |   | Y | Y | s |   | s | s  | Y  | s  |    |    | Y  |    |    |
   +------------------------+---+---+---+---+---+---+---+---+---+----+----+----+----+----+----+----+----+
   | sdirk4                 |   | Y |   | Y |   |   |   |   |   |    |    |    |    |    | Y  |    |    |
   +------------------------+---+---+---+---+---+---+---+---+---+----+----+----+----+----+----+----+----+
   | sdirk_adj              |   | Y | Y | Y | Y | Y | s | s |   |    |    |    |    |    | Y  |    |    |
   +------------------------+---+---+---+---+---+---+---+---+---+----+----+----+----+----+----+----+----+
   | sdirk                  |   | Y | Y | Y | Y | Y |   |   |   |    |    |    |    |    | Y  |    |    |
   +------------------------+---+---+---+---+---+---+---+---+---+----+----+----+----+----+----+----+----+
   | sdirk_tlm              |   | Y | Y | Y | Y | Y | s |   | s |    |    | s  |    |    | Y  |    |    |
   +------------------------+---+---+---+---+---+---+---+---+---+----+----+----+----+----+----+----+----+
   | seulex                 | Y | Y |   | Y |   |   |   |   |   | s  | s  | s  | s  | s  | Y  |    |    |
   +------------------------+---+---+---+---+---+---+---+---+---+----+----+----+----+----+----+----+----+
   | tau_leap               |   |   |   |   |   |   |   |   |   |    |    |    |    |    |    |    |    |
   +------------------------+---+---+---+---+---+---+---+---+---+----+----+----+----+----+----+----+----+

ICNTRL(1)
~~~~~~~~~

Specifies the time-dependence of the function :math:`F` to be integrated:

.. table::
   :align: left

   +------------+----------------------------------------+
   | ICNTRL(1)  | Time dependence of :math:`F`           |
   +============+========================================+
   | 0          | :math:`F = F(t,y)` aka non-autonomous  |
   +------------+----------------------------------------+
   | 1          | :math:`F = F(y)` aka autonomous        |
   +------------+----------------------------------------+

ICNTRL(2)
~~~~~~~~~

Specifies the dimensionality of the absolute (:code:`ATOL`) and
relative (:code:`RTOL`) tolerances.  These can be expressed by
either a scalar or individually for each species in a vector.

.. table::
   :align: left

   +------------+-------------------------------------------+
   | ICNTRL(2)  | Dimensionality of ``ATOL`` and ``RTOL``   |
   +============+===========================================+
   | 0          | :code:`NVAR` -dimensional vector          |
   +------------+-------------------------------------------+
   | 1          | scalar                                    |
   +------------+-------------------------------------------+

ICNTRL(3)
~~~~~~~~~

Selects a specific integration method.

.. table::
   :align: left

   +-----------------------+------------+-----------------------+
   | Integrator            | ICNTRL(3)  | Method selected       |
   +=======================+============+=======================+
   | rosenbrock            | 0 or 4     | Rodas3 (default)      |
   |                       +------------+-----------------------+
   | rosenbrock_adj        | 1          | Ros2                  |
   |                       +------------+-----------------------+
   | rosenbrock_autoreduce | 2          | Ros3                  |
   |                       +------------+-----------------------+
   | rosenbrock_tlm        | 3          | Ros4                  |
   |                       +------------+-----------------------+
   |                       | 5          | Rodas4                |
   |                       +------------+-----------------------+
   |                       | 6          | Rang                  |
   |                       +------------+-----------------------+
   |                       | 7          | For future use        |
   +-----------------------+------------+-----------------------+
   | runge_kutta           | 0 or 1     | Radau-2A (default)    |
   |                       +------------+-----------------------+
   | runge_kutta_adj       | 2          | Lobatto-3C            |
   |                       +------------+-----------------------+
   | runge_kutta_tlm       | 3          | Gauss                 |
   |                       +------------+-----------------------+
   |                       | 4          | Radau-1A              |
   |                       +------------+-----------------------+
   |                       | 5          | Lobatto-3A            |
   +-----------------------+------------+-----------------------+
   | sdirk                 | 0 or 1     | Sdirk-2A (default)    |
   |                       +------------+-----------------------+
   | sdirk_adj             | 2          | Sdirk-2B              |
   |                       +------------+-----------------------+
   | sdirk_tlm             | 3          | Sdirk-3A              |
   |                       +------------+-----------------------+
   |                       | 4          | Sdirk-4A              |
   |                       +------------+-----------------------+
   |                       | 5          | Sdirk-4B              |
   |                       +------------+-----------------------+
   |                       | 6          | Backward Euler        |
   +-----------------------+------------+-----------------------+


ICNTRL(4)
~~~~~~~~~

Specifies the maximum number of integration steps.

ICNTRL(5)
~~~~~~~~~

Specifies the maximum number of Newton iterations.

ICNTRL(6)
~~~~~~~~~

Selects integrator-specific settings.

.. table::
   :align: left

   +-----------------+------------+-----------------------------------+
   | Integrator      | ICNTRL(6)  | Option selected                   |
   +=================+============+===================================+
   | rosenbrock_adj  | 0 thru 6   | Selection of a particular         |
   |                 |            | Rosenbrock method for the         |
   |                 |            | continuous adjoint integration    |
   |                 |            | (see ``ICNTRL(3)``)               |
   +-----------------+------------+-----------------------------------+
   | radau5          | 0          | Starting values for Newton        |
   |                 |            | iterations are interpolated       |
   |                 |            | (default)                         |
   | runge_kutta_adj +------------+-----------------------------------+
   |                 | 1          | Starting values for Newton        |
   | runge_kutta_tlm |            | iterations are zero               |
   |                 |            |                                   |
   | sdirk           |            |                                   |
   |                 |            |                                   |
   | sdirk_adj       |            |                                   |
   |                 |            |                                   |
   | sdirk_tlm       |            |                                   |
   +-----------------+------------+-----------------------------------+

ICNTRL(7)
~~~~~~~~~

Selects options for adjoint integrators.

.. table:: Selection of adjoint algorithm
   :align: left

   +-----------------+------------+----------------------------------+
   | Integrator      | ICNTRL(7)  | Adjoint algorithm selected       |
   +=================+============+==================================+
   | rosenbrock_adj  | 0 or 2     | Discrete adjoint with method     |
   |                 |            | ``ICNTRL(3)`` (default)          |
   |                 +------------+----------------------------------+
   |                 | 1          | No adjoint                       |
   |                 +------------+----------------------------------+
   |                 | 3          | Fully adaptive continous adjoint |
   |                 |            | with method ``ICNTRL(6)``        |
   |                 +------------+----------------------------------+
   |                 | 4          | Simplified continuous adjoint    |
   |                 |            | with method ``ICNTRL(6)``        |
   +-----------------+------------+----------------------------------+

.. table:: Method to solve the linear Adj equations
   :align: left

   +-----------------+------------+----------------------------------+
   | Integrator      | ICNTRL(7)  | Method selected                  |
   +=================+============+==================================+
   | runge_kutta_adj | 0 or 1     | Modified Newton re-using LU      |
   |                 |            | (default)                        |
   |                 +------------+----------------------------------+
   |                 | 2          | Direct solution (additional one  |
   |                 |            | LU factorizsation of 3Nx3N       |
   |                 |            | matrix per step; good for        |
   |                 |            | debugging                        |
   |                 +------------+----------------------------------+
   |                 | 3          | Adaptive solution (if Newton     |
   |                 |            | does not converge, switch to     |
   |                 |            | direct)                          |
   +-----------------+------------+----------------------------------+
   | sdirk_adj       | 0          | Modified Newton re-using LU      |
   |                 |            | (default)                        |
   | sdirk_tlm       +------------+----------------------------------+
   |                 | 1          | Direct solution (additional one  |
   |                 |            | LU factorizsation per stage)     |
   +-----------------+------------+----------------------------------+

ICNTRL(8)
~~~~~~~~

Determines if LU factorization will be checkpointed at each step
(for adjoint integrators only).

.. table::
   :align: left

   +-----------------+------------+------------------------------------+
   | Integrator      | ICNTRL(8)  | Option selected                    |
   +=================+============+====================================+
   | rosenbrock_adj  | 0          | Do not save LU factorization at    |
   |                 |            | each step (default)                |
   | runge_kutta_adj +------------+------------------------------------+
   |                 | 1          | Save LU factorization at each step |
   | sdirk_adj       |            |                                    |
   +-----------------+------------+------------------------------------+

ICNTRL(9)
~~~~~~~~~

Selects options for adjoint and tangent linear method (TLM)
integrators.

.. table:: Selection of adjoint algorithm
   :align: left

   +-----------------+------------+----------------------------------+
   | Integrator      | ICNTRL(9)  | Adjoint algorithm selected       |
   +=================+============+==================================+
   | runge_kutta_adj | 0 or 2     | Discrete adjoint with method     |
   |                 |            | ``ICNTRL(3)`` (default)          |
   |                 +------------+----------------------------------+
   |                 | 1          | No adjoint                       |
   |                 +------------+----------------------------------+
   |                 | 3          | Fully adaptive continous adjoint |
   |                 |            | with method ``ICNTRL(6)``        |
   |                 +------------+----------------------------------+
   |                 | 4          | Simplified continuous adjoint    |
   |                 |            | with method ``ICNTRL(6)``        |
   +-----------------+------------+----------------------------------+

.. table:: Selection of tangent linear method (TLM) error
	      estimation strategy
   :align: left

   +-----------------+------------+-----------------------------------+
   | Integrator      | ICNTRL(9)  | Strategy selected                 |
   +=================+============+===================================+
   | runge_kutta_tlm | 0          | Base number of iterations as      |
   |                 |            | forward solution                  |
   | sdirk_tlm       +------------+-----------------------------------+
   |                 | 1          | USE ``ATOL_tlm`` and ``RTOL_tlm`` |
   |                 |            | to calculate error estimation     |
   |                 |            | for TLM at Newton stages          |
   +-----------------+------------+-----------------------------------+

ICNTRL(10)
~~~~~~~~~~

Selects integrator-specific options.

.. table::
   :align: left

   +-----------------+------------+-----------------------------------+
   | Integrator      | ICNTRL(10) | Option selected                   |
   +=================+============+===================================+
   | lsode           | user       | Maximum order of the integration  |
   |                 | supplied   | formula allowed                   |
   +-----------------+------------+-----------------------------------+
   | runge_kutta     | 0          | Error estimation strategy:        |
   |                 |            | one additional stage at ``c=0``   |
   | runge_kutta_adj |            | (default)                         |
   |                 +------------+-----------------------------------+
   | runge_kutta_tlm | 1          | Error estimation strategy:        |
   |                 |            | Two additional stages at ``c=0``  |
   |                 |            | and SDIRK at ``c=1``, stiffly     |
   |                 |            | accurate                          |
   |                 +------------+-----------------------------------+
   |                 |            |                                   |
   +-----------------+------------+-----------------------------------+
   | seulex          | 0          | No verbose output                 |
   |                 +------------+-----------------------------------+
   |                 | 1          | Dense verbose output              |
   +-----------------+------------+-----------------------------------+

ICNTRL(11)
~~~~~~~~~~

Selects integrator-specific settings.

.. table::
   :align: left

   +-----------------+------------+-----------------------------------+
   | Integrator      | ICNTRL(11) | Option selected                   |
   +=================+============+===================================+
   | seulex          | user       | Maximum number of columns in the  |
   |                 | supplied   | extrapolation. Default is 12.     |
   +-----------------+------------+-----------------------------------+
   | radau5          | 0          | Gustaffson step size control      |
   |                 +------------+-----------------------------------+
   |                 | 1          | Classical step size control       |
   | runge_kutta     |            |                                   |
   |                 |            |                                   |
   | runge_kutta_adj |            |                                   |
   |                 |            |                                   |
   | runge_kutta_tlm |            |                                   |
   |                 |            |                                   |
   |                 |            |                                   |
   +-----------------+------------+-----------------------------------+

ICNTRL(12)
~~~~~~~~~~

Selects integrator-specific settings.

.. table::
   :align: left

   +-----------------------+------------+-----------------------------------+
   | Integrator            | ICNTRL(12) | Option selected                   |
   +=======================+============+===================================+
   | rosenbrock_autoreduce | 0          | Disable mechanism auto-reduction  |
   |                       |            | (i.e. acts in the same way as the |
   |                       |            | ``rosenbrock`` integrator)        |
   |                       +------------+-----------------------------------+
   |                       | 1          | Enables mechanism auto-reduction, |
   |                       |            | set threshold in ``RCNTRL(`12)``  |
   |                       +------------+-----------------------------------+
   |                       |            |                                   |
   +-----------------------+------------+-----------------------------------+
   | rosenbrock_tlm        | 0          | TLM truncation error is not used  |
   |                       +------------+-----------------------------------+
   | runge_kutta_tlm       | 1          | TLM truncation error is computed  |
   |                       |            | and used                          |
   | sdirk_tlm             |            |                                   |
   |                       +------------+-----------------------------------+
   |                       |            |                                   |
   +-----------------------+------------+-----------------------------------+
   | seulex                | 0          | Nsequence = 2 (default)           |
   |                       +------------+-----------------------------------+
   |                       | 1          | Nsequence =                       |
   |                       |            | 1,2,3,4,6,8,12,16,24,32,48,...    |
   |                       +------------+-----------------------------------+
   |                       | 2          | Nsequence =                       |
   |                       |            | 2,3,4,6,8,12,16,24,32,48,64,...   |
   |                       +------------+-----------------------------------+
   |                       | 3          | Nsequence =                       |
   |                       |            | 1,2,3,4,5,6,7,8,9,10,...          |
   |                       +------------+-----------------------------------+
   |                       | 4          | Nsequence =                       |
   |                       |            | 2,3,4,5,6,7,8,9,10,11,...         |
   +-----------------------+------------+-----------------------------------+

ICNTRL(13)
~~~~~~~~~~

Selects integrator-specific settings.

.. table::
   :align: left

   +-----------------------+------------+-----------------------------------+
   | Integrator            | ICNTRL(13) | Option selected                   |
   +=======================+============+===================================+
   | rosenbrock_autoreduce | 0          | In auto-reduction, disables       |
   |                       |            | scanning species P and L rates    |
   |                       |            | throughout the internal timesteps |
   |                       |            | of the integrator.                |
   |                       +------------+-----------------------------------+
   |                       | 1          | In auto-reduction, ensables       |
   |                       |            | scanning species P and L rates    |
   |                       |            | throughout the internal timesteps |
   |                       |            | of the integrator, for            |
   |                       |            | repartitioning.                   |
   +-----------------------+------------+-----------------------------------+
   | seulex                | 0, 1       | Sets the ``lambda`` parameter     |
   |                       |            | for verbose output.               |
   +-----------------------+------------+-----------------------------------+

ICNTRL(14)
~~~~~~~~~~

(Solver-specific for :code:`rosenbrock_autoreduce`) If set to
:code:`> 0`, then the threshold is calculated based on the max of
production and loss rate of the species ID specified in
:code:`ICNTRL(14)` multiplied by :code:`RCNTRL(14)`.

ICNTRL(15)
~~~~~~~~~~

Determines which :code:`Update_*` subroutines are called within the
integrator.

.. list-table::
   :align: left
   :header-rows: 1

   * - ICNTRL(15)
     - Option selcted
   * - -1
     - Do not call any :code:`Update_*` subroutines
   * - 0
     - Use the integrator-specific default values
   * - 1
     - Call  :code:`Update_RCONST` from within the integrator.
   * - 2
     - Call :code:`Update_PHOTO` from within the integrator
   * - 3
     - Call :code:`Update_RCONST` and :code:`Update_PHOTO` from within
       the integrator.
   * - 4
     - Call :code:`Update_SUN` from within the integrator
   * - 5
     - Call :code:`Update_SUN` and :code:`Update_RCONST` from within
       the integrator
   * - 6
     - Call :code:`Update_SUN` and :code:`Update_PHOTO` from within
       the integrator.
   * - 7
     - Call :code:`Update_SUN`, :code:`Update_PHOTO`, and
       :code:`Update_RCONST` from within the integrator.

Calling :code:`Update_RCONST` may be necessary when reaction rate
coefficients depend on the concentration of a specific species, e.
g.:

.. code-block:: console

   HSO3m + HSO5m + Hp = 2 HSO4m + Hp : k_aqueous( C(ind_Hp) );

This ensures that the concentration :code:`C(ind_Hp)` at the specific
integration time is used when the reaction rate coefficient is
updated within the integrator.

ICNTRL(16)
~~~~~~~~~~

Specifies how negative values should be handled.

.. list-table::
   :align: left
   :header-rows: 1

   * - ICNTRL(16)
     - Option selcted
   * - 0
     - Leave negative values unchanged
   * - 1
     - Set negative values to zero
   * - 2
     - Print warning and continue
   * - 3
     - Print error message and stop

ICNTRL(17)
~~~~~~~~~~

Selects the amount of verbose output that will be generated.

.. list-table::
   :align: left
   :header-rows: 1

   * - ICNTRL(17)
     - Option selcted
   * - 0
     - Only return error number
   * - 1
     - Full verbose error output

ICNTRL(18)
~~~~~~~~~~

Currently not used.

ICNTRL(19)
~~~~~~~~~~

Currently not used.

ICNTRL(20)
~~~~~~~~~~

Currently not used.

RCNTRL
------

.. _table-rcntrl:

.. table:: Summary of RCNTRL usage in the f90 integrators.
           Here, Y = used, and s = solver-specific usage.
   :align: center

   +------------------------+---+---+---+---+---+---+---+---+---+----+----+----+----+----+----+----+----+----+----+
   | RCNTRL                 | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9 | 10 | 11 | 12 | 13 | 14 | 15 | 16 | 17 | 18 | 19 |
   +========================+===+===+===+===+===+===+===+===+===+====+====+====+====+====+====+====+====+====+====+
   | beuler                 | Y | Y | Y | Y | Y | Y | Y | Y | Y | Y  | Y  |    |    |    |    |    |    |    |    |
   +------------------------+---+---+---+---+---+---+---+---+---+----+----+----+----+----+----+----+----+----+----+
   | dvode                  |   |   |   |   |   |   |   |   |   |    |    |    |    |    |    |    |    |    |    |
   +------------------------+---+---+---+---+---+---+---+---+---+----+----+----+----+----+----+----+----+----+----+
   | exponential            |   |   |   |   |   |   |   |   |   |    |    |    |    |    |    |    |    |    |    |
   +------------------------+---+---+---+---+---+---+---+---+---+----+----+----+----+----+----+----+----+----+----+
   | feuler                 |   |   |   |   |   |   |   |   |   |    |    |    |    |    |    |    |    |    |    |
   +------------------------+---+---+---+---+---+---+---+---+---+----+----+----+----+----+----+----+----+----+----+
   | gillespie              |   |   |   |   |   |   |   |   |   |    |    |    |    |    |    |    |    |    |    |
   +------------------------+---+---+---+---+---+---+---+---+---+----+----+----+----+----+----+----+----+----+----+
   | lsode                  | Y | Y | Y |   |   |   |   |   |   |    |    |    |    |    |    |    |    |    |    |
   +------------------------+---+---+---+---+---+---+---+---+---+----+----+----+----+----+----+----+----+----+----+
   | radau5                 |   | Y |   | Y | Y | Y | Y | Y | Y | Y  | Y  |    |    |    |    |    |    |    |    |
   +------------------------+---+---+---+---+---+---+---+---+---+----+----+----+----+----+----+----+----+----+----+
   | rosenbrock_adj         | Y | Y | Y | Y | Y | Y | Y |   |   |    |    |    |    |    |    |    |    |    |    |
   +------------------------+---+---+---+---+---+---+---+---+---+----+----+----+----+----+----+----+----+----+----+
   | rosenbrock             | Y | Y | Y | Y | Y | Y | Y |   |   |    |    |    |    |    |    |    |    |    |    |
   +------------------------+---+---+---+---+---+---+---+---+---+----+----+----+----+----+----+----+----+----+----+
   | rosenbrock_tlm         | Y | Y | Y | Y | Y | Y | Y |   |   |    |    |    |    |    |    |    |    |    |    |
   +------------------------+---+---+---+---+---+---+---+---+---+----+----+----+----+----+----+----+----+----+----+
   | rosenbrock_autoreduce  | Y | Y | Y | Y | Y | Y | Y |   |   |    |    | s  |    | s  |    |    |    |    |    |
   +------------------------+---+---+---+---+---+---+---+---+---+----+----+----+----+----+----+----+----+----+----+
   | runge_kutta_adj        | Y | Y | Y | Y | Y | Y | Y | Y | Y | Y  | Y  |    |    |    |    |    |    |    |    |
   +------------------------+---+---+---+---+---+---+---+---+---+----+----+----+----+----+----+----+----+----+----+
   | runge_kutta            | Y | Y | Y | Y | Y | Y | Y | Y | Y | Y  | Y  |    |    |    |    |    |    |    |    |
   +------------------------+---+---+---+---+---+---+---+---+---+----+----+----+----+----+----+----+----+----+----+
   | runge_kutta_tlm        | Y | Y | Y | Y | Y | Y | Y | Y | Y | Y  | Y  |    |    |    |    |    |    |    |    |
   +------------------------+---+---+---+---+---+---+---+---+---+----+----+----+----+----+----+----+----+----+----+
   | sdirk4                 | Y | Y | Y | Y | Y | Y | Y | Y | Y | Y  | Y  |    |    |    |    |    |    |    |    |
   +------------------------+---+---+---+---+---+---+---+---+---+----+----+----+----+----+----+----+----+----+----+
   | sdirk_adj              | Y | Y | Y | Y | Y | Y | Y | Y | Y | Y  | Y  |    |    |    |    |    |    |    |    |
   +------------------------+---+---+---+---+---+---+---+---+---+----+----+----+----+----+----+----+----+----+----+
   | sdirk                  | Y | Y | Y | Y | Y | Y | Y | Y | Y | Y  | Y  |    |    |    |    |    |    |    |    |
   +------------------------+---+---+---+---+---+---+---+---+---+----+----+----+----+----+----+----+----+----+----+
   | sdirk_tlm              | Y | Y | Y | Y | Y | Y | Y | Y | Y | Y  | Y  |    |    |    |    |    |    |    |    |
   +------------------------+---+---+---+---+---+---+---+---+---+----+----+----+----+----+----+----+----+----+----+
   | seulex                 | Y | Y | Y | Y | Y | Y | Y | Y |   | s  | s  | s  | s  | s  | s  | s  | s  | s  | s  |
   +------------------------+---+---+---+---+---+---+---+---+---+----+----+----+----+----+----+----+----+----+----+
   | tau_leap               |   |   |   |   |   |   |   |   |   |    |    |    |    |    |    |    |    |    |    |
   +------------------------+---+---+---+---+---+---+---+---+---+----+----+----+----+----+----+----+----+----+----+

RCNTRL(1)
~~~~~~~~~

:code:`Hmin`, the lower bound of the integration step size. It is
not recommended to change the default value of zero.

RCNTRL(2)
~~~~~~~~~

:code:`Hmax`, the upper bound of the integration step size.

RCNTRL(3)
~~~~~~~~~

:code:`Hstart`, the starting value of the integration step size.

RCNTRL(4)
~~~~~~~~~

:code:`FacMin`, lower bound on step decrease factor.

RCNTRL(5)
~~~~~~~~~

:code:`FacMax`, upper bound on step increase factor.

RCNTRL(6)
~~~~~~~~~

:code:`FacRej`, step decrease factor after multiple rejections.

RCNTRL(7)
~~~~~~~~~

:code:`FacSafe`, the factor by which the new step is slightly
   smaller than the predicted value.

RCNTRL(8)
~~~~~~~~~

:code:`ThetaMin`. If the Newton convergence rate is smaller than
ThetaMin, the Jacobian is not recomputed.

RCNTRL(9)
~~~~~~~~~

:code:`NewtonTol`, the stopping criterion for Newton’s method.

RCNTRL(10)
~~~~~~~~~~

:code:`Qmin`

RCNTRL(11)
~~~~~~~~~~

:code:`Qmax`. If :code:`Qmin < Hnew/Hold < Qmax`, then the step
size is kept constant and the LU factorization is reused.

RCNTRL(12)
~~~~~~~~~~

(Solver-specific for :code:`rosenbrock_autoreduce`) Used to specify
the threshold for auto-reduction partitioning, if :code:`ICNTRL(12) = 1`,
and :code:`ICNTRL(14) = 0`. Will be ignored if :code:`ICNTRL(14) > 0`.

RCNTRL(14)
~~~~~~~~~~

(Solver-specific for :code:`rosenbrock_autoreduce`) Used to specify
the multiplier for threshold for auto-reduction partitioning, if
:code:`ICNTRL(12) = 1`, and :code:`ICNTRL(14) > 0`, :code:`RCNTRL(14)`
is multiplied against max  of production and loss rates of species
:code:`ICNTRL(14)` to produce the partitioning threshold, ignoring
:code:`RCNTRL(12)`.

RCNTRL(15) - RCNTRL(19)
~~~~~~~~~~~~~~~~~~~~~~~

Solver-specific settings for :code:`seulex`.

RCNTRL(20)
~~~~~~~~~~

currently not used.
