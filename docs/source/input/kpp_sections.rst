.. _kpp-sections:

############
KPP sections
############

A :literal:`#` sign at the beginning of a line followed by a section
name starts a new KPP section. Then a list of items separated by
semicolons follows. A section ends when another KPP section or command
occurs, i.e. when another :literal:`#` sign occurs at the beginning of
a line. The syntax of an item definition is different for each
particular section.

.. _atoms:

======
#ATOMS
======

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

======
#CHECK
======

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

===================
#DEFVAR and #DEFFIX
===================

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

==========
#EQUATIONS
==========

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

=========
#FAMILIES
=========

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

===========
#INITVALUES
===========

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

====================
#LOOKAT and #MONITOR
====================

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

===================
#SETVAR and #SETFIX
===================

The commands :command:`#SETVAR` and :command:`#SETFIX` change the type of an
already defined species. Then, depending on the integration method,
one may or may not use the initial classification, or can easily move
one species from one category to another. The use of the generic
species :code:`VAR_SPEC`, :code:`FIX_SPEC`, and :code:`ALL_SPEC` is
also allowed. Examples for these sections are:

.. code-block:: console

   #SETVAR ALL_SPEC;
   #SETFIX H2O; CO2;
