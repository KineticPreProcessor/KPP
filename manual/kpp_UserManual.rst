.. container:: center

   | |image|
   | **KPP- User’s Manual**
   | *The Kinetic PreProcessor KPP.*
   | *An Environment for the*
   | *Simulation of Chemical Kinetic Systems*
   | **Adrian Sandu\ :math:`^1` & Rolf Sander\ :math:`^2`
     **((TODO: UPDATE LIST OF AUTHORS))****
   | :math:`^1` Department of Computer Science
   | Virginia Polytechnic Institute and State University
   | Blacksburg, Virginia 24060, USA
   | `sandu@cs.vt.edu <sandu@cs.vt.edu>`__
   | :math:`^2` Air Chemistry Department
   | Max-Planck Institute of Chemistry
   | PO Box 3060, 55020 Mainz, Germany
   | `rolf.sander@mpic.de <rolf.sander@mpic.de>`__

   This manual is available at:
   https://github.com/KineticPreProcessor/KPP/blob/main/manual/kpp_UserManual.pdf

.. container:: center

   Date: 2022-05-04

.. _`sec:install`:

Installation
============

This section can be skipped if KPP is already installed on your system.
If you work under Linux, you can probably use the precompiled executable
file that is in the directory of the distribution. Then you only have to
define the ``$KPP_HOME`` environment variable.

#. Define the ``$KPP_HOME`` environment variable to point to the
   complete path where KPP is installed. Also, add the path of the KPP
   executable to the ``$PATH`` environment variable. If, for example,
   KPP is installed in ``$HOME/kpp``, under the C shell you have to edit
   the file ``$HOME/.cshrc`` and add:

   ::

      setenv KPP_HOME $HOME/kpp
      setenv PATH ${PATH}:$KPP_HOME/bin

   If you use the bash shell, edit and add:

   ::

      export KPP_HOME=$HOME/kpp
      export PATH=$PATH:$KPP_HOME/bin

   After editing or , start a new shell to make sure these changes are
   in effect.

#. Make sure that is installed on your machine. Type to test this.

#. Make sure that is installed on your machine. Type to test this.

#. Make sure that the fast lexical analyzer generator is installed on
   your machine. Type to test this. Enter the path where the flex
   library ( or ) is located into , e.g. .

#. Change to the KPP directory:

   ::

      cd $KPP_HOME

#. To clean the KPP installation, delete the KPP object files and all
   the examples with:

   ::

      make clean

   To delete the KPP executable as well, type:

   ::

      make distclean

#. If necessary, edit and enter the name of your C compiler. The default
   setting is .

#. Create the kpp executable with:

   ::

      make

.. _`sec:model`:

Running KPP with an Example Stratospheric Mechanism
===================================================

Here we consider as an example a very simple Chapman-like mechanism for
stratospheric chemistry:

.. math::

   \begin{aligned}
   \chem{O_2}                 & \TOHV & 2 \chem{O}\\
   \chem{O} + \chem{O_2}      & \TO   & \chem{O_3}\\
   \chem{O_3}                 & \TOHV & \chem{O} + \chem{O_2}\\
   \chem{O} + \chem{O_3}      & \TO   & 2 \chem{O_2}\\
   \chem{O_3}                 & \TOHV & \chem{O(^1D)} + \chem{O_2}\\
   \chem{O(^1D)} + \chem{M}   & \TO   & \chem{O} + \chem{M}\\
   \chem{O(^1D)} + \chem{O_3} & \TO   & 2 \chem{O_2}\\
   \chem{NO} + \chem{O_3}     & \TO   & \chem{NO_2} + \chem{O_2}\\
   \chem{NO_2} + \chem{O}     & \TO   & \chem{NO} + \chem{O_2}\\
   \chem{NO_2}                & \TOHV & \chem{NO} + \chem{O}\end{aligned}

We use the mechanism with the purpose of illustrating the KPP
capabilities. However, the software tools are general and can be applied
to virtually any kinetic mechanism.

We focus on Fortran90. Particularities of the C, Fortran77, and Matlab
languages are discussed in Sections `4.2 <#sec:c>`__,
`4.3 <#sec:f77>`__, `4.4 <#sec:matlab>`__, respectively.

The KPP input files (with suffix ) specify the model, the target
language, the precision, the integrator and the driver, etc. The file
name (without the suffix ) serves as the root name for the simulation.
In this paper we will refer to this name as root. Since the root name
will be incorporated into Fortran90 module names, it can only contain
valid Fortran90 characters, i.e. letters, numbers, and the underscore.
To specify a KPP model, write a root file with the following lines:

::

   #MODEL      small_strato
   #LANGUAGE   Fortran90
   #DOUBLE     ON
   #INTEGRATOR rosenbrock
   #DRIVER     general
   #JACOBIAN   SPARSE_LU_ROW
   #HESSIAN    ON
   #STOICMAT   ON

The target language Fortran90 (i.e. the language of the code generated
by KPP) is selected with the command:

::

   #LANGUAGE Fortran90

Here, we have chosen Fortran90. See
Sect. `3.2.11 <#sec:command-language>`__ for other options.

The data type of the generated model can be switched between
single/double precision with the command . The command selects a
specific numerical integration routine (from the templates provided by
KPP or implemented by the user) and the command selects a specific main
program. The command selects a specific kinetic mechanism. In our
example the model definition file includes the species and the equation
files,

::

   #INCLUDE small_strato.spc
   #INCLUDE small_strato.eqn

The species file lists all the species in the model. Some of them are
variable (defined with ), meaning that their concentrations change
according to the law of mass action kinetics. Others are fixed (defined
with ), with the concentrations determined by physical and not chemical
factors. For each species its atomic composition is given (unless the
user chooses to ignore it). The atom file lists the periodic table of
elements in an section. The equation file contains the description of
the equations in an section.

::

   #INCLUDE atoms
   #DEFVAR
     O   = O;
     O1D = O;
     O3  = O + O + O;
     NO  = N + O;
     NO2 = N + O + O;
   #DEFFIX
     M   = IGNORE;
     O2  = O + O;

The chemical kinetic mechanism is specified in the KPP language in the
file . Each reaction is described as “the sum of reactants equals the
sum of products” and is followed by its rate coefficient. is the
normalized sunlight intensity, equal to one at noon and zero at night.

::

   #EQUATIONS { Stratospheric Mechanism }
   <R1>  O2  + hv = 2O       : 2.643E-10*SUN;
   <R2>  O   + O2 = O3       : 8.018E-17;
   <R3>  O3  + hv = O   + O2 : 6.120E-04*SUN;
   <R4>  O   + O3 = 2O2      : 1.576E-15;
   <R5>  O3  + hv = O1D + O2 : 1.070E-03*SUN;
   <R6>  O1D + M  = O   + M  : 7.110E-11;
   <R7>  O1D + O3 = 2O2      : 1.200E-10;
   <R8>  NO  + O3 = NO2 + O2 : 6.062E-15;
   <R9>  NO2 + O  = NO  + O2 : 1.069E-11;
   <R10> NO2 + hv = NO  + O  : 1.289E-02*SUN;

To run the model, type:

::

   kpp small_strato.kpp

Next, compile and run the Fortran90 code:

::

   make -fMakefile_small_strato
   ./small_strato.exe

.. _`sec:input`:

Input for KPP
=============

KPP basically handles two types of files: Kinetic description files and
auxiliary files. Kinetic description files are in KPP syntax and
described in the following sections. Auxiliary files are described in
Sect. `3.4 <#sec:substitution-preproc>`__. KPP kinetic description files
specify the chemical equations, the initial values of each of the
species involved, the integration parameters, and many other options.
The KPP preprocessor parses the kinetic description files and generates
several output files. Files that are written in KPP syntax have one of
the suffixes , , , or . An exception is the file , which has no suffix.

The following general rules define the structure of a kinetic
description file:

-  A KPP program is composed of KPP sections, KPP commands and inlined
   code. Their syntax is presented in the appendix.

-  Comments are either enclosed between the curly braces “``{``” and
   “``}``”, or written in a line starting with two slashes .

-  Any name given by the user to denote an atom or a species is
   restricted to be less than 32 character in length and can only
   contain letters, numbers, or the underscore character. The first
   character cannot be a number. All names are case insensitive.

The kinetic description files contain a detailed specification of the
chemical model, information about the integration method and the desired
type of results. KPP accepts only one of these files as input, but using
the command, code from separate files cn be combined. The include files
can be nested up to 10 levels. KPP will parse these files as if they
were a single big file. By carefully splitting the chemical description,
KPP can be configured for a broad range of users. In this way the users
can have direct access to that part of the model that they are
interested in, and all the other details can be hidden inside several
include files. Often, a file with atom definitions is included first,
then species definitions, and finally the equations of the chemical
mechanism.

KPP sections
------------

.. container:: center

   .. container::
      :name: tab:sections

      .. table:: KPP sections

         ==== =======================================
         name see Sect.                               
         \    `3.1.1 <#sec:section-atoms>`__          
         \    `3.1.2 <#sec:section-check>`__          
         \    `3.1.3 <#sec:section-defvar-deffix>`__  
         \    `3.1.3 <#sec:section-defvar-deffix>`__  
         \    `3.1.4 <#sec:section-equations>`__      
         \    `3.1.5 <#sec:section-initvalues>`__     
         \    `3.1.6 <#sec:section-lookat-monitor>`__ 
         \    `3.1.7 <#sec:section-lump>`__           
         \    `3.1.6 <#sec:section-lookat-monitor>`__ 
         \    `3.1.8 <#sec:section-setvar-setfix>`__  
         \    `3.1.8 <#sec:section-setvar-setfix>`__  
         \    `3.1.9 <#sec:section-transport>`__      
         \                                            
         ==== =======================================

A sign at the beginning of a line followed by a section name starts a
new KPP section. Then a list of items separated by semicolons follows. A
section ends when another KPP section or command occurs, i.e. when
another sign occurs at the beginning of a line. The syntax of an item
definition is different for each particular section.
Table `1 <#tab:sections>`__ shows all the sections defined in the KPP
language. Each of them will be described separately in the following
subsections.

.. _`sec:section-atoms`:

Atom definitions (``#ATOMS``)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The atoms that will be further used to specify the components of a
species must be declared in an section, e.g.:

::

   #ATOMS N; O; Na; Br;

Usually, the names of the atoms are the ones specified in the periodic
table of elements. For this table there is a predefined file containing
all definitions that can be used by the command:

::

   #INCLUDE atoms

This should be the first line in a KPP input file, because it allows to
use any atom in the periodic table of elements throughout the kinetic
description file.

.. _`sec:section-check`:

Mass balance checking (``#CHECK``)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

KPP is able to do a mass balance checking for all equations. Some
chemical equations are not balanced for all atoms, and this might still
be correct from a chemical point of view. To accommodate for this, KPP
can perform mass balance checking only for the list of atoms specified
in the section, e.g.:

::

   #CHECK N; C; O;

The balance checking for all atoms can be enabled by using the command.
Without or , no checking is performed. The atom can also be used to
control mass balance checking.

.. _`sec:section-defvar-deffix`:

Species definitions (``#DEFVAR`` and ``#DEFFIX``)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

There are two ways to declare new species together with their atom
composition: and . These sections define all the species that will be
used in the chemical mechanism. Species can be variable or fixed. The
type is implicitly specified by defining the species in the appropriate
sections. A species can be considered fixed if its concentration does
not vary too much. The variable species are medium or short lived
species and their concentrations vary in time. This division of species
into different categories is helpful for integrators that benefit from
treating them differently.

For each species the user has to declare the atom composition. This
information is used for mass balance checking. If the species is a
lumped species without an exact composition, it can be ignored. To do
this one can declare the predefined atom as being part of the species
composition. Examples for these sections are:

::

   #DEFVAR
     NO2 = N + 2O;
     CH3OOH = C + 4H + 2O;
     HSO4m = IGNORE;
     RCHO = IGNORE;
   #DEFFIX
     CO2 = C + 2O;

.. _`sec:section-equations`:

Equations (``#EQUATIONS``)
~~~~~~~~~~~~~~~~~~~~~~~~~~

The chemical mechanism is specified in the section. Each equation is
written in the natural way in which a chemist would write it, e.g.:

::

   #EQUATIONS
     NO2 + hv = NO + O : 0.533*SUN;
     OH + NO2 = HNO3 : k_3rd(temp,
       cair,2.E-30,3.,2.5E-11,0.,0.6);

Only the names of already defined species can be used. The rate
coefficient has to be placed at the end of each equation, separated by a
colon. The rate coefficient does not necessarily need to be a numerical
value. Instead, it can be a valid expression in the target language. If
there are several sections in the input, their contents will be
concatenated.

A minus sign in an equation shows that a species is consumed in a
reaction but it does not affect the reaction rate. For example, the
oxidation of methane can be written as:

::

   CH4 + OH = CH3OO + H2O - O2 : k_CH4_OH;

However, it should be noted that using negative products may lead to
numerical instabilities.

Often, the stoichiometric factors are integers. However, it is also
possible to have non-integer yields, which is very useful to
parameterize organic reactions that branch into several side reactions:

::

   CH4 + O1D = .75 CH3O2 + .75 OH + .25 HCHO
               + .4 H + .05 H2 : k_CH4_O1D;

**((TODO: CHECK IF THE FOLLOWING DESCRIPTION OF PROD AND HV IS
CORRECT))**

KPP provides two pre-defined dummy species: and . Using dummy species
does not affect the numerics of the integrators. It only serves to
improve the readability of the equations. For photolysis reactions, can
be specified as one of the reagents to indicate that light
(:math:`h\nu`) is needed for this reaction, e.g.:

::

     NO2 + hv = NO + O : J_NO2;

When the products of a reaction are not known oder not important, the
dummy species should be used as a product. This is necessary because the
KPP syntax does not allow an empty list of products. For example, the
dry deposition of atmospheric ozone to the surface can be written as:

::

   O3 = PROD : v_d_O3;

The same equation must not occur twice in the section. For example, you
may have both the gas-phase reaction of with water in your mechanism and
also the heterogeneous reaction on aerosols:

::

   N2O5 + H2O = 2 HNO3 : k_gas;
   N2O5 + H2O = 2 HNO3 : k_aerosol;

These reactions must be merged by adding the rate coefficients:

::

   N2O5 + H2O = 2 HNO3 : k_gas+k_aerosol;

.. _`sec:section-initvalues`:

Initial values (``#INITVALUES``)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The initial concentration values for all species can be defined in the
section, e.g.:

::

   #INITVALUES
     CFACTOR = 2.5E19;
     NO2 = 1.4E-9;
     CO2 = MyCO2Func();
     ALL_SPEC = 0.0;

If no value is specified for a particular species, the default value
zero is used. One can set the default values using the generic species
names: , , and . In order to use coherent units for concentration and
rate coefficients, it is sometimes necessary to multiply each value by a
constant factor. This factor can be set by using the generic name . Each
of the initial values will be multiplied by this factor before being
used. If is omitted, it defaults to one.

The information gathered in this section is used to generate the
subroutine (see Sect. `4.1.3 <#sec:output-init>`__). In more complex 3D
models, the initial values are usually taken from some input files or
some global data structures. In this case, may not be needed.

.. _`sec:section-lookat-monitor`:

Output data selection (``#LOOKAT`` and ``#MONITOR``)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

There are two sections in this category: and .

The section instructs the preprocessor what are the species for which
the evolution of the concentration, should be saved in a data file. By
default, if no section is present, all the species are saved. If an atom
is specified in the list then the total mass of the particular atom is
reported. This allows to check how the mass of a specific atom was
conserved by the integration method. The command can be used to specify
all the species. Output of can be directed to the file root using the
utility subroutines described in Sect. `4.1.16 <#sec:output-utility>`__.

The section defines a different list of species and atoms. This list is
used by the driver to display the concentration of the elements in the
list during the integration. This may give us a feedback of the
evolution in time of the selected species during the integration. The
syntax is similar to the section. With the driver , output of goes to
the screen (STDOUT). The order of the output is: first variable species,
then fixes species, finally atoms. It is not the order in the command.

Examples for these sections are:

::

   #LOOKAT NO2; CO2; O3; N;
   #MONITOR O3; N;

.. _`sec:section-lump`:

Lump species definitions (``#LUMP``)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

To reduce the stiffness of some models, various lumping of species may
be defined in the section. The example below shows how the concentration
of can be replaced by the sum of concentrations for and which is
considered to be a single variable. At the end of the integration, the
concentration of is computed by substraction from the lumped variable.

::

   #LUMP NO2 + NO : NO2

.. _`sec:section-setvar-setfix`:

Redefining species definitions (``#SETVAR`` and ``#SETFIX``)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The commands and change the type of an already defined species. Then,
depending on the integration method, one may or may not use the initial
classification, or can easily move one species from one category to
another. The use of the generic species , , and is also allowed.
Examples for these sections are:

::

   #SETVAR ALL_SPEC;
   #SETFIX H2O; CO2;

.. _`sec:section-transport`:

Transport (``#TRANSPORT``)
~~~~~~~~~~~~~~~~~~~~~~~~~~

The section is only used for transport chemistry models. It specifies
the list of species that needs to be included in the transport model,
e.g.:

::

   #TRANSPORT NO2; CO2; O3; N;

One may use a more complex chemical model from which only a couple of
species are considered for the transport calculations. The command is
also available as a shorthand for specifying that all the species used
in the chemical model have to be included in the transport calculations.

KPP commands
------------

.. container:: center

   .. container::
      :name: tab:commands

      .. table:: KPP commands

         ==== ===========================================
         name see Sect.                                   
         \    `3.2.18 <#sec:command-shorthand>`__         
         \    `3.2.1 <#sec:command-declare>`__            
         \    `3.2.2 <#sec:command-double>`__             
         \    `3.2.3 <#sec:command-driver>`__             
         \    `3.2.4 <#sec:command-dummyindex>`__         
         \    `3.2.5 <#sec:command-eqntags>`__            
         \    `3.2.6 <#sec:command-function>`__           
         \    `3.2.7 <#sec:command-hessian>`__            
         \    `3.2.8 <#sec:command-include>`__            
         \    `3.2.9 <#sec:command-integrator-intfile>`__ 
         \    `3.2.9 <#sec:command-integrator-intfile>`__ 
         \    `3.2.10 <#sec:command-jacobian>`__          
         \    `3.2.11 <#sec:command-language>`__          
         \    `3.2.18 <#sec:command-shorthand>`__         
         \    `3.2.12 <#sec:command-mex>`__               
         \    `3.2.13 <#sec:command-minversion>`__        
         \    `3.2.14 <#sec:command-model>`__             
         \    `3.2.15 <#sec:command-reorder>`__           
         \    `3.2.16 <#sec:command-stochastic>`__        
         \    `3.2.17 <#sec:command-stoicmat>`__          
         \    `3.2.18 <#sec:command-shorthand>`__         
         \    `3.2.19 <#sec:command-uppercasef90>`__      
         \                                                
         ==== ===========================================

A command begins on a new line with a sign, followed by a command name
and one or more parameters. A summary of the commands available in KPP
is shown in Table `2 <#tab:commands>`__. Details about each command are
given in the following subsections.

.. _`sec:command-declare`:

Referencing constants (``#DECLARE``)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**((TODO: CHECK IF THIS SECTION IS OKAY))**

The command determines how constants like ``dp``, ``NSPEC``, ``NVAR``,
``NFIX``, and ``NREACT`` are inserted into the KPP-generated code. (the
default) means that they are referenced by their names
(e.g. “C(NSPEC)”), whereas means that their values are inserted
(e.g. “C(7)”).

.. _`sec:command-double`:

Precision control (``#DOUBLE``)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The command selects single or double precision arithmetic. (the default)
means use double precision, means use single precision (see
Sect. `4.1.6 <#sec:output-precision>`__).

.. _`sec:command-driver`:

Driver selection (``#DRIVER``)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The command selects the driver, i.e., the file from which the main
function is to be taken. The parameter is a file name, without suffix.
The appropriate suffix (, , , or ) is automatically appended.

Normally, KPP tries to find the selected driver file in the directory
``$KPP_HOME/drv/``. However, if the supplied file name contains a slash,
it is assumed to be absolute. To access a driver in the current
directory, the prefix can be used, e.g.:

::

   #DRIVER ./mydriver

It is possible to choose the empty dummy driver , if the user wants to
include the KPP generated modules into a larger model (e.g. a general
circulation or a chemical transport model) instead of creating a
stand-alone version of the chemical integrator. The driver is also
selected when the command is missing. If the command occurs twice, the
second replaces the first.

.. _`sec:command-dummyindex`:

Dummy indices (``#DUMMYINDEX``)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

It is possible to declare species in the section that are not used in
the section. If your model needs to check at run-time if a certain
species is included in the current mechanism, you can set to . Then, KPP
will set the indices to zero for all species that do not occur in any
reaction. With (the default), those are undefined variables. For
example, if you frequently switch between mechanisms with and without
sulfuric acid, you can use this code:

::

   IF (ind_H2SO4=0) THEN
     PRINT *, 'no H2SO4 in current mechanism'
   ELSE
     PRINT *, 'c(H2SO4) =', C(ind_H2SO4)
   ENDIF

.. _`sec:command-eqntags`:

Generation of equation tags (``#EQNTAGS``)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Each reaction in the section may start with an equation tag which is
enclosed in angle brackets, e.g.:

::

   <J1> NO2 + hv = NO + O : 0.533*SUN;

With set to , this equation tag can be used to refer to a specific
equation, as described in Sect. `4.1.5 <#sec:output-monitor>`__. The
default for is .

.. _`sec:command-function`:

The function generation (``#FUNCTION``)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The command controls which functions are generated to compute the
production/destruction terms for variable species. generates one
function that computes the normal derivatives. generates two functions
for the derivatives in production and destruction forms.

.. _`sec:command-hessian`:

Generation of Hessian (``#HESSIAN``)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The option (the default) of the command turns the Hessian generation on
(see Sect. `4.1.12 <#sec:output-ode-hess>`__). With it is switched off.

.. _`sec:command-include`:

File include command (``#INCLUDE``)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The command instructs KPP to look for the file specified as a parameter
and parse the content of this file before proceeding to the next line.
This allows the atoms definition, the species definition and even the
equation definition to be shared between several models. Moreover this
allows for custom configuration of KPP to accommodate various classes of
users. Include files can be either in one of the KPP directories or in
the current directory.

.. _`sec:command-integrator-intfile`:

Integrator selection (``#INTEGRATOR`` and ``#INTFILE``)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The command selects the integrator definition file. The parameter is the
file name of an integrator, without suffix. The effect of:

*integrator*

is similar to:

``$KPP_HOME/int/``\ *integrator*

The integrator definition file selects an integrator file with and also
defines some suitable options for it. The command selects the file that
contains the integrator routine. This command allows the use of
different integration techniques on the same model. The parameter of the
command is a file name, without suffix. The appropriate suffix (, , , or
) is appended and the result selects the file from which the integrator
is taken. This file will be copied into the code file in the appropriate
place. All integrators have to conform to the same specific calling
sequence. Normally, KPP tries to find the selected integrator file in
the directory ``$KPP_HOME/int/``. However, if the supplied file name
contains a slash, it is assumed to be absolute. To access an integrator
in the current directory, the prefix can be used, e.g.:

::

   #INTEGRATOR ./mydeffile
   #INTFILE ./myintegrator

If the command occurs twice, the second replaces the first.

.. _`sec:command-jacobian`:

The Jacobian (``#JACOBIAN``)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The command controls which functions are generated to compute the
Jacobian. The option inhibits the generation of the Jacobian routine.
The option generates the Jacobian as a square (:math:`\times`) matrix.
It should be used if the integrator needs the whole Jacobians. The
options and (the default) both generate the Jacobian in sparse
(compressed on rows) format. They should be used if the integrator needs
the whole Jacobian, but in a sparse form. The format used is compressed
on rows. With , KPP extends the number of nonzeros to account for the
fill-in due to the LU decomposition.

.. _`sec:command-language`:

Target language selection (``#LANGUAGE``)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The command selects the target language in which the code file is to be
generated. Available options are , , , or .

.. _`sec:command-mex`:

Mex files (``#MEX``)
~~~~~~~~~~~~~~~~~~~~

Mex is a *m*\ atlab *ex*\ tension that allows to call functions written
in Fortran and C directly from within the Matlab environment. KPP
generates the mex interface routines for the ODE function, Jacobian, and
Hessian, for the target languages C, Fortran77, and Fortran90. The
default is . With , no Mex files are generated.

.. _`sec:command-minversion`:

Specify a minimum KPP version to use (``#MINVERSION``)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

You may restrict a chemical mechanism to use a given version of KPP or
later. To do this, add . The version number adheres to the Semantic
Versioning style (https://semver.org), where is the major version
number, is the minor version number, and is the bugfix (aka “patch”)
version number.

For example, if is specified, then KPP will quit with an error message
unless you are using KPP 2.4.0 or later.

.. _`sec:command-model`:

Selcting a chemical model (``#MODEL``)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The chemical model contains the description of the atoms, species, and
chemical equations. It also contains default initial values for the
species and default options including the best integrator for the model.
In the simplest case, the main kinetic description file, i.e. the one
passed as parameter to KPP, can contain just a single line selecting the
model. KPP tries to find a file with the name of the model and the
suffix in the ``$KPP_HOME/models`` subdirectory. This file is then
parsed. The content of the model definition file is written in the KPP
language. The model definition file points to a species file and an
equation file. The species file includes further the atom definition
file. All default values regarding the model are automatically selected.
For convenience, the best integrator and driver for the given model are
also automatically selected.

The command is optional, and intended for using a predefined model.
Users who supply their own reaction mechanism do not need it.

.. _`sec:command-reorder`:

Reordering (``#REORDER``)
~~~~~~~~~~~~~~~~~~~~~~~~~

Reordering of the species is performed in order to minimize the fill-in
during the LU factorization, and therefore preserve the sparsity
structure and increase efficiency. The reordering is done using a
diagonal markowitz algorithm. The details are explained in
:raw-latex:`\citet{IMPLEMENTATION}`. The default is . means that KPP
does not reorder the species. The order of the variables is the order in
which the species are declared in the section.

.. _`sec:command-stochastic`:

Stochastic simulation (``#STOCHASTIC``)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The option of the command turns on the generation of code for stochastic
kinetic simulations (see Sect. `4.1.15 <#sec:output-stochastic>`__). The
default option is .

.. _`sec:command-stoicmat`:

The Stoichiometric Formulation (``#STOICMAT``)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Unless this command is set to , KPP generates code for the
stoichiometric matrix, the vector of reactant products in each reaction,
and the partial derivative of the time derivative function with respect
to rate coefficients. These elements are discussed in
Sect. `4.1.14 <#sec:output-stoichiom>`__.

.. _`sec:command-shorthand`:

Shorthand commands (``#CHECKALL``, ``#LOOKATALL`` and ``#TRANSPORTALL``)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

KPP defines a couple of shorthand commands. The commands that fall into
this category are , and . All of them have been described in the
previous sections.

.. _`sec:command-uppercasef90`:

Generate Fortran90 files with the ``.F90`` suffix (``#UPPERCASEF90``)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

If you have selected the option, KPP will generate source code ending in
by default. Setting ``#UPPERCASEF90 ON`` will tell KPP to generate
Fortran90 code ending in instead.

Inlined code
------------

In order to offer maximum flexibility, KPP allows the user to include
pieces of code in the kinetic description file. Inlined code begins on a
new line with and the *inline_type*. Next, one or more lines of code
follow, written in the target language (Fortran90, Fortran77, C, or
Matlab) as specified by the *inline_type*. The inlined code ends with .
The code is inserted into the KPP output at a position which is also
determined by *inline_type* as explained in
Table `[tab:inlining] <#tab:inlining>`__. If two inline commands with
the same inline type are declared, then the contents of the second is
appended to the first one. In this manual, we show the inline types for
Fortran90. The inline types for the other languages are produced by
replacing by , , or , respectively.

.. container:: table*

   .. container:: center

      +---------------+------+---------------------+---------------------+
      | *inline_type* | file | placement           | usage               |
      +---------------+------+---------------------+---------------------+
      |               | root | specification       | (obsolete)          |
      |               |      | section             |                     |
      +---------------+------+---------------------+---------------------+
      |               | root | specification       | global variables    |
      |               |      | section             |                     |
      +---------------+------+---------------------+---------------------+
      |               | root | subroutine          | integration         |
      |               |      |                     | parameters          |
      +---------------+------+---------------------+---------------------+
      |               | root | executable section  | rate law functions  |
      +---------------+------+---------------------+---------------------+
      |               | root | subroutine          | statements and      |
      |               |      |                     | definitions of rate |
      |               |      |                     | coefficients        |
      +---------------+------+---------------------+---------------------+
      |               | root | executable section  | utility functions   |
      +---------------+------+---------------------+---------------------+
      |               |      |                     |                     |
      +---------------+------+---------------------+---------------------+

.. _`sec:f90-data`:

Inline type ``F90_DATA``
~~~~~~~~~~~~~~~~~~~~~~~~

This inline type was introduced in a previous version of KPP to
initialize variables. It is now obsolete but kept for compatibility. For
Fortran90, should be used instead.

.. _`sec:f90-global`:

Inline type ``F90_GLOBAL``
~~~~~~~~~~~~~~~~~~~~~~~~~~

can be used to declare global variables, e.g. for a special rate
coefficient:

::

   #INLINE F90_GLOBAL
     REAL(dp) :: k_DMS_OH
   #ENDINLINE

.. _`sec:f90-init`:

Inline type ``F90_INIT``
~~~~~~~~~~~~~~~~~~~~~~~~

can be used to define initial values before the start of the
integartion, e.g.:

::

   #INLINE F90_INIT
     TSTART = (12.*3600.)
     TEND = TSTART + (3.*24.*3600.)
     DT = 0.25*3600.
     TEMP = 270.
   #ENDINLINE

.. _`sec:f90-rates`:

Inline type ``F90_RATES``
~~~~~~~~~~~~~~~~~~~~~~~~~

can be used to add new subroutines to calculate rate coefficients, e.g.:

::

   #INLINE F90_RATES
     REAL FUNCTION k_SIV_H2O2(k_298,tdep,cHp,temp)
       ! special rate function for S(IV) + H2O2
       REAL, INTENT(IN) :: k_298, tdep, cHp, temp
       k_SIV_H2O2 = k_298 &
         * EXP(tdep*(1./temp-3.3540E-3)) &
         * cHp / (cHp+0.1)
     END FUNCTION k_SIV_H2O2
   #ENDINLINE

.. _`sec:f90-rconst`:

Inline type ``F90_RCONST``
~~~~~~~~~~~~~~~~~~~~~~~~~~

can be used to define time-dependent values of rate coefficients that
were declared with :

::

   #INLINE F90_RCONST
     k_DMS_OH = 1.E-9*EXP(5820./temp)*C(ind_O2)/ &
       (1.E30+5.*EXP(6280./temp)*C(ind_O2))
   #ENDINLINE

.. _`sec:f90-util`:

Inline type ``F90_UTIL``
~~~~~~~~~~~~~~~~~~~~~~~~

can be used to define utility subroutines.

.. _`sec:substitution-preproc`:

Auxiliary files and the substitution preprocessor
-------------------------------------------------

.. container:: table*

   .. container:: center

      ==== ==========================================
      File Contents                                   
      \    derivatives with respect to reaction rates 
      \    derivatives with respect to reaction rates 
      \    unix makefiles                             
      \    mex files                                  
      \    mex files                                  
      \    mex files                                  
      \    Sparse utility functions                   
      \    Function related to equation tags          
      \    Function related to solar zenith angle     
      \    User-defined rate law functions            
      \    Input/output utilities                     
      \                                               
      ==== ==========================================

.. container:: table*

   .. container:: center

      =========== ==================================================== =======
      Placeholder Replaced by                                          Example
      \           the root name                                        
      \           the real data type                                   
      \           number of species                                    7
      \           number of variable species                           5
      \           number of fixed species                              2
      \           number of chemical reactions                         10
      \           number of Jacobian nonzero elements                  18
      \           number of Jacobian nonzero elements, with LU fill-in 19
      \           number of Hessian nonzero elements                   10
      \                                                                
      =========== ==================================================== =======

The auxiliary files (listed in Table `[tab:aux] <#tab:aux>`__) are
templates for integrators, drivers, and utilities. They are inserted
into the KPP output after being run through the substitution
preprocessor. This preprocessor replaces several placeholders (listed in
Table `[tab:substitution] <#tab:substitution>`__) in the template files
with their particular values in the model at hand. Usually, only and are
needed because the other values can also be obtained via the variables
listed in Tab. `3 <#tab:parameters>`__.

is replaced by the appropriate single or double precision declaration
type. Depending on the target language KPP will select the correct
declaration type. For example if one needs to declare an array BIG of
size 1000, a declaration like the following must be used:

When used with the option , the above line will be automatically
translated into:

and when used with the option , the same line will become:

in the resulting Fortran90 output file.

is replaced by the root file name of the main kinetic description file.
In our example where we are processing , a line in an auxiliary
Fortran90 file like

::

   USE KPP_ROOT_Monitor

will be translated into

::

   USE small_strato_Monitor

in the generated Fortran90 output file.

.. _`sec:output`:

Output from KPP
===============

The Fortran90 Code
------------------

.. container:: table*

   .. container:: center

      +------+---------------------+---------------+---------------------+
      | File | Description         | Only if…      | see Sect.           |
      +------+---------------------+---------------+---------------------+
      | root | Driver              | :math:`\neq`  | `4.1.1 <#           |
      |      |                     |               | sec:output-main>`__ |
      +------+---------------------+---------------+---------------------+
      | root | ODE function        |               | `4.1.10 <#sec       |
      |      |                     |               | :output-ode-fun>`__ |
      +------+---------------------+---------------+---------------------+
      | root | Global data headers |               | `4.1.9 <#se         |
      |      |                     |               | c:output-global>`__ |
      +------+---------------------+---------------+---------------------+
      | root | Initialization      |               | `4.1.3 <#           |
      |      |                     |               | sec:output-init>`__ |
      +------+---------------------+---------------+---------------------+
      | root | Numerical           |               | `4.1.4 <#sec:ou     |
      |      | integration         |               | tput-integrator>`__ |
      +------+---------------------+---------------+---------------------+
      | root | Sparse linear       |               | `4.1.13             |
      |      | algebra             |               | <#sec:output-la>`__ |
      +------+---------------------+---------------+---------------------+
      | root | Summary of modules  |               | `4.1.2 <#s          |
      |      |                     |               | ec:output-model>`__ |
      +------+---------------------+---------------+---------------------+
      | root | Equation info       |               | `4.1.5 <#sec        |
      |      |                     |               | :output-monitor>`__ |
      +------+---------------------+---------------+---------------------+
      | root | Model parameters    |               | `4.1.8 <#sec:ou     |
      |      |                     |               | tput-parameters>`__ |
      +------+---------------------+---------------+---------------------+
      | root | Parameterized types |               | `4.1.6 <#sec:o      |
      |      |                     |               | utput-precision>`__ |
      +------+---------------------+---------------+---------------------+
      | root | User-defined rate   |               | `4.1.7 <#s          |
      |      | laws                |               | ec:output-rates>`__ |
      +------+---------------------+---------------+---------------------+
      | root | Utility             |               | `4.1.16 <#sec       |
      |      | input-output        |               | :output-utility>`__ |
      +------+---------------------+---------------+---------------------+
      | root | ODE Jacobian        |               | `4.1.11 <#sec       |
      |      |                     |               | :output-ode-jac>`__ |
      +------+---------------------+---------------+---------------------+
      | root | Jacobian sparsity   | :math:`*`     | `4.1.11 <#sec       |
      |      |                     |               | :output-ode-jac>`__ |
      +------+---------------------+---------------+---------------------+
      | root | ODE Hessian         |               | `4.1.12 <#sec:      |
      |      |                     |               | output-ode-hess>`__ |
      +------+---------------------+---------------+---------------------+
      | root | Sparse Hessian data | and :math:`*` | `4.1.12 <#sec:      |
      |      |                     |               | output-ode-hess>`__ |
      +------+---------------------+---------------+---------------------+
      | root | Stochastic          |               | `4.1.15 <#sec:ou    |
      |      | functions           |               | tput-stochastic>`__ |
      +------+---------------------+---------------+---------------------+
      | root | Stoichiometric      |               | `4.1.14 <#sec:o     |
      |      | model               |               | utput-stoichiom>`__ |
      +------+---------------------+---------------+---------------------+
      | root | Stoichiometric      | and :math:`*` | `4.1.14 <#sec:o     |
      |      | matrix              |               | utput-stoichiom>`__ |
      +------+---------------------+---------------+---------------------+
      | root | Matlab interface    |               | `4.1.17 <#sec       |
      |      | Fun                 |               | :output-mexcode>`__ |
      +------+---------------------+---------------+---------------------+
      | root | Matlab interface    | and :math:`*` | `4.1.17 <#sec       |
      |      | Jac                 |               | :output-mexcode>`__ |
      +------+---------------------+---------------+---------------------+
      | root | Matlab interface    | and           | `4.1.17 <#sec       |
      |      | Hess                |               | :output-mexcode>`__ |
      +------+---------------------+---------------+---------------------+
      |      | Makefile            |               | `4.1.18 <#sec:      |
      |      |                     |               | output-makefile>`__ |
      +------+---------------------+---------------+---------------------+
      | root | Human-readable info |               | `4.5 <              |
      |      |                     |               | #sec:output-map>`__ |
      +------+---------------------+---------------+---------------------+
      |      |                     |               |                     |
      +------+---------------------+---------------+---------------------+

.. container:: figure*

   .. image:: Figures/kpp2_use_diagr.pdf
      :alt: image

.. container:: table*

   .. container:: center

      ========== ==================================================== ====
      Subroutine Description                                          File
      \          ODE function                                         root
      \          ODE Jacobian in sparse format                        root
      \          sparse multiplication                                root
      \          sparse multiplication                                root
      \          ODE Jacobian in full format                          root
      \          ODE Hessian in sparse format                         root
      \          Hessian action on vectors                            root
      \          Transposed Hessian action on vectors                 root
      \          Derivatives of Fun with respect to rate coefficients root
      \          Derivatives of Jac with respect to rate coefficients root
      \          Reactant products                                    root
      \          Jacobian of reactant products                        root
      \          Sparse LU decomposition                              root
      \          Sparse back substitution                             root
      \          Update photolysis rate coefficients                  root
      \          Update all rate coefficients                         root
      \          Update solar intensity                               root
      \          Set initial values                                   root
      \          Integrate one time step                              root
      \          Check mass balance for selected atoms                root
      \          Shuffle concentration vector                         root
      \          Shuffle concentration vector                         root
      \          Utility for command                                  root
      \          Utility for command                                  root
      \          Utility for command                                  root
      \          Calculate reaction number from equation tag          root
      \                                                               
      ========== ==================================================== ====

The code generated by KPP is organized in a set of separate files. Each
has a time stamp and a complete description of how it was generated at
the begining of the file. The files associated with root are named with
a corresponding prefix “root”. The list of files and a short description
is shown in Table `[tab:generated_files] <#tab:generated_files>`__. All
subroutines and functions, global parameters, variables, and sparsity
data structures are encapsulated in modules. There is exactly one module
in each file, and the name of the module is identical to the file name
but without the suffix . Fig. `[fig:use_diagr] <#fig:use_diagr>`__ shows
how these modules are related to each other. A concise list of the main
subroutines generated by KPP is shown in
Table `[tab:functions] <#tab:functions>`__. The generated code is
consistent with the Fortran90 standard. It will not exceed the maximum
number of 39 continuation lines. If KPP cannot properly split an
expression to keep the number of continuation lines below the threshold
then it will generate a warning message pointing to the location of this
expression.

.. _`sec:output-main`:

root\ ``_Main.f90``
~~~~~~~~~~~~~~~~~~~

root is the main Fortran90 program. It contains the driver after
modifications by the substitution preprocessor. The name of the file is
computed by KPP by appending the suffix to the root name.

.. _`sec:output-model`:

root\ ``_Model.f90``
~~~~~~~~~~~~~~~~~~~~

The file root\ ``_Model.f90`` completely defines the model by using all
the associated modules.

.. _`sec:output-init`:

root\ ``_Initialize.f90``
~~~~~~~~~~~~~~~~~~~~~~~~~

The file root\ ``_Initialize.f90`` contains the subroutine which defines
initial values of the chemical species. The driver calls the subroutine
once before the time integration loop starts.

.. _`sec:output-integrator`:

root\ ``_Integrator.f90``
~~~~~~~~~~~~~~~~~~~~~~~~~

The file root\ ``_Integrator.f90`` contains the subroutine which is
called every time step during the integration. The integrator that was
chosen with is also included in root\ ``_Integrator.f90``. In case of an
unsuccessful integration, the module root provides a short error message
in the public variable .

.. _`sec:output-monitor`:

root\ ``_Monitor.f90``
~~~~~~~~~~~~~~~~~~~~~~

The file root\ ``_Monitor.f90`` contains arrays with information about
the chemical mechanism. The names of all species are included in and the
names of all equations are included in .

It was shown above (Sect. `3.2.5 <#sec:command-eqntags>`__) that each
reaction in the section may start with an equation tag which is enclosed
in angle brackets, e.g.:

::

   <J1> NO2 + hv = NO + O : 0.533*SUN;

If the equation tags are switched on, KPP also generates the array . In
combination with and the function that converts the equation tag to the
KPP-internal equation number, this can be used to describe a reaction:

::

     PRINT *,'Reaction J1 is:', &
       EQN_NAMES(tag2num('J1'))

.. _`sec:output-precision`:

root\ ``_Precision.f90``
~~~~~~~~~~~~~~~~~~~~~~~~

Fortran90 code uses parameterized real types. root contains the
following real kind definitions:

::

   ! KPP_SP - Single precision kind
     INTEGER, PARAMETER :: &
       SP = SELECTED_REAL_KIND(6,30)
   ! KPP_DP - Double precision kind
     INTEGER, PARAMETER :: &
       DP = SELECTED_REAL_KIND(12,300)

Depending on the choice of the command, the real variables are of type
double () or single precision (). Changing the parameters of the
function in this module will cause a change in the working precision for
the whole model.

.. _`sec:output-rates`:

root\ ``_Rates.f90``
~~~~~~~~~~~~~~~~~~~~

The code to update the rate constants is in root. The user defined rate
law functions are also placed here.

.. _`sec:output-parameters`:

root\ ``_Parameters.f90``
~~~~~~~~~~~~~~~~~~~~~~~~~

The global parameters (Table `3 <#tab:parameters>`__) are defined and
initialized in root.

KPP orders the variable species such that the sparsity pattern of the
Jacobian is maintained after an LU decomposition. For our example there
are five variable species (=5) ordered as

::

   ind_O1D=1, ind_O=2, ind_O3=3,
   ind_NO=4, ind_NO2=5

and two fixed species (=2)

::

   ind_M = 6, ind_O2 = 7.

KPP defines a complete set of simulation parameters, including the
numbers of variable and fixed species, the number of chemical reactions,
the number of nonzero entries in the sparse Jacobian and in the sparse
Hessian, etc. Some important simulation parameters generated by KPP are
presented in Table `3 <#tab:parameters>`__.

.. container::
   :name: tab:parameters

   .. table:: List of important simulation parameters and their values
   for the ``small_strato`` example

      ========= ================================ =====
      Parameter Represents                       Value
      \         No. chemical species             7
      \         No. variable species             5
      \         No. fixed species                2
      \         No. reactions                    10
      \         No. nonzero entries Jacobian     18
      \         As above, after LU factorization 19
      \         Length, sparse Hessian           10
      \         Length, sparse Jacobian          13
      \         Length, stoichiometric matrix    22
      \         Index of species *spc* in        
      \         Index of fixed species *spc* in  
      \                                          
      ========= ================================ =====

.. _`sec:output-global`:

root\ ``_Global.f90``
~~~~~~~~~~~~~~~~~~~~~

The global variables (Table `4 <#tab:global>`__) are declared in root.
Global variables are presented in Table `4 <#tab:global>`__.

.. container::
   :name: tab:global

   .. table:: List of important global variables

      =============== ================================
      Global variable Represents
      \               Concentrations, all species
      \               Concentrations, variable species
      \               Concentrations, fixed species
      \               Rate coefficient values
      \               Current integration time
      \               Sun intensity between 0 and 1
      \               Temperature
      \               Simulation start/end time
      \               Simulation step
      \               Absolute tolerances
      \               Relative tolerances
      \               Lower bound for time step
      \               Upper bound for time step
      \               Conversion factor
      \               Names of chemical species
      \               Names of chemical equations
      \               
      =============== ================================

Both variable and fixed species are stored in the one-dimensional array
. The first part (indices from 1 to ) contains the variable species, and
the second part (indices from to ) the fixed species. The total number
of species is the sum of the and . The parts can also be accessed
separately through the arrays and :

::

   VAR(1:NVAR) = C(1:NVAR)
   FIX(1:NFIX) = C(NVAR+1:NSPEC)

.. _`sec:output-ode-fun`:

root\ ``_Function.f90``
~~~~~~~~~~~~~~~~~~~~~~~

The chemical ODE system for our example is:

.. math::

   \begin{aligned}
     \frac{{\mathrm{d}}\, [\chem{O(^1D)}]}{{\mathrm{d}}t} & = & k_{5}\, [\chem{O_3}] - k_{6}\,
     [\chem{O(^1D)}]\, [\chem{M}] - k_{7}\, [\chem{O(^1D)}]\, [\chem{O_3}]\\
     \frac{{\mathrm{d}}\, [\chem{O}]}{{\mathrm{d}}t} & = & 2\, k_{1}\, [\chem{O_2}] - k_{2}\,
     [\chem{O}]\, [\chem{O_2}] + k_{3}\, [\chem{O_3}]\\
     & & - k_{4}\, [\chem{O}]\, [\chem{O_3}]+ k_{6}\, [\chem{O(^1D)}]\,
     [\chem{M}]\\
     & & - k_{9}\, [\chem{O}]\, [\chem{NO_2}] + k_{10}\, [\chem{NO_2}]\\
     \frac{{\mathrm{d}}\, [\chem{O_3}]}{{\mathrm{d}}t} & = & k_{2}\, [\chem{O}]\, [\chem{O_2}] -
     k_{3}\, [\chem{O_3}] - k_{4}\, [\chem{O}]\, [\chem{O_3}] - k_{5}\,
     [\chem{O_3}]\\
     & & - k_{7}\, [\chem{O(^1D)}]\, [\chem{O_3}] - k_{8}\, [\chem{O_3}]\, 
     [\chem{NO}]\\
     \frac{{\mathrm{d}}\, [\chem{NO}]}{{\mathrm{d}}t} & = & - k_{8}\, [\chem{O_3}]\, [\chem{NO}] +
     k_{9}\, [\chem{O}]\, [\chem{NO_2}] + k_{10}\, [\chem{NO_2}]\\
     \frac{{\mathrm{d}}\, [\chem{NO_2}]}{{\mathrm{d}}t} & = & k_{8}\, [\chem{O_3}]\, [\chem{NO}] -
     k_{9}\, [\chem{O}]\, [\chem{NO_2}] - k_{10}\, [\chem{NO_2}]\end{aligned}

where square brackets denote concentrations of the species. The code for
the ODE function is in root. The chemical reaction mechanism represents
a set of ordinary differential equations (ODEs) of dimension . The
concentrations of fixed species are parameters in the derivative
function. The subroutine computes first the vector of reaction rates and
then the vector of variable species time derivatives. The input
arguments , , and are the concentrations of variable species, fixed
species, and the rate coefficients, respectively. Below is the Fortran90
code generated by KPP for the ODE function of our example.

::

   SUBROUTINE Fun (V, F, RCT, Vdot )
      REAL(kind=DP) ::  V(NVAR), &
            F(NFIX), RCT(NREACT), &
            Vdot(NVAR), A(NREACT) &
   ! Computation of equation rates
      A(1) = RCT(1)*F(2)
      A(2) = RCT(2)*V(2)*F(2)
      A(3) = RCT(3)*V(3)
      A(4) = RCT(4)*V(2)*V(3)
      A(5) = RCT(5)*V(3)
      A(6) = RCT(6)*V(1)*F(1)
      A(7) = RCT(7)*V(1)*V(3)
      A(8) = RCT(8)*V(3)*V(4)
      A(9) = RCT(9)*V(2)*V(5)
      A(10) = RCT(10)*V(5)
   ! Aggregate function
      Vdot(1) = A(5)-A(6)-A(7)
      Vdot(2) = 2*A(1)-A(2)+A(3) &
                -A(4)+A(6)-A(9)+A(10)
      Vdot(3) = A(2)-A(3)-A(4)-A(5) &
                -A(7)-A(8)
      Vdot(4) = -A(8)+A(9)+A(10)
      Vdot(5) = A(8)-A(9)-A(10)
   END SUBROUTINE Fun

.. _`sec:output-ode-jac`:

root\ ``_Jacobian.f90`` and root\ ``_JacobianSP.f90``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. figure:: Figures/small_jac.pdf
   :alt: The sparsity pattern of the Jacobian for the ``small_strato``
   example. All non-zero elements are marked with a bullet. Note that
   even though ``J(3,5)=0``, it is also included here because of the
   fill-in. **((TODO: ADRIAN, IS THIS EXPLANATION OF J(3,5) CORRECT?))**
   :name: fig:jac

   The sparsity pattern of the Jacobian for the ``small_strato``
   example. All non-zero elements are marked with a bullet. Note that
   even though ``J(3,5)=0``, it is also included here because of the
   fill-in. **((TODO: ADRIAN, IS THIS EXPLANATION OF J(3,5) CORRECT?))**

The Jacobian matrix for our example contains 18 non-zero elements:

.. math::

   \begin{aligned}
     \mathbf{J}(1,1) & = & - k_{6}\, [\chem{M}] - k_{7}\, [\chem{O_3}]\\
     \mathbf{J}(1,3) & = & k_{5} - k_{7}\, [\chem{O(^1D)}]\\
     \mathbf{J}(2,1) & = & k_{6}\, [\chem{M}]\\
     \mathbf{J}(2,2) & = & - k_{2}\, [\chem{O_2}] - k_{4}\, [\chem{O_3}] 
                           - k_{9}\, [\chem{NO_2}]\\
     \mathbf{J}(2,3) & = & k_{3} - k_{4}\, [\chem{O}]\\
     \mathbf{J}(2,5) & = & - k_{9}\, [\chem{O}] + k_{10}\\
     \mathbf{J}(3,1) & = & - k_{7}\, [\chem{O_3}]\\
     \mathbf{J}(3,2) & = & k_{2}\, [\chem{O_2}] - k_{4}\, [\chem{O_3}]\\
     \mathbf{J}(3,3) & = & - k_{3} - k_{4}\, [\chem{O}] - k_{5} - k_{7}\, 
                           [\chem{O(^1D)}] - k_{8}\, [\chem{NO}]\\
     \mathbf{J}(3,4) & = & - k_{8}\, [\chem{O_3}]\\
     \mathbf{J}(4,2) & = & k_{9}\, [\chem{NO_2}]\\
     \mathbf{J}(4,3) & = & - k_{8}\, [\chem{NO}]\\
     \mathbf{J}(4,4) & = & - k_{8}\, [\chem{O_3}]\\
     \mathbf{J}(4,5) & = & k_{9}\, [\chem{O}] + k_{10}\\
     \mathbf{J}(5,2) & = & - k_{9}\, [\chem{NO_2}]\\
     \mathbf{J}(5,3) & = & k_{8}\, [\chem{NO}]\\
     \mathbf{J}(5,4) & = & k_{8}\, [\chem{O_3}]\\
     \mathbf{J}(5,5) & = & - k_{9}\, [\chem{O}] - k_{10}\\\end{aligned}

It defines how the temporal change of each chemical species depends on
all other species. For example, :math:`\mathbf{J}(5,2)` shows that
(species number 5) is affected by (species number 2) via reaction number
R9. The sparse data structures for the Jacobian are declared and
initialized in root. The code for the ODE Jacobian and sparse
multiplications is in root. The Jacobian of the ODE function is
automatically constructed by KPP. KPP generates the Jacobian subroutine
or where the latter is generated when the sparse format is required.
Using the variable species , the fixed species , and the rate
coefficients as input, the subroutine calculates the Jacobian . The
default data structures for the sparse compressed on rows Jacobian
representation are shown in Table `5 <#tab:sparse-jac>`__ (for the case
where the LU fill-in is accounted for). stores the elements of the
Jacobian in row order. Each row starts at position , and . The location
of the -th diagonal element is . The sparse element is the Jacobian
entry in row and column . For the example KPP generates the following
Jacobian sparse data structure:

::

   LU_ICOL = (/ 1,3,1,2,3,5,1,2,3,4, &
               5,2,3,4,5,2,3,4,5 /)
   LU_IROW = (/ 1,1,2,2,2,2,3,3,3,3, &
               3,4,4,4,4,5,5,5,5 /)
   LU_CROW = (/ 1,3,7,12,16,20 /)
   LU_DIAG = (/ 1,4,9,14,19,20 /)

This is visualized in Fig. `1 <#fig:jac>`__. The sparsity coordinate
vectors are computed by KPP and initialized statically. These vectors
are constant as the sparsity pattern of the Jacobian does not change
during the computation.

.. container::
   :name: tab:sparse-jac

   .. table:: Sparse Jacobian Data Structures

      =============== =========================
      Global variable Represents
      \               Jacobian nonzero elements
      \               Row indices
      \               Column indices
      \               Start of rows
      \               Diagonal entries
      \               
      =============== =========================

Two other KPP-generated routines, and are useful for direct and adjoint
sensitivity analysis. They perform sparse multiplication of (or its
transpose for ) with the user-supplied vector without any indirect
addressing.

.. _`sec:output-ode-hess`:

root\ ``_Hessian.f90`` and root\ ``_HessianSP.f90``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. container:: figure*

   .. image:: Figures/small_hess1.pdf
      :alt: image
      :width: 80.0%

The sparse data structures for the Hessian are declared and initialized
in root. The Hessian function and associated sparse multiplications are
in root. The Hessian contains the second order derivatives of the time
derivative functions. More exactly, the Hessian is a 3-tensor such that

.. math::

   H_{i,j,k} = \frac{\partial^2 ({\mathrm{d}}c/{\mathrm{d}}t)_i}{\partial c_j \,\partial c_k}~,
     \qquad 1 \le i,j,k \le N_{\rm var}~.
   \label{eqn:Hessian1}

KPP generates the routine . Using the variable species , the fixed
species , and the rate coefficients as input, the subroutine calculates
the Hessian. The Hessian is a very sparse tensor. The sparsity of the
Hessian for our example is visualized in
Fig. `[fig:hess1] <#fig:hess1>`__. KPP computes the number of nonzero
Hessian entries and saves it in the variable . The Hessian itself is
represented in coordinate sparse format. The real vector holds the
values, and the integer vectors , , and hold the indices of nonzero
entries as illustrated in Table `6 <#tab:sparse-hess>`__. Since the time
derivative function is smooth, these Hessian matrices are symmetric,
:math:`_{i,j,k}`\ =\ :math:`_{i,k,j}`. KPP stores only those entries
:math:`_{i,j,k}` with :math:`j \le k`. The sparsity coordinate vectors ,
, and are computed by KPP and initialized statically. They are constant
as the sparsity pattern of the Hessian does not change during the
computation.

.. container::
   :name: tab:sparse-hess

   .. table:: Sparse Hessian Data

      ======== ===========================================
      Variable Represents
      \        Hessian nonzero elements :math:`_{i,j,k}`
      \        Index :math:`i` of element :math:`_{i,j,k}`
      \        Index :math:`j` of element :math:`_{i,j,k}`
      \        Index :math:`k` of element :math:`_{i,j,k}`
      \        
      ======== ===========================================

The routines and compute the action of the Hessian (or its transpose) on
a pair of user-supplied vectors and . Sparse operations are employed to
produce the result vector.

.. _`sec:output-la`:

root\ ``_LinearAlgebra.f90``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Sparse linear algebra routines are in the file root. To numerically
solve for the chemical concentrations one must employ an implicit
timestepping technique, as the system is usually stiff. Implicit
integrators solve systems of the form

.. math:: P\, x = (I - h \gamma J)\, x = b

where the matrix :math:`P=I - h \gamma J` is refered to as the
“prediction matrix”. :math:`I` the identity matrix, :math:`h` the
integration time step, :math:`\gamma` a scalar parameter depending on
the method, and :math:`J` the system Jacobian. The vector :math:`b` is
the system right hand side and the solution :math:`x` typically
represents an increment to update the solution.

The chemical Jacobians are typically sparse, i.e. only a relatively
small number of entries are nonzero. The sparsity structure of :math:`P`
is given by the sparsity structure of the Jacobian, and is produced by
KPP (with account for the fill-in) as discussed above.

KPP generates the sparse linear algebra subroutine which performs an
in-place, non-pivoting, sparse LU decomposition of the prediction matrix
:math:`P`. Since the sparsity structure accounts for fill-in, all
elements of the full LU decomposition are actually stored. The output
argument returns a value that is nonzero if singularity is detected.

The subroutines and use the in-place LU factorization :math:`P` as
computed by and perform sparse backward and forward substitutions (using
:math:`P` or its transpose). The sparse linear algebra routines and are
extremely efficient, as shown by :raw-latex:`\citep{IMPLEMENTATION}`.

.. _`sec:output-stoichiom`:

root\ ``_Stoichiom.f90`` and root\ ``_StoichiomSP.f90``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

These files contain a description of the chemical mechanism in
stoichiometric form. The file root contains the functions for reactant
products and its Jacobian, and derivatives with respect to rate
coefficients. The declaration and initialization of the stoichiometric
matrix and the associated sparse data structures is done in root.

The stoichiometric matrix is constant sparse. For our example the matrix
has 22 nonzero entries out of 50 entries. KPP produces the
stoichiometric matrix in sparse, column-compressed format, as shown in
Table `7 <#tab:sparse-stoicm>`__. Elements are stored in columnwise
order in the one-dimensional vector of values . Their row and column
indices are stored in and respectively. The vector contains pointers to
the start of each column. For example column starts in the sparse vector
at position and ends at . The last value simplifies the handling of
sparse data structures.

.. container::
   :name: tab:sparse-stoicm

   .. table:: Sparse Stoichiometric Matrix

      =============== =====================
      Global variable Represents
      \               Stoichiometric matrix
      \               Row indices
      \               Column indices
      \               Start of columns
      \               
      =============== =====================

The subroutine computes the reactant products for each reaction, and the
subroutine computes the Jacobian of reactant products vector, i.e.:

.. math:: {\tt JVRP} = \partial {\tt ARP} / \partial {\tt V}

The matrix is sparse and is computed and stored in row compressed sparse
format, as shown in Table `8 <#tab:sparse-jvrp>`__. The parameter holds
the number of nonzero elements. For our example:

::

   NJVRP = 13
   CROW_JVRP = (/ 1,1,2,3,5,6,7,9,11,13,14 /)
   ICOL_JVRP = (/ 2,3,2,3,3,1,1,3,3,4,2,5,4 /)

.. container::
   :name: tab:sparse-jvrp

   .. table:: Sparse Data for Jacobian of Reactant Products

      =============== ========================
      Global variable Represents
      \               Nonzero elements of JVRP
      \               Column indices in JVRP
      \               Row indices in JVRP
      \               Start of rows in JVRP
      \               
      =============== ========================

If is set to , the stoichiometric formulation allows a direct
computation of the derivatives with respect to rate coefficients.

The subroutine computes the partial derivative of the ODE function with
respect to a subset of reaction coefficients, whose indices are
specifies in the array

.. math:: {\tt DFDR} = \partial {\tt Vdot} / \partial {\tt RCT(JCOEFF)}

Similarly one can obtain the partial derivative of the Jacobian with
respect to a subset of the rate coefficients. More exactly, KPP
generates the subroutine which calculates , the product of this partial
derivative with a user-supplied vector:

.. math:: {\tt DJDR} = [ \partial {\tt JVS} / \partial {\tt RCT(JCOEFF)} ] \times {\tt U}

.. _`sec:output-stochastic`:

root\ ``_Stochastic.f90``
~~~~~~~~~~~~~~~~~~~~~~~~~

If the generation of stochastic functions is switched on, KPP produces
the file root with the following functions:

calculates the propensity vector. The propensity function uses the
number of molecules of variable () and fixed () species, as well as the
stochastic rate coefficients () to calculate the vector of propensity
rates (). The propensity :math:`_j` defines the probability that the
next reaction in the system is the :math:`j^{th}` reaction.

converts deterministic rates to stochastic. The stochastic rate
coefficients () are obtained through a scaling of the deterministic rate
coefficients (). The scaling depends on the of the reaction container
and on the number of molecules which react.

calculates changes in the number of molecules. When the reaction with
index takes place, the number of molecules of species involved in that
reaction changes. The total number of molecules is updated by the
function.

These functions are used by the Gillespie numerical integrators (direct
stochastic simulation algorithm). These integrators are provided in both
Fortran90 and C implementations (the template file name is ). Drivers
for stochastic simulations are also implemented (the template file name
is ).

.. _`sec:output-utility`:

root\ ``_Util.f90``
~~~~~~~~~~~~~~~~~~~

The utility and input/output functions are in root. In addition to the
chemical system description routines discussed above, KPP generates
several utility routines, some of which are summarized in
Table `[tab:functions] <#tab:functions>`__.

The subroutines , , and can be used to print the concentration of the
species that were selected with to the file root.

.. _`sec:output-mexcode`:

root\ ``_mex_Fun.f90``, root\ ``_mex_Jac_SP.f90``, and root\ ``_mex_Hessian.f90``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Mex is a *M*\ atlab *ex*\ tension. KPP generates the mex routines for
the ODE function, Jacobian, and Hessian, for the target languages C,
Fortran77, and Fortran90. After compilation (using Matlab’s mex
compiler) the mex functions can be called instead of the corresponding
Matlab m-functions. Since the calling syntaxes are identical, the user
only has to insert the string within the corresponding function name.
Replacing m-functions by mex-functions gives the same numerical results,
but the computational time could be considerably smaller, especially for
large kinetic systems.

If possible we recommend to build mex files using the C language, as
Matlab offers most mex interface options for the C language. Moreover,
Matlab distributions come with a native C compiler (lcc) for building
executable functions from mex files. Fortran77 mex files work well on
most platforms without additional efforts. However, the mex files built
using Fortran90 may require further platform-specific tuning of the mex
compiler options.

.. _`sec:output-makefile`:

The Makefile
~~~~~~~~~~~~

KPP produces a Makefile that allows for an easy compilation of all
KPP-generated source files. The file name is . The Makefile assumes that
the selected driver contains the main program. However, if no driver was
selected (i.e. ), it is necessary to add the name of the main program
file manually to the Makefile.

.. _`sec:c`:

The C Code
----------

The driver file root contains the main (driver) and numerical integrator
functions, as well as declarations and initializations of global
variables. The generated C code includes three header files which are -d
in other files as appropriate. The global parameters
(Table `3 <#tab:parameters>`__) are -d in the header file root. The
global variables (Table `4 <#tab:global>`__) are extern-declared in
root, and declared in the driver file root. The header file root
contains extern declarations of sparse data structures for the Jacobian
(Table `5 <#tab:sparse-jac>`__), Hessian
(Table `6 <#tab:sparse-hess>`__), stoichiometric matrix
(Table `7 <#tab:sparse-stoicm>`__), and the Jacobian of reaction
products (Table `8 <#tab:sparse-jvrp>`__). The actual declarations of
each data structures is done in the corresponding files.

The code for the ODE function (Sect. `4.1.10 <#sec:output-ode-fun>`__)
is in root. The code for the ODE Jacobian and sparse multiplications
(Sect. `4.1.11 <#sec:output-ode-jac>`__) is in root, and the declaration
and initialization of the Jacobian sparse data structures
(Table `5 <#tab:sparse-jac>`__) is in the file root. Similarly, the
Hessian function and associated sparse multiplications
(Section `4.1.12 <#sec:output-ode-hess>`__) are in root, and the
declaration and initialization of Hessian sparse data structures (Table
`6 <#tab:sparse-hess>`__) in root.

The file root contains the functions for reactant products and its
Jacobian, and derivatives with respect to rate coefficients
(Sect. `4.1.14 <#sec:output-stoichiom>`__). The declaration and
initialization of the stoichiometric matrix and the associated sparse
data structures (Tables `7 <#tab:sparse-stoicm>`__ and
`8 <#tab:sparse-jvrp>`__) is done in root.

Sparse linear algebra routines (Sect. `4.1.13 <#sec:output-la>`__) are
in the file root. The code to update the rate constants and user defined
code for rate laws is in root.

Various utility and input/output functions
(Sect. `4.1.16 <#sec:output-utility>`__) are in root and root.

Finally, mex gateway routines that allow the C implementation of the ODE
function, Jacobian, and Hessian to be called directly from Matlab
(Sect. `4.1.17 <#sec:output-mexcode>`__) are also generated (in the
files root, root, and root).

.. _`sec:f77`:

The Fortran77 Code
------------------

The general layout of the Fortran77 code is similar to the layout of the
C code. The driver file root contains the main (driver) and numerical
integrator functions.

The generated Fortran77 code includes three header files. The global
parameters (Table `3 <#tab:parameters>`__) are defined as parameters and
initialized in the header file root. The global variables
(Table `4 <#tab:global>`__) are declared in root as common block
variables. There are global common blocks for real (), integer (), and
character () global data. They can be accessed from within each program
unit that includes the global header file.

The header file root contains common block declarations of sparse data
structures for the Jacobian (Table `5 <#tab:sparse-jac>`__), Hessian
(Table `6 <#tab:sparse-hess>`__), stoichiometric matrix
(Table `7 <#tab:sparse-stoicm>`__), and the Jacobian of reaction
products (Table `8 <#tab:sparse-jvrp>`__). These sparse data structures
are initialized in four named block data statements: (in root), (in
root), and (in root).

The code for the ODE function (Sect. `4.1.10 <#sec:output-ode-fun>`__)
is in root. The code for the ODE Jacobian and sparse multiplications
(Sect. `4.1.11 <#sec:output-ode-jac>`__) is in root. The Hessian
function and associated sparse multiplications
(Sect. `4.1.12 <#sec:output-ode-hess>`__) are in root.

The file root contains the functions for reactant products and its
Jacobian, and derivatives with respect to rate coefficients
(Sect. `4.1.14 <#sec:output-stoichiom>`__). The declaration and
initialization of the stoichiometric matrix and the associated sparse
data structures (Tables `7 <#tab:sparse-stoicm>`__ and
`8 <#tab:sparse-jvrp>`__) is done in the block data statement.

Sparse linear algebra routines (Sect. `4.1.13 <#sec:output-la>`__) are
in the file root. The code to update the rate constants is in root, and
the utility and input/output functions
(Sect. `4.1.16 <#sec:output-utility>`__) are in root and root.

Matlab-mex gateway routines for the ODE function, Jacobian, and Hessian
are discussed in Sect. `4.1.17 <#sec:output-mexcode>`__.

.. _`sec:matlab`:

The Matlab Code
---------------

Matlab (http://www.mathworks.com/products/matlab/) provides a high-level
programming environment that allows algorithm development, numerical
computations, and data analysis and visualization. The KPP-generated
Matlab code allows for a rapid prototyping of chemical kinetic schemes,
and for a convenient analysis and visualization of the results.
Differences between different kinetic mechanisms can be easily
understood. The Matlab code can be used to derive reference numerical
solutions, which are then compared against the results obtained with
user-supplied numerical techniques. Last but not least Matlab is an
excellent environment for educational purposes. KPP/Matlab can be used
to teach students fundamentals of chemical kinetics and chemical
numerical simulations.

Each Matlab function has to reside in a separate m-file. Function calls
use the m-function-file names to reference the function. Consequently,
KPP generates one m-function-file for each of the functions discussed in
Sections `4.1.10 <#sec:output-ode-fun>`__,
`4.1.11 <#sec:output-ode-jac>`__, `4.1.12 <#sec:output-ode-hess>`__,
`4.1.13 <#sec:output-la>`__, `4.1.14 <#sec:output-stoichiom>`__, and
`4.1.16 <#sec:output-utility>`__. The names of the m-function-files are
the same as the names of the functions (prefixed by the model name
root).

The Matlab syntax for calling each function is

::

   [Vdot] = Fun    (V, F, RCT);
   [JVS ] = Jac_SP (V, F, RCT);
   [HESS] = Hessian(V, F, RCT);

The global parameters (Table `3 <#tab:parameters>`__) are defined as
Matlab variables and initialized in the file root. The variables of
Table `4 <#tab:global>`__ are declared as Matlab variables in the file
root. They can be accessed from within each Matlab function by using
declarations of the variables of interest.

The sparse data structures for the Jacobian
(Table `5 <#tab:sparse-jac>`__), the Hessian
(Table `6 <#tab:sparse-hess>`__), the stoichiometric matrix
(Table `7 <#tab:sparse-stoicm>`__), and the Jacobian of reaction
products (Table `8 <#tab:sparse-jvrp>`__) are declared as Matlab
variables in the file root. They are initialized in separate m-files,
namely root root, and root respectively.

Two wrappers (root and root) are provided for interfacing the ODE
function and the sparse ODE Jacobian with Matlab’s suite of ODE
integrators. Specifically, the syntax of the wrapper calls matches the
syntax required by Matlab’s integrators like ode15s. Moreover, the
Jacobian wrapper converts the sparse KPP format into a Matlab sparse
matrix.

.. container:: table*

   .. container:: center

      ==== ============================================================
      File Description
      root driver
      root Global parameters
      root Global variables
      root Global monitor variables
      root Global sparsity data
      root Template for ODE function
      root ODE function
      root Template for ODE Jacobian
      root ODE Jacobian in sparse format
      root Sparsity data structures
      root ODE Hessian in sparse format
      root Sparsity data structures
      root Hessian action on vectors
      root Transposed Hessian action on vectors
      root Derivatives of Fun and Jac with respect to rate coefficients
      root Sparse data
      root Reactant products
      root Jacobian of reactant products
      root User-defined reaction rate laws
      root Update photolysis rate coefficients
      root Update all rate coefficients
      root Update solar intensity
      root Check mass balance for selected atoms
      root Set initial values
      root Shuffle concentration vector
      root Shuffle concentration vector
      \    
      ==== ============================================================

.. _`sec:output-map`:

The map file
------------

The map file root contains a summary of all the functions, subroutines
and data structures defined in the code file, plus a summary of the
numbering and category of the species involved.

This file contains supplementary information for the user. Several
statistics are listed here, like the total number equations, the total
number of species, the number of variable and fixed species. Each
species from the chemical mechanism is then listed followed by its type
and numbering.

Furthermore it contains the complete list of all the functions generated
in the target source file. For each function, a brief description of the
computation performed is attached containing also the meaning of the
input and output parameters.

.. _`sec:developer-info`:

Information for KPP Developers
==============================

This chapter is meant for KPP Developers. It describes the internal
architecture of the KPP preprocessor, the basic modules and their
functionalities, and the preprocessing analysis performed on the input
files. KPP can be very easily configured to suit a broad class of users.

.. _`sec:directory-structure`:

KPP directory structure
-----------------------

.. container:: center

   .. container::
      :name: tab:source

      .. table:: Source code files

         ==== =================================
         File Description
         \    main program
         \    generic code generation functions
         \    header file
         \    generation of C code
         \    generation of Fortran77 code
         \    generation of Fortran90 code
         \    generation of matlab code
         \    debugging output
         \    header file
         \    header file
         \    generic code generation functions
         \    flex/bison-generated file
         \    input for flex and bison
         \    input for flex
         \    input for bison
         \    evaluate parsed input
         \    evaluate parsed input
         \    flex/bison-generated file
         \    flex/bison-generated header file
         \    
         ==== =================================

The KPP distribution will unfold a directory ``$KPP_HOME`` with the
following subdirectories:

-  **src/** Contains the KPP source code files, as listed in
   Table `9 <#tab:source>`__.

-  **bin/** Contains the KPP executable. The path to this directory
   needs to be added to the environment variable.

-  **util/** Contains different function templates useful for the
   simulation. Each template file has a suffix that matches the
   appropriate target language (, , , or ). KPP will run the template
   files through the substitution preprocessor. The user can define
   their own auxiliary functions by inserting them into the files.

-  **models/** Contains the description of the chemical models. Users
   can define their own models by placing the model description files in
   this directory. The KPP distribution contains several models from
   atmospheric chemistry which can be used as templates for model
   definitions.

-  **drv/** Contains driver templates for chemical simulations. Each
   driver has a suffix that matches the appropriate target language (, ,
   , or ). KPP will run the appropriate driver through the substitution
   preprocessor. The driver template provided with the distribution
   works with any example. Users can define here their own driver
   templates.

-  **int/** Contains numerical time stepping (integrator) routines. The
   command “*integrator*” will force KPP to look into this directory for
   a definition file *integrator*. This file selects the numerical
   routine (with the command) and sets the function type, the Jacobian
   sparsity type, the target language, etc. Each integrator template is
   found in a file that ends with the appropriate suffix (, , , or ).
   The selected template is processed by the substitution preprocessor.
   Users can define here their own numerical integration routines.

-  **examples/** Contains several model description examples ( files)
   which can be used as templates for building simulations with KPP.

-  **site-lisp/** Contains the file which provides a KPP mode for emacs
   with color highlighting.

KPP environment variables
-------------------------

.. container:: table*

   .. container:: center

      +----------------+-------------------------+----------------------+
      | Variable       | Description             | Default assumed      |
      +----------------+-------------------------+----------------------+
      | ``$KPP_HOME``  | Required, stores the    | no default           |
      |                | absolute path to the    |                      |
      |                | KPP distribution        |                      |
      +----------------+-------------------------+----------------------+
      | ``$KPP_MODEL`` | Optional, specifies     | ``$KPP_HOME/models`` |
      |                | additional places were  |                      |
      |                | KPP will look for model |                      |
      |                | files before searching  |                      |
      |                | the default             |                      |
      +----------------+-------------------------+----------------------+
      | ``$KPP_INT``   | Optional, specifies     | ``$KPP_HOME/int``    |
      |                | additional places were  |                      |
      |                | KPP will look for       |                      |
      |                | integrator files before |                      |
      |                | searching the default.  |                      |
      +----------------+-------------------------+----------------------+
      | ``$KPP_DRV``   | Optional, specifies     | ``$KPP_HOME/drv``    |
      |                | additional places were  |                      |
      |                | KPP will look for       |                      |
      |                | driver files before     |                      |
      |                | searching the default   |                      |
      +----------------+-------------------------+----------------------+
      |                |                         |                      |
      +----------------+-------------------------+----------------------+

In order for KPP to find its components, it has to know the path to the
location where the KPP distribution is installed. This is achieved by
requiring the ``$KPP_HOME`` environment variable to be set to the path
where KPP is installed.

The PATH variable should be updated to contain the ``$KPP_HOME/bin``
directory.

There are several optional environment variable that control the places
where KPP looks for module files, integrators, and drivers. They are all
summarized in Table `[tab:environment] <#tab:environment>`__.

KPP internal modules
--------------------

Scanner and Parser
~~~~~~~~~~~~~~~~~~

This module is responsible for reading the kinetic description files and
extracting the information necessary in the code generation phase. We
make use of the flex and bison generic tools in implementing our own
scanner and parser. Using these tools this module gathers information
from the input files and fills in the following data structures in
memory:

-  The atom list

-  The species list

-  The left hand side matrix of coefficients

-  The right hand side matrix of coefficients

-  The equation rates

-  The option list

Error checking is performed at each step in the scanner and the parser.
For each syntax error the exact line and input file, along with an
appropriate error message are produced. In most of the cases the exact
cause of the error can be identified, therefore the error messages are
very precise. Some other errors like mass balance, and equation
duplicates, are tested at the end of this phase.

Species reordering
~~~~~~~~~~~~~~~~~~

When parsing the input files, the species list is updated as soon as a
new species is encountered in a chemical equation. Therefore the
ordering of the species is the order in which they appear in the
equation description section. This is not a useful order for subsequent
operations. The species have to be first sorted such that all variable
species and all fixed species are put together. Then if a sparsity
structure of the Jacobian is required, it might be better to reorder the
species in such a way that the factorization of the Jacobian will
preserve the sparsity. This reordering is done using a Markovitz type of
algorithm.

Expression trees computation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This is the core of the preprocessor. This module has to generate the
production/destruction functions the Jacobian and all the data structure
nedeed by these functions. This module has to build a language
independent structure of each function and statement in the target
source file. Instead of using an intermediate format for this as some
other compilers do, KPP generates the intermediate format for just one
statement at a time. The vast majority of the statements in the target
source file are assignments. The expression tree for each assignment is
incrementally build by scanning the coefficient matrices and the rate
constant vector. At the end these expression trees are simplified.
Similar approaches are applied to function declaration and prototypes,
data declaration and initialization.

Code generation
~~~~~~~~~~~~~~~

There are basically two modules, each dealing with the syntax
particularities of the target language. For example, the C module
includes a function that generates a valid C assignment when given an
expression tree. Similarly there are functions for data declaration,
initializations, comments, function prototypes, etc. Each of these
functions produce the code into an output buffer. A language specific
routine reads from this buffer and splits the statements into lines to
improve readability of the generated code.

Adding new KPP commands
-----------------------

To add a new KPP command, the source code has to be edited at several
locations. A short summary is presented here, using the new command as
an example:

-  Add to several files in the directory:

   -  : ``void CmdNEWCMD( char *cmd );``

   -  : ``{ "NEWCMD", PRM_STATE, NEWCMD },``

   -  : ``void CmdNEWCMD( char *cmd )``

   -  :

      -  ``%token NEWCMD``

      -  ``NEWCMD PARAMETER``

      -  ``{ CmdNEWCMD( $2 ); }``

-  Maybe add a CI-test:

   -  Create a new directory

   -  Add new CI-test to and in the directory

-  Other:

   -  Explain in user manual :

      -  Add to table

      -  Add a section

      -  Add to BNF description table

   -  Add to

Numerical methods
=================

.. container:: table*

   .. container:: center

      +----------------------------------+----------------------------------+
      | Symbol                           | Description                      |
      +----------------------------------+----------------------------------+
      | :math:`s`                        | Number of stages                 |
      +----------------------------------+----------------------------------+
      | :math:`t^n`                      | Discrete time moment             |
      +----------------------------------+----------------------------------+
      | :math:`h`                        | Time step :math:`h=t^{n+1}-t^n`  |
      +----------------------------------+----------------------------------+
      | :math:`y^n`                      | Numerical solution               |
      |                                  | (concentration) at :math:`t^n`   |
      +----------------------------------+----------------------------------+
      | :math:`\delta y^n`               | tangent linear solution at       |
      |                                  | :math:`t^n`                      |
      +----------------------------------+----------------------------------+
      | :math:`\lambda^n`                | Adjoint numerical solution at    |
      |                                  | :math:`t^n`                      |
      +----------------------------------+----------------------------------+
      | :math:`f(\cdot,\cdot)`           | The ODE derivative function:     |
      |                                  | :math:`y'=f(t,y)`                |
      +----------------------------------+----------------------------------+
      | :math:`f_t(\cdot,\cdot)`         | Partial time derivative          |
      |                                  | :math:`f_t(                      |
      |                                  | t,y)=\partial f(t,y)/\partial t` |
      +----------------------------------+----------------------------------+
      | :math:`J(\cdot,\cdot)`           | The Jacobian                     |
      |                                  | :math:`J(                        |
      |                                  | t,y)=\partial f(t,y)/\partial y` |
      +----------------------------------+----------------------------------+
      | :math:`J_t(\cdot,\cdot)`         | Partial time derivative of       |
      |                                  | Jacobian                         |
      |                                  | :math:`J_t(                      |
      |                                  | t,y)=\partial J(t,y)/\partial t` |
      +----------------------------------+----------------------------------+
      | :math:`A`                        | The system matrix                |
      +----------------------------------+----------------------------------+
      | :math:`H(\cdot,\cdot)`           | The Hessian                      |
      |                                  | :math:`H(t,y)                    |
      |                                  | =\partial^2 f(t,y)/\partial y^2` |
      +----------------------------------+----------------------------------+
      | :math:`T_i`                      | Internal stage time moment for   |
      |                                  | Runge-Kutta and Rosenbrock       |
      |                                  | methods                          |
      +----------------------------------+----------------------------------+
      | :math:`Y_i`                      | Internal stage solution for      |
      |                                  | Runge-Kutta and Rosenbrock       |
      |                                  | methods                          |
      +----------------------------------+----------------------------------+
      | :math:`k_i`, :math:`\ell_i`,     | Internal stage vectors for       |
      | :math:`u_i`, :math:`v_i`         | Runge-Kutta and Rosenbrock       |
      |                                  | methods, their tangent linear    |
      |                                  | and adjoint models               |
      +----------------------------------+----------------------------------+
      | :math:`\alpha_i`,                | Method coefficients              |
      | :math:`\alpha_{ij}`,             |                                  |
      | :math:`a_{ij}`, :math:`b_i`,     |                                  |
      | :math:`c_i`, :math:`c_{ij}`,     |                                  |
      | :math:`e_i`, :math:`m_i`         |                                  |
      +----------------------------------+----------------------------------+
      |                                  |                                  |
      +----------------------------------+----------------------------------+

The KPP numerical library contains a set of numerical integrators
selected to be very efficient in the low to medium accuracy regime
(relative errors :math:`\sim 10^{-2} \dots 10^{-5}`). In addition, the
KPP numerical integrators preserve the linear invariants (i.e., mass) of
the chemical system.

KPP implements several Rosenbrock methods: ROS–2
:raw-latex:`\citep{Verwer99}`, ROS–3 :raw-latex:`\citep{BENCHMARK-2}`,
RODAS–3 :raw-latex:`\citep{BENCHMARK-2}`, ROS–4
:raw-latex:`\citep{k:HW2}`, and RODAS–4 :raw-latex:`\citep{k:HW2}`. For
each of them KPP implements the tangent linear model (direct decoupled
sensitivity) and the adjoint models. The implementations distinguish
between sensitivities with respect to initial values and sensitivities
with respect to parameters for efficiency.

Note that KPP produces the building blocks for the simulation and also
for the sensitivity calculations. It also provides application
programming templates. Some minimal programming may be required from the
users in order to construct their own application from the KPP building
blocks.

.. container:: table*

   .. container:: center

      .. container:: tabular

         | lp11cm Variable & Description
         | & = 1: :math:`F = F(y)`, i.e. independent of t (autonomous)
         | & = 0: :math:`F = F(t,y)`, i.e. depends on t (non-autonomous)
         | & (only available for some of the integrators)
         | & The absolute () and relative () tolerances can be expressed
           by either a scalar or individually for each species in a
           vector:
         | & 0: -dimensional vector
         | & 1: scalar
         | & Selection of a specific method (only available for some of
           the integrators).
         | & Maximum number of integration steps.
         | & Maximum number of Newton iterations (only available for
           some of the integrators).
         | & Starting values of Newton iterations:
         | & 0: interpolated
         | & 1: zero
         | & (only available for some of the integrators)
         | & This determines which subroutines are called within the
           integrator:

         +---------------+-----------------------------------------------------+
         | -1:           | Do not call any subroutines                         |
         +---------------+-----------------------------------------------------+
         | 0:            | Use the integrator-specific default values          |
         +---------------+-----------------------------------------------------+
         | :math:`>`\ 1: | A number between 1 and 7, derived by adding up the  |
         |               | bits with values 4, 2, and 1. The first digit (4)   |
         |               | activates . The second digit (2) activates . The    |
         |               | last digit (1) activates .                          |
         +---------------+-----------------------------------------------------+

         | For example, (4+2) will activate the calls to and but not the
           call to .

         | & , the lower bound of the integration step size. It is not
           recommended to change the default value of zero.
         | & , the upper bound of the integration step size.
         | & , the starting value of the integration step size.
         | & , lower bound on step decrease factor.
         | & , upper bound on step increase factor.
         | & , step decrease factor after multiple rejections.
         | & , the factor by which the new step is slightly smaller than
           the predicted value.
         | & . If the Newton convergence rate is smaller than , the
           Jacobian is not recomputed (only available for some of the
           integrators).
         | & , the stopping criterion for Newton’s method (only
           available for some of the integrators).
         | & (only available for some of the integrators).
         | & . If :math:`<` / :math:`<` , then the step size is kept
           constant and the LU factorization is reused (only available
           for some of the integrators).

.. container:: table*

   .. container:: center

      +----------+----------------------------------------------------------+
      | Variable | Description                                              |
      +----------+----------------------------------------------------------+
      |          | Number of function calls.                                |
      +----------+----------------------------------------------------------+
      |          | Number of Jacobian calls.                                |
      +----------+----------------------------------------------------------+
      |          | Number of steps.                                         |
      +----------+----------------------------------------------------------+
      |          | Number of accepted steps.                                |
      +----------+----------------------------------------------------------+
      |          | Number of rejected steps (except at very beginning).     |
      +----------+----------------------------------------------------------+
      |          | Number of LU decompositions.                             |
      +----------+----------------------------------------------------------+
      |          | Number of forward/backward substitutions.                |
      +----------+----------------------------------------------------------+
      |          | Number of singular matrix decompositions.                |
      +----------+----------------------------------------------------------+
      |          | , the time corresponding to the computed upon return.    |
      +----------+----------------------------------------------------------+
      |          | , the last accepted step before exit.                    |
      +----------+----------------------------------------------------------+
      |          | , the last predicted step (not yet taken). For multiple  |
      |          | restarts, use as in the subsequent run.                  |
      +----------+----------------------------------------------------------+
      |          |                                                          |
      +----------+----------------------------------------------------------+

In order to offer more control over the integrator, the KPP-generated
subroutine provides the optional input parameters and . Each of them is
an array of 20 elements that allow the fine-tuning of the integrator, as
shown in Table `[tab:control] <#tab:control>`__. Similarly, to obtain
more information about the integration, the subroutine provides the
optional output parameters and . They are both arrays of 20 elements, as
shown in Table `[tab:status] <#tab:status>`__.

In the following sections we introduce the numerical methods implemented
in KPP. The symbols used in the formulas are explained in
Table `[tab:symbols] <#tab:symbols>`__.

Rosenbrock Methods
------------------

An :math:`s`-stage Rosenbrock method
:raw-latex:`\cite[Section IV.7]{k:HW2}` computes the next-step solution
by the formulas

.. math::

   \begin{aligned}
   \label{eqn:altRosenbrock}
   y^{n+1} &=& y^n + \sum_{i=1}^s m_i k_i~,
   \quad {\rm Err}^{n+1} = \sum_{i=1}^s e_i k_i\\
   \nonumber
   T_i &=& t^n + \alpha_i h~, \quad
   Y_i =y^n + \sum_{j=1}^{i-1} a_{ij} k_j~,\\
   \nonumber
   A &=& \left[ \frac{1}{h \gamma} - J^T(t^n,y^n) \right]\\
   \nonumber
   A \cdot k_i &=&  f\left( \, T_i,
   \, Y_i \,\right) + \sum_{j=1}^{i-1} \frac{c_{ij}}{h} k_j + h \gamma_i
   f_t\left(t^n,y^n\right)~.\end{aligned}

where :math:`s` is the number of stages,
:math:`\alpha_i = \sum_j \alpha_{ij}` and
:math:`\gamma_i = \sum_j \gamma_{ij}`. The formula coefficients
(:math:`a_{ij}` and :math:`\gamma_{ij}`) give the order of consistency
and the stability properties. :math:`A` is the system matrix (in the
linear systems to be solved during implicit integration, or in the
Newton’s method used to solve the nonlinear systems). It is the scaled
identity matrix minus the Jacobian.

The coefficients of the methods implemented in KPP are shown in
Table `[tab:Rosenbrock] <#tab:Rosenbrock>`__.

.. container:: table*

   .. container:: center

      +---------+-----------+----------+-------+-----------+-----------+
      |         | Stages    | Function | Order | Stability | Method    |
      +---------+-----------+----------+-------+-----------+-----------+
      | name    | (:        | calls    |       | p         | coe       |
      |         | math:`s`) |          |       | roperties | fficients |
      +---------+-----------+----------+-------+-----------+-----------+
      | ROS–2   | 2         | 2        | 2(1)  | L-stable  | :math     |
      |         |           |          |       |           | :`\gamma  |
      |         |           |          |       |           | = 1 + 1/\ |
      |         |           |          |       |           | sqrt{2}`, |
      |         |           |          |       |           | :math:`a_ |
      |         |           |          |       |           | {2,1} = 1 |
      |         |           |          |       |           | /\gamma`, |
      |         |           |          |       |           | :         |
      |         |           |          |       |           | math:`c_{ |
      |         |           |          |       |           | 2,1} = -2 |
      |         |           |          |       |           | /\gamma`, |
      |         |           |          |       |           | :math:`m  |
      |         |           |          |       |           | _1 = 3/(2 |
      |         |           |          |       |           | \gamma)`, |
      |         |           |          |       |           | :math:`m  |
      |         |           |          |       |           | _2 = 1/(2 |
      |         |           |          |       |           | \gamma)`, |
      |         |           |          |       |           | :math:`e  |
      |         |           |          |       |           | _1 = 1/(2 |
      |         |           |          |       |           | \gamma)`, |
      |         |           |          |       |           | :math:`e  |
      |         |           |          |       |           | _2 = 1/(2 |
      |         |           |          |       |           | \gamma)`, |
      |         |           |          |       |           | :ma       |
      |         |           |          |       |           | th:`\alph |
      |         |           |          |       |           | a_1 = 0`, |
      |         |           |          |       |           | :ma       |
      |         |           |          |       |           | th:`\alph |
      |         |           |          |       |           | a_2 = 1`, |
      |         |           |          |       |           | :math:`\  |
      |         |           |          |       |           | gamma_1 = |
      |         |           |          |       |           |  \gamma`, |
      |         |           |          |       |           | :math:`\  |
      |         |           |          |       |           | gamma_2 = |
      |         |           |          |       |           |  -\gamma` |
      +---------+-----------+----------+-------+-----------+-----------+
      | ROS–3   | 3         | 2        | 3(2)  | L-stable  | :m        |
      |         |           |          |       |           | ath:`a_{2 |
      |         |           |          |       |           | ,1} = 1`, |
      |         |           |          |       |           | :m        |
      |         |           |          |       |           | ath:`a_{3 |
      |         |           |          |       |           | ,1} = 1`, |
      |         |           |          |       |           | :m        |
      |         |           |          |       |           | ath:`a_{3 |
      |         |           |          |       |           | ,2} = 0`, |
      |         |           |          |       |           | :math:`   |
      |         |           |          |       |           | c_{2,1} = |
      |         |           |          |       |           |  -1.015`, |
      |         |           |          |       |           | :math:    |
      |         |           |          |       |           | `c_{3,1}  |
      |         |           |          |       |           | = 4.075`, |
      |         |           |          |       |           | :math:    |
      |         |           |          |       |           | `c_{3,2}  |
      |         |           |          |       |           | = 9.207`, |
      |         |           |          |       |           | :math:`   |
      |         |           |          |       |           | m_1 = 1`, |
      |         |           |          |       |           | :m        |
      |         |           |          |       |           | ath:`m_2  |
      |         |           |          |       |           | = 6.169`, |
      |         |           |          |       |           | :ma       |
      |         |           |          |       |           | th:`m_3 = |
      |         |           |          |       |           |  -0.427`, |
      |         |           |          |       |           | :math:`e_ |
      |         |           |          |       |           | 1 = 0.5`, |
      |         |           |          |       |           | :ma       |
      |         |           |          |       |           | th:`e_2 = |
      |         |           |          |       |           |  -2.908`, |
      |         |           |          |       |           | :m        |
      |         |           |          |       |           | ath:`e_3  |
      |         |           |          |       |           | = 0.223`, |
      |         |           |          |       |           | :ma       |
      |         |           |          |       |           | th:`\alph |
      |         |           |          |       |           | a_1 = 0`, |
      |         |           |          |       |           | :math:`   |
      |         |           |          |       |           | \alpha_2  |
      |         |           |          |       |           | = 0.436`, |
      |         |           |          |       |           | :math:`   |
      |         |           |          |       |           | \alpha_3  |
      |         |           |          |       |           | = 0.436`, |
      |         |           |          |       |           | :math:`   |
      |         |           |          |       |           | \gamma_1  |
      |         |           |          |       |           | = 0.436`, |
      |         |           |          |       |           | :math:`   |
      |         |           |          |       |           | \gamma_2  |
      |         |           |          |       |           | = 0.243`, |
      |         |           |          |       |           | :math:    |
      |         |           |          |       |           | `\gamma_3 |
      |         |           |          |       |           |  = 2.185` |
      +---------+-----------+----------+-------+-----------+-----------+
      | ROS–4   | 4         | 3        | 4(3)  | L-stable  | :m        |
      |         |           |          |       |           | ath:`a_{2 |
      |         |           |          |       |           | ,1} = 2`, |
      |         |           |          |       |           | :math:    |
      |         |           |          |       |           | `a_{3,1}  |
      |         |           |          |       |           | = 1.868`, |
      |         |           |          |       |           | :math:    |
      |         |           |          |       |           | `a_{3,2}  |
      |         |           |          |       |           | = 0.234`, |
      |         |           |          |       |           | :math:`a  |
      |         |           |          |       |           | _{4,1} =  |
      |         |           |          |       |           | a_{3,1}`, |
      |         |           |          |       |           | :math:`a  |
      |         |           |          |       |           | _{4,2} =  |
      |         |           |          |       |           | a_{3,2}`, |
      |         |           |          |       |           | :m        |
      |         |           |          |       |           | ath:`a_{4 |
      |         |           |          |       |           | ,3} = 0`, |
      |         |           |          |       |           | :math:`   |
      |         |           |          |       |           | c_{2,1} = |
      |         |           |          |       |           |  -7.137`, |
      |         |           |          |       |           | :math:    |
      |         |           |          |       |           | `c_{3,1}  |
      |         |           |          |       |           | = 2.581`, |
      |         |           |          |       |           | :math:    |
      |         |           |          |       |           | `c_{3,2}  |
      |         |           |          |       |           | = 0.652`, |
      |         |           |          |       |           | :math:`   |
      |         |           |          |       |           | c_{4,1} = |
      |         |           |          |       |           |  -2.137`, |
      |         |           |          |       |           | :math:`   |
      |         |           |          |       |           | c_{4,2} = |
      |         |           |          |       |           |  -0.321`, |
      |         |           |          |       |           | :math:`   |
      |         |           |          |       |           | c_{4,3} = |
      |         |           |          |       |           |  -0.695`, |
      |         |           |          |       |           | :m        |
      |         |           |          |       |           | ath:`m_1  |
      |         |           |          |       |           | = 2.256`, |
      |         |           |          |       |           | :m        |
      |         |           |          |       |           | ath:`m_2  |
      |         |           |          |       |           | = 0.287`, |
      |         |           |          |       |           | :m        |
      |         |           |          |       |           | ath:`m_3  |
      |         |           |          |       |           | = 0.435`, |
      |         |           |          |       |           | :m        |
      |         |           |          |       |           | ath:`m_4  |
      |         |           |          |       |           | = 1.094`, |
      |         |           |          |       |           | :ma       |
      |         |           |          |       |           | th:`e_1 = |
      |         |           |          |       |           |  -0.282`, |
      |         |           |          |       |           | :ma       |
      |         |           |          |       |           | th:`e_2 = |
      |         |           |          |       |           |  -0.073`, |
      |         |           |          |       |           | :ma       |
      |         |           |          |       |           | th:`e_3 = |
      |         |           |          |       |           |  -0.108`, |
      |         |           |          |       |           | :ma       |
      |         |           |          |       |           | th:`e_4 = |
      |         |           |          |       |           |  -1.093`, |
      |         |           |          |       |           | :ma       |
      |         |           |          |       |           | th:`\alph |
      |         |           |          |       |           | a_1 = 0`, |
      |         |           |          |       |           | :math:`   |
      |         |           |          |       |           | \alpha_2  |
      |         |           |          |       |           | = 1.146`, |
      |         |           |          |       |           | :math:`   |
      |         |           |          |       |           | \alpha_3  |
      |         |           |          |       |           | = 0.655`, |
      |         |           |          |       |           | :         |
      |         |           |          |       |           | math:`\al |
      |         |           |          |       |           | pha_4 = \ |
      |         |           |          |       |           | alpha_3`, |
      |         |           |          |       |           | :math:`   |
      |         |           |          |       |           | \gamma_1  |
      |         |           |          |       |           | = 0.573`, |
      |         |           |          |       |           | :math:`\  |
      |         |           |          |       |           | gamma_2 = |
      |         |           |          |       |           |  -1.769`, |
      |         |           |          |       |           | :math:`   |
      |         |           |          |       |           | \gamma_3  |
      |         |           |          |       |           | = 0.759`, |
      |         |           |          |       |           | :math:`   |
      |         |           |          |       |           | \gamma_4  |
      |         |           |          |       |           | = -0.104` |
      +---------+-----------+----------+-------+-----------+-----------+
      | RODAS–3 | 4         | 3        | 3(2)  | Stiffly   | :m        |
      |         |           |          |       |           | ath:`a_{2 |
      |         |           |          |       | accurate  | ,1} = 0`, |
      |         |           |          |       |           | :m        |
      |         |           |          |       |           | ath:`a_{3 |
      |         |           |          |       |           | ,1} = 2`, |
      |         |           |          |       |           | :m        |
      |         |           |          |       |           | ath:`a_{3 |
      |         |           |          |       |           | ,2} = 0`, |
      |         |           |          |       |           | :m        |
      |         |           |          |       |           | ath:`a_{4 |
      |         |           |          |       |           | ,1} = 2`, |
      |         |           |          |       |           | :m        |
      |         |           |          |       |           | ath:`a_{4 |
      |         |           |          |       |           | ,2} = 0`, |
      |         |           |          |       |           | :m        |
      |         |           |          |       |           | ath:`a_{4 |
      |         |           |          |       |           | ,3} = 1`, |
      |         |           |          |       |           | :m        |
      |         |           |          |       |           | ath:`c_{2 |
      |         |           |          |       |           | ,1} = 4`, |
      |         |           |          |       |           | :m        |
      |         |           |          |       |           | ath:`c_{3 |
      |         |           |          |       |           | ,1} = 1`, |
      |         |           |          |       |           | :ma       |
      |         |           |          |       |           | th:`c_{3, |
      |         |           |          |       |           | 2} = -1`, |
      |         |           |          |       |           | :m        |
      |         |           |          |       |           | ath:`c_{4 |
      |         |           |          |       |           | ,1} = 1`, |
      |         |           |          |       |           | :ma       |
      |         |           |          |       |           | th:`c_{4, |
      |         |           |          |       |           | 2} = -1`, |
      |         |           |          |       |           | :math     |
      |         |           |          |       |           | :`c_{4,3} |
      |         |           |          |       |           |  = -8/3`, |
      |         |           |          |       |           | :math:`   |
      |         |           |          |       |           | m_1 = 2`, |
      |         |           |          |       |           | :math:`   |
      |         |           |          |       |           | m_2 = 0`, |
      |         |           |          |       |           | :math:`   |
      |         |           |          |       |           | m_3 = 1`, |
      |         |           |          |       |           | :math:`   |
      |         |           |          |       |           | m_4 = 1`, |
      |         |           |          |       |           | :math:`   |
      |         |           |          |       |           | e_1 = 0`, |
      |         |           |          |       |           | :math:`   |
      |         |           |          |       |           | e_2 = 0`, |
      |         |           |          |       |           | :math:`   |
      |         |           |          |       |           | e_3 = 0`, |
      |         |           |          |       |           | :math:`   |
      |         |           |          |       |           | e_4 = 1`, |
      |         |           |          |       |           | :ma       |
      |         |           |          |       |           | th:`\alph |
      |         |           |          |       |           | a_1 = 0`, |
      |         |           |          |       |           | :ma       |
      |         |           |          |       |           | th:`\alph |
      |         |           |          |       |           | a_2 = 0`, |
      |         |           |          |       |           | :ma       |
      |         |           |          |       |           | th:`\alph |
      |         |           |          |       |           | a_3 = 1`, |
      |         |           |          |       |           | :ma       |
      |         |           |          |       |           | th:`\alph |
      |         |           |          |       |           | a_4 = 1`, |
      |         |           |          |       |           | :math     |
      |         |           |          |       |           | :`\gamma_ |
      |         |           |          |       |           | 1 = 0.5`, |
      |         |           |          |       |           | :math     |
      |         |           |          |       |           | :`\gamma_ |
      |         |           |          |       |           | 2 = 1.5`, |
      |         |           |          |       |           | :ma       |
      |         |           |          |       |           | th:`\gamm |
      |         |           |          |       |           | a_3 = 0`, |
      |         |           |          |       |           | :m        |
      |         |           |          |       |           | ath:`\gam |
      |         |           |          |       |           | ma_4 = 0` |
      +---------+-----------+----------+-------+-----------+-----------+
      | RODAS–4 | 6         | 5        | 4(3)  | Stiffly   | :ma       |
      |         |           |          |       |           | th:`\alph |
      |         |           |          |       | accurate  | a_1 = 0`, |
      |         |           |          |       |           | :math:`   |
      |         |           |          |       |           | \alpha_2  |
      |         |           |          |       |           | = 0.386`, |
      |         |           |          |       |           | :math:`   |
      |         |           |          |       |           | \alpha_3  |
      |         |           |          |       |           | = 0.210`, |
      |         |           |          |       |           | :math:`   |
      |         |           |          |       |           | \alpha_4  |
      |         |           |          |       |           | = 0.630`, |
      |         |           |          |       |           | :ma       |
      |         |           |          |       |           | th:`\alph |
      |         |           |          |       |           | a_5 = 1`, |
      |         |           |          |       |           | :ma       |
      |         |           |          |       |           | th:`\alph |
      |         |           |          |       |           | a_6 = 1`, |
      |         |           |          |       |           | :math:    |
      |         |           |          |       |           | `\gamma_1 |
      |         |           |          |       |           |  = 0.25`, |
      |         |           |          |       |           | :math:`\  |
      |         |           |          |       |           | gamma_2 = |
      |         |           |          |       |           |  -0.104`, |
      |         |           |          |       |           | :math:`   |
      |         |           |          |       |           | \gamma_3  |
      |         |           |          |       |           | = 0.104`, |
      |         |           |          |       |           | :math:`\  |
      |         |           |          |       |           | gamma_4 = |
      |         |           |          |       |           |  -0.036`, |
      |         |           |          |       |           | :ma       |
      |         |           |          |       |           | th:`\gamm |
      |         |           |          |       |           | a_5 = 0`, |
      |         |           |          |       |           | :ma       |
      |         |           |          |       |           | th:`\gamm |
      |         |           |          |       |           | a_6 = 0`, |
      |         |           |          |       |           | :math:    |
      |         |           |          |       |           | `a_{2,1}  |
      |         |           |          |       |           | = 1.544`, |
      |         |           |          |       |           | :math:    |
      |         |           |          |       |           | `a_{3,1}  |
      |         |           |          |       |           | = 0.946`, |
      |         |           |          |       |           | :math:    |
      |         |           |          |       |           | `a_{3,2}  |
      |         |           |          |       |           | = 0.255`, |
      |         |           |          |       |           | :math:    |
      |         |           |          |       |           | `a_{4,1}  |
      |         |           |          |       |           | = 3.314`, |
      |         |           |          |       |           | :math:    |
      |         |           |          |       |           | `a_{4,2}  |
      |         |           |          |       |           | = 2.896`, |
      |         |           |          |       |           | :math:    |
      |         |           |          |       |           | `a_{4,3}  |
      |         |           |          |       |           | = 0.998`, |
      |         |           |          |       |           | :math:    |
      |         |           |          |       |           | `a_{5,1}  |
      |         |           |          |       |           | = 1.221`, |
      |         |           |          |       |           | :math:    |
      |         |           |          |       |           | `a_{5,2}  |
      |         |           |          |       |           | = 6.019`, |
      |         |           |          |       |           | :math:`   |
      |         |           |          |       |           | a_{5,3} = |
      |         |           |          |       |           |  12.537`, |
      |         |           |          |       |           | :math:`   |
      |         |           |          |       |           | a_{5,4} = |
      |         |           |          |       |           |  -0.687`, |
      |         |           |          |       |           | :math:`a  |
      |         |           |          |       |           | _{6,1} =  |
      |         |           |          |       |           | a_{5,1}`, |
      |         |           |          |       |           | :math:`a  |
      |         |           |          |       |           | _{6,2} =  |
      |         |           |          |       |           | a_{5,2}`, |
      |         |           |          |       |           | :math:`a  |
      |         |           |          |       |           | _{6,3} =  |
      |         |           |          |       |           | a_{5,3}`, |
      |         |           |          |       |           | :math:`a  |
      |         |           |          |       |           | _{6,4} =  |
      |         |           |          |       |           | a_{5,4}`, |
      |         |           |          |       |           | :m        |
      |         |           |          |       |           | ath:`a_{6 |
      |         |           |          |       |           | ,5} = 1`, |
      |         |           |          |       |           | :math:`   |
      |         |           |          |       |           | c_{2,1} = |
      |         |           |          |       |           |  -5.668`, |
      |         |           |          |       |           | :math:`   |
      |         |           |          |       |           | c_{3,1} = |
      |         |           |          |       |           |  -2.430`, |
      |         |           |          |       |           | :math:`   |
      |         |           |          |       |           | c_{3,2} = |
      |         |           |          |       |           |  -0.206`, |
      |         |           |          |       |           | :math:`   |
      |         |           |          |       |           | c_{4,1} = |
      |         |           |          |       |           |  -0.107`, |
      |         |           |          |       |           | :math:`   |
      |         |           |          |       |           | c_{4,2} = |
      |         |           |          |       |           |  -9.594`, |
      |         |           |          |       |           | :math:`   |
      |         |           |          |       |           | c_{4,3} = |
      |         |           |          |       |           |  -20.47`, |
      |         |           |          |       |           | :math:    |
      |         |           |          |       |           | `c_{5,1}  |
      |         |           |          |       |           | = 7.496`, |
      |         |           |          |       |           | :math:`   |
      |         |           |          |       |           | c_{5,2} = |
      |         |           |          |       |           |  -0.124`, |
      |         |           |          |       |           | :mat      |
      |         |           |          |       |           | h:`c_{5,3 |
      |         |           |          |       |           | } = -34`, |
      |         |           |          |       |           | :math:`   |
      |         |           |          |       |           | c_{5,4} = |
      |         |           |          |       |           |  11.708`, |
      |         |           |          |       |           | :math:    |
      |         |           |          |       |           | `c_{6,1}  |
      |         |           |          |       |           | = 8.083`, |
      |         |           |          |       |           | :math:`   |
      |         |           |          |       |           | c_{6,2} = |
      |         |           |          |       |           |  -7.981`, |
      |         |           |          |       |           | :math:`c  |
      |         |           |          |       |           | _{6,3} =  |
      |         |           |          |       |           | -31.521`, |
      |         |           |          |       |           | :math:`   |
      |         |           |          |       |           | c_{6,4} = |
      |         |           |          |       |           |  16.319`, |
      |         |           |          |       |           | :math:`   |
      |         |           |          |       |           | c_{6,5} = |
      |         |           |          |       |           |  -6.058`, |
      |         |           |          |       |           | :mat      |
      |         |           |          |       |           | h:`m_1 =  |
      |         |           |          |       |           | a_{5,1}`, |
      |         |           |          |       |           | :mat      |
      |         |           |          |       |           | h:`m_2 =  |
      |         |           |          |       |           | a_{5,2}`, |
      |         |           |          |       |           | :mat      |
      |         |           |          |       |           | h:`m_3 =  |
      |         |           |          |       |           | a_{5,3}`, |
      |         |           |          |       |           | :mat      |
      |         |           |          |       |           | h:`m_4 =  |
      |         |           |          |       |           | a_{5,4}`, |
      |         |           |          |       |           | :math:`   |
      |         |           |          |       |           | m_5 = 1`, |
      |         |           |          |       |           | :math:`   |
      |         |           |          |       |           | m_6 = 1`, |
      |         |           |          |       |           | :math:`   |
      |         |           |          |       |           | e_1 = 0`, |
      |         |           |          |       |           | :math:`   |
      |         |           |          |       |           | e_2 = 0`, |
      |         |           |          |       |           | :math:`   |
      |         |           |          |       |           | e_3 = 0`, |
      |         |           |          |       |           | :math:`   |
      |         |           |          |       |           | e_4 = 0`, |
      |         |           |          |       |           | :math:`   |
      |         |           |          |       |           | e_5 = 0`, |
      |         |           |          |       |           | :math:    |
      |         |           |          |       |           | `e_6 = 1` |
      +---------+-----------+----------+-------+-----------+-----------+
      |         |           |          |       |           |           |
      +---------+-----------+----------+-------+-----------+-----------+

.. container:: table*

   .. container:: center

      +-------------+--------------------------+--------------------------+
      |             | File(s)                  | Description              |
      +-------------+--------------------------+--------------------------+
      | Runge-Kutta | ``runge_kutta.f90``      | Fully implicit 3-stage   |
      |             |                          | Runge-Kutta methods.     |
      |             |                          | Several variants are     |
      |             |                          | available:               |
      +-------------+--------------------------+--------------------------+
      |             |                          | RADAU-2A: order 5        |
      +-------------+--------------------------+--------------------------+
      |             |                          | RADAU-1A: order 5        |
      +-------------+--------------------------+--------------------------+
      |             |                          | Lobatto-3C: order 4      |
      +-------------+--------------------------+--------------------------+
      |             |                          | Gauss: order 6           |
      +-------------+--------------------------+--------------------------+
      | RADAU5      | ``atm_radau5.f``,        | This Runge-Kutta method  |
      |             | ``kpp_radau5.f90``       | of order 5 based on      |
      |             |                          | RADAU-IIA quadrature     |
      |             |                          | :raw-latex:`\citep` is   |
      |             |                          | stiffly accurate. The    |
      |             |                          | KPP implementation       |
      |             |                          | follows the original     |
      |             |                          | implementation of        |
      |             |                          | :raw-latex:`\citet`.     |
      |             |                          | While RADAU5 is          |
      |             |                          | relatively expensive     |
      |             |                          | (when compared to the    |
      |             |                          | Rosenbrock methods), it  |
      |             |                          | is more robust and is    |
      |             |                          | useful to obtain         |
      |             |                          | accurate reference       |
      |             |                          | solutions.               |
      +-------------+--------------------------+--------------------------+
      | SDIRK       | ``sdirk.f``,             | SDIRK is an L-stable,    |
      |             | ``sdirk.f90``            | si                       |
      |             |                          | ngly-diagonally-implicit |
      |             |                          | Runge-Kutta method. The  |
      |             |                          | implementation is based  |
      |             |                          | on :raw-latex:`\citet`.  |
      |             |                          | Several variants are     |
      |             |                          | available:               |
      +-------------+--------------------------+--------------------------+
      |             |                          | Sdirk 2a, 2b: 2 stages,  |
      |             |                          | order 2                  |
      +-------------+--------------------------+--------------------------+
      |             |                          | Sdirk 3a: 3 stages,      |
      |             |                          | order 2, and             |
      +-------------+--------------------------+--------------------------+
      |             |                          | Sdirk 4a, 4b: 5 stages,  |
      |             |                          | order 4                  |
      +-------------+--------------------------+--------------------------+
      | SDIRK4      | ``kpp_sdirk4.f``,        | SDIRK4 is an L-stable,   |
      |             | ``kpp_sdirk4.f90``       | si                       |
      |             |                          | ngly-diagonally-implicit |
      |             |                          | Runge-Kutta method of    |
      |             |                          | order 4. The             |
      |             |                          | implementation is based  |
      |             |                          | on :raw-latex:`\citet`.  |
      +-------------+--------------------------+--------------------------+
      | SEULEX      | ``kpp_seulex.f``,        | SEULEX is a variable     |
      |             | ``kpp_seulex.f90``       | order stiff              |
      |             |                          | extrapolation code able  |
      |             |                          | to produce highly        |
      |             |                          | accurate solutions. The  |
      |             |                          | KPP implementation is    |
      |             |                          | based on the             |
      |             |                          | implementation of        |
      |             |                          | :raw-latex:`\citet`.     |
      +-------------+--------------------------+--------------------------+
      |             |                          |                          |
      +-------------+--------------------------+--------------------------+

Tangent Linear Model
~~~~~~~~~~~~~~~~~~~~

The method (`[eqn:altRosenbrock] <#eqn:altRosenbrock>`__) is combined
with the sensitivity equations. One step of the method reads

.. math::

   \begin{aligned}
   \label{eqn:altRosenbrock-sen}
   %y^{n+1} &=& y^n + \sum_{i=1}^s m_i k_i, \qquad
   \delta y^{n+1} &=& \delta y^n + \sum_{i=1}^s m_i \ell_i\\
   \nonumber
   T_i &=& t^n + \alpha_i h~, %\quad Y_i =y^n + \sum_{j=1}^{i-1} a_{ij} k_j~,
   \quad \delta Y_i = \delta y^n + \sum_{j=1}^{i-1} a_{ij} \ell_j\\
   %A &=& \left[ \frac{1}{h \gamma} - J^T(t^n,y^n) \right]\\
   %\nonumber
   %A \cdot k_i &=&
   %           f\left( \, T_i,\, Y_i \,\right)
   %           + \sum_{j=1}^{i-1} \frac{c_{ij}}{h} k_j
   %          + h \gamma_i f_t\left(t^n,y^n\right)~,\\
   \nonumber
   A \cdot \ell_i &=&
           J\left( \, T_i,\, Y_i \,\right)
                 \cdot \delta Y_i
                 + \sum_{j=1}^{i-1} \frac{c_{ij}}{h} \ell_j\\
   \nonumber
   && +
   \left( H( t^n, y^n )\times  k_i \right) \cdot \delta y^n
      + h \gamma_i J_t\left(t^n,y^n\right) \cdot \delta y^n\end{aligned}

The method requires a single :math:`n \times n` LU decomposition per
step to obtain both the concentrations and the sensitivities.

KPP contains tangent linear models (for direct decoupled sensitivity
analysis) for each of the Rosenbrock methods (ROS–2, ROS–3, ROS–4,
RODAS–3, and RODAS–4). The implementations distinguish between
sensitivities with respect to initial values and sensitivities with
respect to parameters for efficiency.

The Discrete Adjoint
~~~~~~~~~~~~~~~~~~~~

To obtain the adjoint we first differentiate the method with respect to
:math:`y_n`. Here :math:`J` denotes the Jacobian and :math:`H` the
Hessian of the derivative function :math:`f`. The discrete adjoint of
the (non-autonomous) Rosenbrock method is

.. math::

   \begin{aligned}
   \label{Ros_disc_adj}
   %A &=& \left[ \frac{1}{h \gamma} - J^T(t^n,y^n) \right]\\
   %\nonumber
   A \cdot u_i
   &=& m_i \lambda^{n+1} + \sum_{j=i+1}^s \left( a_{ji} v_j + \frac{c_{ji}}{h}
   u_j \right)~,\\
   \nonumber
   v_i &=& J^T(T_i,Y_i)\cdot u_i~, \quad i = s,s-1,\cdots,1~,\\
   \nonumber
   \lambda^n &=& \lambda^{n+1} + \sum_{i=1}^s \left( H(t^n,y^n) \times
   k_i\right)^T
   \cdot u_i\\
   \nonumber
   && + h J^T_t(t^n,y^n) \cdot \sum_{i=1}^s \gamma_i u_i+  \sum_{i=1}^s v_i\end{aligned}

KPP contains adjoint models (for direct decoupled sensitivity analysis)
for each of the Rosenbrock methods (ROS–2, ROS–3, ROS–4, RODAS–3, and
RODAS–4).

Runge-Kutta methods
-------------------

A general :math:`s`-stage Runge-Kutta method is defined as
:raw-latex:`\cite[Section II.1]{k:HW1}`

.. math::

   \begin{aligned}
   \label{eqn:RungeKutta}
   y^{n+1} &=& y^n + h \sum_{i=1}^s b_i k_i~,\\
   \nonumber
   T_i &=& t^n + c_i h~, \quad
   Y_i = y^n + h \sum_{j=1}^{s} a_{ij} k_j~,\\
   \nonumber
   k_i &=& f\left( \, T_i, \, Y_i \,\right)~,\end{aligned}

where the coefficients :math:`a_{ij}`, :math:`b_i` and :math:`c_i` are
prescribed for the desired accuracy and stability properties. The stage
derivative values :math:`k_i` are defined implicitly, and require
solving a (set of) nonlinear system(s). Newton-type methods solve
coupled linear systems of dimension (at most) :math:`n \times s`.

The Runge-Kutta methods implemented in KPP are summarized in
Table `[tab:Runge-Kutta] <#tab:Runge-Kutta>`__.

.. _tangent-linear-model-1:

Tangent Linear Model
~~~~~~~~~~~~~~~~~~~~

The tangent linear method associated with the Runge-Kutta method is

.. math::

   \begin{aligned}
   \label{eqn:RK-TLM}
   %y^{n+1} &=& y^n + h \sum_{i=1}^s b_i k_i~,\\
   \delta y^{n+1} &=& \delta y^n + h \sum_{i=1}^s b_i \ell_i~,\\
   \nonumber
   %Y_i &=& y^n + h \sum_{j=1}^{s} a_{ij} k_j~,\\
   \delta Y_i& =& \delta y^n + h \sum_{j=1}^{s} a_{ij} \ell_j~,\\
   \nonumber
   %k_i &=& f\left( \, T_i, \, Y_i \,\right)~,\\
   \ell_i &=& J\left(T_i, \, Y_i \right) \cdot \delta Y_i ~.\end{aligned}

The system (`[eqn:RK-TLM] <#eqn:RK-TLM>`__) is linear and does not
require an iterative procedure. However, even for a SDIRK method
(:math:`a_{ij}=0` for :math:`i>j` and :math:`a_{ii}=\gamma`) each stage
requires the LU factorization of a different matrix.

Discrete Adjoint Model
~~~~~~~~~~~~~~~~~~~~~~

The first order Runge-Kutta adjoint is

.. math::

   \begin{aligned}
   \label{RK-adj}
   u_i &=& h \, J^T(T_i,Y_i)\cdot
   \left( b_i \lambda^{n+1} + \sum_{j=1}^s a_{ji} u_j \right)\\ %\quad i = 1 \cdots s\\
   \nonumber
   \lambda^{n} &=& \lambda^{n+1} +\sum_{j=1}^s u_j~.\end{aligned}

For :math:`b_i \ne 0` the Runge-Kutta adjoint can be rewritten as
another Runge-Kutta method:

.. math::

   \begin{aligned}
   \label{RK-adj-2}
   u_i &=& h \, J^T(T_i,Y_i)\cdot
   \left( \lambda^{n+1} + \sum_{j=1}^s \frac{b_j \,
   a_{ji}}{b_i} u_j \right)\\ %~, \quad i = 1 \cdots s\\
   \nonumber
   \lambda^{n} &=& \lambda^{n+1} +\sum_{j=1}^s b_j \, u_j~.\end{aligned}

Backward Differentiation Formulas
---------------------------------

Backward differentiation formulas (BDF) are linear multistep methods
with excellent stability properties for the integration of chemical
systems :raw-latex:`\citep[Section V.1]{k:HW2}`. The :math:`k`-step BDF
method reads

.. math::

   \sum_{i=0}^k \alpha_i y^{n-i} = h_n \beta\; f\left(t^{n},y^{n}\right)
   \label{BDF}

where the coefficients :math:`\alpha_i` and :math:`\beta` are chosen
such that the method has order of consistency :math:`k`.

The KPP library contains two off-the-shelf, highly popular
implementations of BDF methods, described in
Table `[tab:BDF] <#tab:BDF>`__.

.. container:: table*

   .. container:: center

      +--------+-------------------+---------------------------------------+
      |        | File(s)           | Description                           |
      +--------+-------------------+---------------------------------------+
      | LSODE  | ``kpp_lsode.f90`` | LSODE, the Livermore ODE solver       |
      |        |                   | :raw-latex:`\citep`, implements       |
      |        |                   | backward differentiation formula      |
      |        |                   | (BDF) methods for stiff problems.     |
      |        |                   | LSODE has been translated to          |
      |        |                   | Fortran90 for the incorporation into  |
      |        |                   | the KPP library.                      |
      +--------+-------------------+---------------------------------------+
      | LSODES | ``atm_lsodes.f``  | LSODES :raw-latex:`\citep`, the       |
      |        |                   | sparse version of the Livermore ODE   |
      |        |                   | solver LSODE, is modified to          |
      |        |                   | interface directly with the KPP       |
      |        |                   | generated code                        |
      +--------+-------------------+---------------------------------------+
      | VODE   | ``kpp_dvode.f``   | VODE :raw-latex:`\citep` uses another |
      |        |                   | formulation of backward               |
      |        |                   | differentiation formulas. The version |
      |        |                   | of VODE present in the KPP library    |
      |        |                   | uses directly the KPP sparse linear   |
      |        |                   | algebra routines.                     |
      +--------+-------------------+---------------------------------------+
      | ODESSA | ``atm_odessa.f``  | The BDF-based direct-decoupled        |
      |        |                   | sensitivity integrator Odessa         |
      |        |                   | :raw-latex:`\citep` has been modified |
      |        |                   | to use the KPP sparse linear algebra  |
      |        |                   | routines.                             |
      +--------+-------------------+---------------------------------------+
      |        |                   |                                       |
      +--------+-------------------+---------------------------------------+

Revision history
================

Only the major new features are listed here. For a detailed description
of the changes, read the file ``$KPP_HOME/CHANGELOG``.

KPP-2.5.0
---------

**((TODO: ADD TEXT HERE))**

KPP-2.4.0
---------

**((TODO: ADD TEXT HERE))**

KPP-2.2.3
---------

New features of KPP-2.2.3
~~~~~~~~~~~~~~~~~~~~~~~~~

-  A new function called ``k_3rd_iupac`` is available, calculating
   third-order rate coefficients using the formula used by IUPAC
   :raw-latex:`\citep{1610}`.

-  While previous versions of KPP were using ``yacc`` (yet another
   compiler compiler), the current version has been modified to be
   compatible with the parser generator ``bison``, which is the
   successor of yacc.

-  The new Runge-Kutta integrators ``kpp_sdirk4``, ``runge_kutta``, and
   ``sdirk`` were added.

-  The new KPP command ``#DECLARE`` was added (see
   Sect. `3.2.1 <#sec:command-declare>`__).

KPP-2.1
-------

New features of KPP-2.1
~~~~~~~~~~~~~~~~~~~~~~~

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
   `4.3 <#sec:f77>`__).

-  An automatically generated Makefile facilitates the compilation of
   the KPP-generated code (see Sect. `4.1.18 <#sec:output-makefile>`__).

-  Equation tags provide a convenient way to refer to specific chemical
   reactions (see Sect. `4.1.5 <#sec:output-monitor>`__).

-  The dummy index allows to test if a certain species occurs in the
   current chemistry mechanism. (see
   Sect. `3.2.4 <#sec:command-dummyindex>`__).

-  Lines starting with are comment lines.

Acknowledgements
================

Parts of this user manual are based on the thesis of
:raw-latex:`\citet{1693}`. We thank Jason Lander for his suggestions how
to migrate from to .

.. _`sec:bnf`:

BNF Description of the KPP Language
===================================

Following is the BNF-like specification of the language:

+-------------------+------------------------------------+-----------+
| *program* ::=     | *module* :math:`|` *module*        |           |
|                   | *program*                          |           |
+-------------------+------------------------------------+-----------+
| *module* ::=      | *section* :math:`|` *command*      |           |
|                   | :math:`|` *inline_code*            |           |
+-------------------+------------------------------------+-----------+
| *section* ::=     | *atom_definition_list*             | :math:`|` |
+-------------------+------------------------------------+-----------+
|                   | *atom_list*                        | :math:`|` |
+-------------------+------------------------------------+-----------+
|                   | *species_definition_list*          | :math:`|` |
+-------------------+------------------------------------+-----------+
|                   | *species_definition_list*          | :math:`|` |
+-------------------+------------------------------------+-----------+
|                   | *equation_list*                    | :math:`|` |
+-------------------+------------------------------------+-----------+
|                   | *initvalues_list*                  | :math:`|` |
+-------------------+------------------------------------+-----------+
|                   | *species_list* *atom_list*         | :math:`|` |
+-------------------+------------------------------------+-----------+
|                   | *lump_list*                        | :math:`|` |
+-------------------+------------------------------------+-----------+
|                   | *species_list* *atom_list*         | :math:`|` |
+-------------------+------------------------------------+-----------+
|                   | *species_list_plus*                | :math:`|` |
+-------------------+------------------------------------+-----------+
|                   | *species_list_plus*                | :math:`|` |
+-------------------+------------------------------------+-----------+
|                   | *species_list*                     |           |
+-------------------+------------------------------------+-----------+
| *command* ::=     |                                    | :math:`|` |
+-------------------+------------------------------------+-----------+
|                   | ``[ SYMBOL | VALUE ]``             | :math:`|` |
+-------------------+------------------------------------+-----------+
|                   | ``[ ON | OFF ]``                   | :math:`|` |
+-------------------+------------------------------------+-----------+
|                   | *driver_name*                      | :math:`|` |
+-------------------+------------------------------------+-----------+
|                   | ``[ ON | OFF ]``                   | :math:`|` |
+-------------------+------------------------------------+-----------+
|                   | ``[ ON | OFF ]``                   | :math:`|` |
+-------------------+------------------------------------+-----------+
|                   | ``[ AGGREGATE | SPLIT ]``          | :math:`|` |
+-------------------+------------------------------------+-----------+
|                   | ``[ ON | OFF ]``                   | :math:`|` |
+-------------------+------------------------------------+-----------+
|                   | *file_name*                        | :math:`|` |
+-------------------+------------------------------------+-----------+
|                   | *integrator_name*                  | :math:`|` |
+-------------------+------------------------------------+-----------+
|                   | *integrator_name*                  | :math:`|` |
+-------------------+------------------------------------+-----------+
|                   | ``[ OFF | FUL                      | :math:`|` |
|                   | L | SPARSE_LU_ROW | SPARSE_ROW ]`` |           |
+-------------------+------------------------------------+-----------+
|                   | ``[ Fort                           | :math:`|` |
|                   | ran90 | Fortran77 | C | Matlab ]`` |           |
+-------------------+------------------------------------+-----------+
|                   |                                    | :math:`|` |
+-------------------+------------------------------------+-----------+
|                   | ``[ ON | OFF ]``                   | :math:`|` |
+-------------------+------------------------------------+-----------+
|                   | *minimum_version_number*           | :math:`|` |
+-------------------+------------------------------------+-----------+
|                   | *model_name*                       | :math:`|` |
+-------------------+------------------------------------+-----------+
|                   | ``[ ON | OFF ]``                   | :math:`|` |
+-------------------+------------------------------------+-----------+
|                   | ``[ ON | OFF ]``                   | :math:`|` |
+-------------------+------------------------------------+-----------+
|                   | ``[ ON | OFF ]``                   | :math:`|` |
+-------------------+------------------------------------+-----------+
|                   |                                    | :math:`|` |
+-------------------+------------------------------------+-----------+
|                   | ``[ ON | OFF ]``                   |           |
+-------------------+------------------------------------+-----------+
| *inline_code* ::= | *inline_type*                      |           |
+-------------------+------------------------------------+-----------+
|                   | *inline_code*                      |           |
+-------------------+------------------------------------+-----------+
|                   |                                    |           |
+-------------------+------------------------------------+-----------+

+---------------------------+---------------------------+-----------+
| *atom_count* ::=          | *integer* *atom_name*     | :math:`|` |
+---------------------------+---------------------------+-----------+
|                           | *atom_name*               |           |
+---------------------------+---------------------------+-----------+
| *atom_definition_list*    | *atom_definition*         | :math:`|` |
| ::=                       |                           |           |
+---------------------------+---------------------------+-----------+
|                           | *atom_definition*         |           |
|                           | *atom_definition_list*    |           |
+---------------------------+---------------------------+-----------+
| *atom_list* ::=           | *atom_name*;              | :math:`|` |
+---------------------------+---------------------------+-----------+
|                           | *atom_name*; *atom_list*  |           |
+---------------------------+---------------------------+-----------+
| *equation* ::=            | :math:`<`\                | :math:`|` |
|                           | *equation_tag*\ :math:`>` |           |
|                           | *expression* =            |           |
|                           | *expression* : *rate*;    |           |
+---------------------------+---------------------------+-----------+
|                           | *expression* =            |           |
|                           | *expression* : *rate*;    |           |
+---------------------------+---------------------------+-----------+
| *equation_list* ::=       | *equation*                | :math:`|` |
+---------------------------+---------------------------+-----------+
|                           | *equation*                |           |
|                           | *equation_list*           |           |
+---------------------------+---------------------------+-----------+
| *equation_tag* ::=        | Alphanumeric expression,  |           |
|                           | also including the        |           |
|                           | underscore.               |           |
+---------------------------+---------------------------+-----------+
|                           | In , it is defined as .   |           |
+---------------------------+---------------------------+-----------+
| *expression* ::=          | *term*                    | :math:`|` |
+---------------------------+---------------------------+-----------+
|                           | *term* :math:`+`          | :math:`|` |
|                           | *expression*              |           |
+---------------------------+---------------------------+-----------+
|                           | *term* :math:`-`          |           |
|                           | *expression*              |           |
+---------------------------+---------------------------+-----------+
| *initvalues_assignment*   | *species_name_plus* =     | :math:`|` |
| ::=                       | *program_expression*;     |           |
+---------------------------+---------------------------+-----------+
|                           | = *program_expression*;   |           |
+---------------------------+---------------------------+-----------+
| *initvalues_list* ::=     | *initvalues_assignment*   | :math:`|` |
+---------------------------+---------------------------+-----------+
|                           | *initvalues_assignment*   |           |
|                           | *initvalues_list*         |           |
+---------------------------+---------------------------+-----------+
| *inline_type* ::=         | :math:`|` :math:`|`       | :math:`|` |
|                           | :math:`|` :math:`|`       |           |
|                           | :math:`|`                 |           |
+---------------------------+---------------------------+-----------+
|                           | :math:`|` :math:`|`       | :math:`|` |
|                           | :math:`|` :math:`|`       |           |
|                           | :math:`|`                 |           |
+---------------------------+---------------------------+-----------+
|                           | :math:`|` :math:`|`       | :math:`|` |
|                           | :math:`|` :math:`|`       |           |
|                           | :math:`|`                 |           |
+---------------------------+---------------------------+-----------+
|                           | :math:`|` :math:`|`       | :math:`|` |
+---------------------------+---------------------------+-----------+
|                           | :math:`|` :math:`|`       |           |
+---------------------------+---------------------------+-----------+
| *lump* ::=                | *lump_sum* :              |           |
|                           | *species_name*;           |           |
+---------------------------+---------------------------+-----------+
| *lump_list* ::=           | *lump*                    | :math:`|` |
+---------------------------+---------------------------+-----------+
|                           | *lump* *lump_list*        |           |
+---------------------------+---------------------------+-----------+
| *lump_sum* ::=            | *species_name*            | :math:`|` |
+---------------------------+---------------------------+-----------+
|                           | *species_name* :math:`+`  |           |
|                           | *lump_sum*                |           |
+---------------------------+---------------------------+-----------+
| *rate* ::=                | *number*                  | :math:`|` |
+---------------------------+---------------------------+-----------+
|                           | *program_expression*      |           |
+---------------------------+---------------------------+-----------+
| *species_composition* ::= | *atom_count*              | :math:`|` |
+---------------------------+---------------------------+-----------+
|                           | *atom_count* :math:`+`    | :math:`|` |
|                           | *species_composition*     |           |
+---------------------------+---------------------------+-----------+
|                           |                           |           |
+---------------------------+---------------------------+-----------+
| *species_definition* ::=  | *species_name* =          |           |
|                           | *species_composition*;    |           |
+---------------------------+---------------------------+-----------+
| *species_definition_list* | *species_definition*      | :math:`|` |
| ::=                       |                           |           |
+---------------------------+---------------------------+-----------+
|                           | *species_definition*      |           |
|                           | *species_definition_list* |           |
+---------------------------+---------------------------+-----------+
| *species_list* ::=        | *species_name*;           | :math:`|` |
+---------------------------+---------------------------+-----------+
|                           | *species_name*;           |           |
|                           | *species_list*            |           |
+---------------------------+---------------------------+-----------+
| *species_list_plus* ::=   | *species_name_plus*;      | :math:`|` |
+---------------------------+---------------------------+-----------+
|                           | *species_name_plus*;      |           |
|                           | *species_list_plus*       |           |
+---------------------------+---------------------------+-----------+
| *species_name* ::=        | Alphanumeric expression,  |           |
|                           | also including the        |           |
|                           | underscore, starting with |           |
|                           | a letter.                 |           |
+---------------------------+---------------------------+-----------+
|                           | In , it is defined as .   |           |
|                           | Its maximum length is 15. |           |
+---------------------------+---------------------------+-----------+
| *species_name_plus* ::=   | *species_name*            | :math:`|` |
+---------------------------+---------------------------+-----------+
|                           |                           | :math:`|` |
+---------------------------+---------------------------+-----------+
|                           |                           | :math:`|` |
+---------------------------+---------------------------+-----------+
|                           |                           |           |
+---------------------------+---------------------------+-----------+
| *term* ::=                | *number* *species_name*   | :math:`|` |
+---------------------------+---------------------------+-----------+
|                           | *species_name*            | :math:`|` |
+---------------------------+---------------------------+-----------+
|                           |                           | :math:`|` |
+---------------------------+---------------------------+-----------+
|                           |                           |           |
+---------------------------+---------------------------+-----------+

.. |image| image:: Figures/kpp-logo.pdf
   :width: 70.0%
