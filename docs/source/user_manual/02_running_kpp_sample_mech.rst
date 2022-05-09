.. _running-kpp-with-an-example-mechanism:

###################################################
Running KPP with an Example Stratospheric Mechanism
###################################################

Here we consider as an example a very simple Chapman-like mechanism for
stratospheric chemistry:

.. math::

   \begin{aligned}
   O_2    + h\nu   & \longrightarrow  & 2 O           & ~~~~~~~~~~ (1)\\
   O      + O_2    & \longrightarrow  & O_3           & ~~~~~~~~~~ (2)\\
   O_3    + h\nu   & \longrightarrow  & O      + O_2  & ~~~~~~~~~~ (3)\\
   O      + O_3    & \longrightarrow  & 2 O_2         & ~~~~~~~~~~ (4)\\
   O_3    + h\nu   & \longrightarrow  & O(^1D) + O_2  & ~~~~~~~~~~ (5)\\
   O(^1D) + M      & \longrightarrow  & O + M         & ~~~~~~~~~~ (6)\\
   O(^1D) + O_3    & \longrightarrow  & 2 O_2         & ~~~~~~~~~~ (7)\\
   NO     + O_3    & \longrightarrow  & NO_2   + O_2  & ~~~~~~~~~~ (8)\\
   NO_2   + O      & \longrightarrow  & NO     + O_2  & ~~~~~~~~~~ (9)\\
   NO_2   + h\nu   & \longrightarrow  & NO     + O    & ~~~~~~~~~~ (10)
   \end{aligned}

We use the mechanism with the purpose of illustrating the KPP
capabilities. However, the software tools are general and can be applied
to virtually any kinetic mechanism.

We focus on Fortran90. Particularities of the C, Fortran77, and Matlab
languages are discussed in `Target language selection
<target-language-selection_>`_. 

The KPP input files (with suffix :file:`.kpp`) specify the model, the
target language, the precision, the integrator and the driver, etc. The file
name (without the suffix :file:`.kpp`) serves as the root name for the
simulation. In this paper we will refer to this name as root. Since
the root name will be incorporated into Fortran90 module names, it can
only contain valid Fortran90 characters, i.e. letters, numbers, and the
underscore. To specify a KPP model, write a root file with the following lines: 

.. code-block:: console

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

.. code-block:: console

   #LANGUAGE Fortran90

Here, we have chosen Fortran90. See `Target language selection
<target-language-selection_>`_ for other options. 

The data type of the generated model can be switched between
single/double precision with the command . The command selects a
specific numerical integration routine (from the templates provided by
KPP or implemented by the user) and the command selects a specific main
program. The command selects a specific kinetic mechanism. In our
example the model definition file includes the species and the equation
files,

.. code-block:: console

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

.. code-block:: console

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

.. code-block:: console

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

.. code-block:: console

   $ kpp small_strato.kpp

Next, compile and run the Fortran90 code:

.. code-block:: console

   $ make -f Makefile_small_strato
   $ ./small_strato.exe
