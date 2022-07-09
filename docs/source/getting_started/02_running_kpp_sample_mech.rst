.. _running-kpp-with-an-example-mechanism:

###################################################
Running KPP with an example stratospheric mechanism
###################################################

Here we consider as an example a very simple Chapman-like mechanism for
stratospheric chemistry:

.. math::

   \begin{aligned}
   O_2    + h\nu   & \longrightarrow  & 2 O           & ~~~~~~~~~~ (R1)\\
   O      + O_2    & \longrightarrow  & O_3           & ~~~~~~~~~~ (R2)\\
   O_3    + h\nu   & \longrightarrow  & O      + O_2  & ~~~~~~~~~~ (R3)\\
   O      + O_3    & \longrightarrow  & 2 O_2         & ~~~~~~~~~~ (R4)\\
   O_3    + h\nu   & \longrightarrow  & O(^1D) + O_2  & ~~~~~~~~~~ (R5)\\
   O(^1D) + M      & \longrightarrow  & O + M         & ~~~~~~~~~~ (R6)\\
   O(^1D) + O_3    & \longrightarrow  & 2 O_2         & ~~~~~~~~~~ (R7)\\
   NO     + O_3    & \longrightarrow  & NO_2   + O_2  & ~~~~~~~~~~ (R8)\\
   NO_2   + O      & \longrightarrow  & NO     + O_2  & ~~~~~~~~~~ (R9)\\
   NO_2   + h\nu   & \longrightarrow  & NO     + O    & ~~~~~~~~~~ (R10)
   \end{aligned}

We use the mechanism with the purpose of illustrating the KPP
capabilities. However, the software tools are general and can be applied
to virtually any kinetic mechanism.

We focus on Fortran90. Particularities of the C and Matlab
languages are discussed in the :ref:`language-cmd` section.

.. important::

   Most of the recent KPP developments described in this manual have
   been added into the Fortran90 language.  We look to members of the
   KPP user community to spearhead development in C, Matlab, and other
   languages.

The KPP input files (with suffix :file:`.kpp`) specify the
:ref:`target model <model-cmd>`, the :ref:`target language
<language-cmd>`, the :ref:`integrator <integrator-cmd>` the
:ref:`driver program <driver-cmd>`. etc. The file name (without the
:file:`.kpp`) serves as the root name for the simulation. Here we will
refer to this name as :code:`ROOT`.  Since the root name will  be
incorporated into Fortran90 module names, it can only contain valid
Fortran90 characters, i.e. letters, numbers, and the underscore.

The sections below outline the steps necessary to build and run a
"box-model" simulation with an example mechanism.

.. _example-step-1:

=====================================
1. Create a directory for the example
=====================================

Create a directory in which to build and run the example mechanism:

.. code-block:: console

   $ cd $HOME
   $ mkdir small_strato_example
   $ cd small_strato_example

In the following sections we will refer to
:file:`$HOME/small_strato_example` as "the example directory".

.. _example-step-2:

===============================
2. Create a KPP Definition File
===============================

Create a KPP definition file in the example directory.  The name
of this file will always be :file:`ROOT.kpp`, where :file:`ROOT` is
the name of the chemical mechanism.

For this example, write the following lines into a file named
:file:`small_strato.kpp` in the example directory:

.. code-block:: console

   #MODEL      small_strato
   #LANGUAGE   Fortran90
   #INTEGRATOR rosenbrock
   #DRIVER     general

.. important::

   KPP will look for the relevant files (e.g. mechanism definition,
   driver, etc.) in the proper subdirectories of :envvar:`KPP_HOME`.
   Therefore you won't need to copy these manually to the example
   directory.

We will now look at the :ref:`kpp-commands` in :file:`small_strato.kpp`.

.. _example-model-ss:

#MODEL small_strato
-------------------

The :ref:`model-cmd` command selects a specific kinetic mechanism (in
this example, :program:`small_strato`). KPP will look in the path
:file:`$KPP_HOME/models/` for the *model definition file*
:file:`small_strato.def` which contains the following code in the
:ref:`KPP language <bnf-description>`:

.. code-block:: console

   #include small_strato.spc       { Includes file w/ species definitons     }
   #include small_strato.eqn       { Includes file w/ chemical equations     }

   #LOOKATALL                      { Output all species to small_strato.dat}
   #MONITOR O3;N;O2;O;NO;O1D;NO2;  { Print selected species to screen        }

   #CHECK O; N;                    { Check Mass Balance of oxygen & nitrogen }

   #INITVALUES                     { Set initial values of species           }
     CFACTOR = 1.    ;             { and set units conversion factor to 1    }
     O1D = 9.906E+01 ;
     O   = 6.624E+08 ;
     O3  = 5.326E+11 ;
     O2  = 1.697E+16 ;
     NO  = 8.725E+08 ;
     NO2 = 2.240E+08 ;
     M   = 8.120E+16 ;

   { Fortran code to be inlined into ROOT_Global }
   #INLINE F90_INIT
     TSTART = (12*3600)
     TEND = TSTART + (3*24*3600)
     DT = 0.25*3600
     TEMP = 270
   #ENDINLINE

   { Matlab code to be inlined into ROOT_Global }
   #INLINE MATLAB_INIT
     global TSTART TEND DT TEMP
     TSTART = (12*3600);
     TEND = TSTART + (3*24*3600);
     DT = 0.25*3600;
     TEMP = 270;
   #ENDINLINE

   { C code to be inlined into ROOT_GLOBAL }
   #INLINE C_INIT
     TSTART = (12*3600);
     TEND = TSTART + (3*24*3600);
     DT = 0.25*3600;
     TEMP = 270;
   #ENDINLINE

The *definition file* :file:`small_strato.def` uses the
:ref:`include-cmd` command to include the *species file* and the
*equation file*. It also specifies parameters for running a "box-model"
simulation, such as :ref:`species initial values <initvalues>`, start
time, stop, time, and timestep (cf. :ref:`inlined-code`).

The *species file* :file:`small_strato.spc` lists all the species in the
model. Some of them are variable, meaning that their concentrations
change according to the law of mass action kinetics. Others are fixed,
with the concentrations determined by physical and not chemical factors
(cf. :ref:`defvar-and-deffix`). For each species its atomic composition
is given (unless the user chooses to ignore it).

.. code-block:: console

   #INCLUDE atoms.kpp
   #DEFVAR
     O   = O;
     O1D = O;
     O3  = O + O + O;
     NO  = N + O;
     NO2 = N + O + O;
   #DEFFIX
     M   = IGNORE;
     O2  = O + O;

The species file also includes the *atoms file* (:file:`atoms.kpp`),
which defines the chemical elements in the :ref:`atoms` section.

The *equation file* :file:`small_strato.eqn` contains the description of
the equations in an :ref:`equations` section. The chemical kinetic
mechanism is specified in the :ref:`KPP language <bnf-description>`.
Each reaction is described as “the sum of reactants equals the sum of
products” and, after a colon, is followed by its rate coefficient.
:code:`SUN` is the normalized sunlight intensity, equal to one at noon
and zero at night. Equation tags, e.g. :code:`<R1>`, are optional.

.. code-block:: console

   #EQUATIONS { Small Stratospheric Mechanism }


   <R1>  O2   + hv = 2O            : (2.643E-10) * SUN*SUN*SUN;
   <R2>  O    + O2 = O3            : (8.018E-17);
   <R3>  O3   + hv = O   + O2      : (6.120E-04) * SUN;
   <R4>  O    + O3 = 2O2           : (1.576E-15);
   <R5>  O3   + hv = O1D + O2      : (1.070E-03) * SUN*SUN;
   <R6>  O1D  + M  = O   + M       : (7.110E-11);
   <R7>  O1D  + O3 = 2O2           : (1.200E-10);
   <R8>  NO   + O3 = NO2 + O2      : (6.062E-15);
   <R9>  NO2  + O  = NO  + O2      : (1.069E-11);
   <R10> NO2  + hv = NO  + O       : (1.289E-02) * SUN;

.. _example-language-f90:

#LANGUAGE Fortran90
-------------------

The :ref:`language-cmd` command selects the language for the
KPP-generated solver code.  In this example we are using Fortran90.

.. _example-double-on:

#INTEGRATOR rosenbrock
----------------------

The :ref:`integrator-cmd` command selects a numerical integration routine
from the templates provided in the :file:`$KPP_HOME/int` directory, or
implemented by the user.

In this example, the :ref:`Rosenbrock integrator <rosenbrock-methods>`
and the Fortran90 language have been been specified. Therefore, the file
:file:`$KPP_HOME/int/rosenbrock.f90` will be used.
      
.. _example-driver-general:


#DRIVER general
---------------

The :ref:`driver-cmd` command selects a specific main program (located
in the :file:`$KPP_HOME/drv` directory):

#. :file:`general_adj.f90` : Used with integrators that use the
   discrete adjoint method
#. :file:`general_tlm.f90` : Used with integrators that use the
   tangent-linear method
#. :file:`general.f90` : Used with all other integrators.

In this example, the :file:`rosenbrock.f90` integrator does not use
either adjoint or tangent-linear methods, so the
:file:`$KPP_HOME/drv/general.f90` will be used.


.. _example-step-3:

===============================
3. Build the mechanism with KPP
===============================

Now that all the necessary files have been copied to the example directory,
the :program:`small_strato` mechanism can be built. Type:

.. code-block:: console

   $ kpp small_strato.kpp

You should see output similar to:

.. code-block:: console

   This is KPP-3.0.0.

   KPP is parsing the equation file.
   KPP is computing Jacobian sparsity structure.
   KPP is starting the code generation.
   KPP is initializing the code generation.
   KPP is generating the monitor data:
       - small_strato_Monitor
   KPP is generating the utility data:
       - small_strato_Util
   KPP is generating the global declarations:
       - small_strato_Main
   KPP is generating the ODE function:
       - small_strato_Function
   KPP is generating the ODE Jacobian:
       - small_strato_Jacobian
       - small_strato_JacobianSP
   KPP is generating the linear algebra routines:
       - small_strato_LinearAlgebra
   KPP is generating the Hessian:
       - small_strato_Hessian
       - small_strato_HessianSP
   KPP is generating the utility functions:
       - small_strato_Util
   KPP is generating the rate laws:
       - small_strato_Rates
   KPP is generating the parameters:
       - small_strato_Parameters
   KPP is generating the global data:
       - small_strato_Global
   KPP is generating the stoichiometric description files:
       - small_strato_Stoichiom
       - small_strato_StoichiomSP
   KPP is generating the driver from general.f90:
       - small_strato_Main
   KPP is starting the code post-processing.

   KPP has succesfully created the model "small_strato".

This will generate the Fortran90 code needed to solve the
:program:`small_strato` mechanism. The file listing should be similar
to:

.. code-block:: console

   atoms.kpp                     small_strato.kpp
   general.f90                   small_strato_LinearAlgebra.f90
   Makefile_small_strato         small_strato_Main.f90
   rosenbrock.def                small_strato_mex_Fun.f90
   rosenbrock.f90                small_strato_mex_Hessian.f90
   small_strato.def              small_strato_mex_Jac_SP.f90
   small_strato.eqn              small_strato_Model.f90
   small_strato_Function.f90     small_strato_Monitor.f90
   small_strato_Global.f90       small_strato_Parameters.f90
   small_strato_Hessian.f90      small_strato_Precision.f90
   small_strato_HessianSP.f90    small_strato_Rates.f90
   small_strato_Initialize.f90   small_strato.spc@
   small_strato_Integrator.f90   small_strato_Stoichiom.f90
   small_strato_Jacobian.f90     small_strato_StoichiomSP.f90
   small_strato_JacobianSP.f90   small_strato_Util.f90

KPP creates Fortran90 beginning with the mechanism name (which is
:file:`ROOT_` = :file:`small_strato_` in this example). KPP also
generates a human-readable summary of the mechanism
(:file:`small_strato.log`) as well as the Makefile
:file:`Makefile_small_strato` that can be used to build the executable.

.. _example_step_4:

=========================================
4. Build and run the small_strato example
=========================================

To compile the Fortran90 code generated by KPP into an executable, type:

.. code-block:: console

   $ make -f Makefile_small_strato

You will see output similar to this:

.. code-block:: console

   gfortran -cpp -O -g  -c small_strato_Precision
   gfortran -cpp -O -g  -c small_strato_Precision.f90
   gfortran -cpp -O -g  -c small_strato_Parameters.f90
   gfortran -cpp -O -g  -c small_strato_Global.f90
   gfortran -cpp -O -g  -c small_strato_Function.f90
   gfortran -cpp -O -g  -c small_strato_JacobianSP.f90
   gfortran -cpp -O -g  -c small_strato_Jacobian.f90
   gfortran -cpp -O -g  -c small_strato_HessianSP.f90
   gfortran -cpp -O -g  -c small_strato_Hessian.f90
   gfortran -cpp -O -g  -c small_strato_StoichiomSP.f90
   gfortran -cpp -O -g  -c small_strato_Stoichiom.f90
   gfortran -cpp -O -g  -c small_strato_Rates.f90
   gfortran -cpp -O -g  -c small_strato_Monitor.f90
   gfortran -cpp -O -g  -c small_strato_Util.f90
   gfortran -cpp -O -g  -c small_strato_LinearAlgebra.f90
   gfortran -cpp -O -g  -c small_strato_Initialize.f90
   gfortran -cpp -O -g  -c small_strato_Integrator.f90
   gfortran -cpp -O -g  -c small_strato_Model.f90
   gfortran -cpp -O -g  -c small_strato_Main.f90
   gfortran -cpp -O -g  small_strato_Precision.o    small_strato_Parameters.o    small_strato_Global.o small_strato_Function.o small_strato_JacobianSP.o small_strato_Jacobian.o small_strato_HessianSP.o small_strato_Hessian.o small_strato_Stoichiom.o small_strato_StoichiomSP.o small_strato_Rates.o   small_strato_Util.o   small_strato_Monitor.o small_strato_LinearAlgebra.o small_strato_Main.o          small_strato_Initialize.o small_strato_Integrator.o    small_strato_Model.o  -o small_strato.exe

Once compilation has finished, you can run the :program:`small_strato`
example by typing:

.. code-block:: console

   $ ./small_strato.exe | tee small_strato.log

This will run a "box-model" simulation forward several steps in time.
You will see the concentrations of selected species at several
timesteps displayed to the screen (aka the Unix stdout stream) as well
as to a log file (:file:`small_strato.log`).

If your simulation results exits abruptly with the :code:`Killed` error,
you probably need to increase your stack memory limit. On most Linux
systems the default stacksize limit is 8 kIb = or 8192 kB. You can max
this out with the following commands:

If you are using bash, type:

.. code-block:: console

   $ ulimit -s unlimited

If you are using csh, type:

.. code-block:: console

   $ limit stacksize unlimited

.. _example-step-5:

==========
5. Cleanup
==========

If you wish to remove the executable (:file:`small_strato.exe`), as
well as the object (:file:`*.o`) and module (:file:`*.mod`)
files generated by the Fortran compiler, type:

.. code-block:: console

   $ make -f Makefile_small_strato clean

If you also wish to remove all the files that were generated by KPP
(i.e. :file:`small_strato.log` and :file:`small_strato_*.f90`), type:

.. code-block:: console

   $ make -f Makefile_small_strato distclean
